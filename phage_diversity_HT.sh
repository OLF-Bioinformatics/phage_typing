#!/bin/bash

###################################
#                                 #
#    Prophage Sequence Typing     #
#                                 #
###################################


author='duceppemo'
version='0.2.0'


##############
#            #
#    Note    #
#            #
##############


# The pipleline will take as sample name everything that is before the first "_"
# in the file name. Thus, it requires that the "_" is present in your file name,
# which is usually the case with Illumina reads.


#################
#               #
#    Options    #
#               #
#################


function displayHelp()
{
    echo -e "\

    Usage:
        pahge_diversity_HT.sh -{1|2|3|4|5} -i /input/folder/ -o /output/folder/ \\
            [-n "$(nproc)"] [-p 2] [-m 8] [-k 21,33,55,77,99,127] -c 0.99,0.99,11

    Subcommands:
        1           Assemble
                    Trim and meege Illimina paired-end reads using BBtools
                    Assemble with SPAdes

        2           Submit
                    Submit samples to PHASTER web server

        3           Check
                    Check if results are available from the PHASTER server.
                    Useful when many samples are submitted or there is a large queue on the server
                    Allows to stop the pipeline and try to get the results later

        4           Cluster
                    Run the clustering (CD-HIT-EST and QIIME)
                    Useful to test multiple clustering settings

        5           All
                    Run assemble, submit, check and cluster

    Note:
        Subcommands should be run in order because they take as input the output of the previous one

    Options:
        -i          Input folder
                    Requires gziped Illumina paired-end files (e.g. "sample_R1.fastq.gz" and "sample_R2.fastq.gz")
                    Mandatory for \"all\" or \"assemble\" subcommands

        -o          Output folder
                    Madatory

        -n          Number of CPUs to use for parallel download
                    Default is max numbers of CPUs avaiable ("$(nproc)")
                    Optional

        -p          Parallelism level
                    Numbe of samples to run in parallel
                    Number of CPUs per sample = \"-n\" value / \"-p\" value
                    Default is 2
                    Optional

        -m          Memory in GB
                    Default is 85% of available memory ("${maxMem}"GB)
                    Optional

        -k          kmer value(s) for SPAdes assembler
                    A comma-separated list of k-mer sizes may be used
                    All values must be odd, less than 128 and listed in ascending order
                    Default is \"21,33,55,77,99,127\"
                    Optional

        -c          CD-HIT-EST settings
                    Format is \"sequence identity threshold,length difference cutoff, word size\"
                    For details, see https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#CDHITEST 
                    Default is \"0.99,0.99,11\"
                    Optional

        -h          Print this help message
    "
}

#Colored error message
BLUE='\033[1;34m'
NC='\033[0m' # No Color

# Default values
baseDir=''
reads=''
subcommand=''
export maxProc=2  # 2 samples parallel
export cpu=$(nproc)  # all available threads
maxMem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
export mem="$maxMem"
export kmer="21,33,55,77,99,127"  # k-mer sizes for SPAdes assembler (must be odd number(s))
clust="0.99,0.99,11"
help=''

# If the very first character of the option-string is a : (colon),
# which would normally be nonsense because there's no option letter preceding it,
# getopts switches to "silent error reporting mode"
# When you want getopts to expect an argument for an option, just place a : (colon) after the proper option flag
options=':12345i:o:m:n:p:k:c:h'

while getopts "$options" opt; do
    case "$opt" in
        1)          subcommand="assemble";; 
        2)          subcommand="submit";; 
        3)          subcommand="check";;    
        4)          subcommand="cluster";;
        5)          subcommand="all";;
        i)          export reads="$OPTARG";;
        o)          export baseDir="$OPTARG";;
        m)          export mem="$OPTARG";;
        n)          export cpu="$OPTARG";;
        p)          export maxProc="$OPTARG";;
        k)          export kmer="$OPTARG";;
        c)          export clust="$OPTARG";;
        h)          displayHelp;;
        \?)         printf ""${BLUE}"Invalid option: -"$OPTARG"\n\n"${NC}"" >&2
                    # echo "Invalid option: -"$OPTARG"" >&2
                    displayHelp && exit 1;;
        :)          printf ""${BLUE}"Option -"$OPTARG" requires an argument.\n\n"${NC}"" >&2
                    # echo "Option -"$OPTARG" requires an argument." >&2
                    displayHelp && exit 1;;
    esac
done

shift $(($OPTIND - 1))


function check_options()
{
    #check subcommand
    if [[ "$subcommand" != "all" ]] && [[ "$subcommand" != "assemble" ]] && [[ "$subcommand" != "submit" ]] \
            && [[ "$subcommand" != "check" ]] && [[ "$subcommand" != "cluster" ]]; then
        printf ""${BLUE}"Please use a valid subcommand"${NC}"\n" >&2
        displayHelp && exit 1
    fi

    # Check input folder
    if [ -z "$reads" ]; then
        printf ""${BLUE}"An input folder is required"${NC}"\n" >&2
        displayHelp && exit 1
    fi

    # check if fastq.gz files are presnet in input folder
    input_files_count=$(find "$reads" -name "*.fastq.gz" | wc -l)
    if [[ "$input_files_count" -lt 2 ]]; then
        printf ""${BLUE}"No Illumina paired-end sequencing files found in "${reads}""${NC}"\n" >&2
        displayHelp && exit 1
    fi

    # check if paired-end
    if [ $(("$input_files_count"%2)) -ne 0 ]; then
        printf ""${BLUE}"Sequencing files do not seems to be paired-end\n"  >&2
        printf "Make sure you have 2 files per sample, i.e. a \"R1\" and a \"R2\""${NC}"\n" >&2
        displayHelp && exit 1
    fi

    # Check output folder
    if [ -z "$baseDir" ]; then
        printf ""${BLUE}"An output folder is required"${NC}"\n" >&2
        displayHelp && exit 1
    fi

    # Check if number of core is exceding number of cores available
    if [ "$cpu" -gt "$(nproc)" ]; then
        printf ""${BLUE}"Number of cores entered ("$cpu") excedes total number of cores available ("$(nproc)")"${NC}"\n" >&2
        printf ""${BLUE}"Number of cores has been set to \""$(nproc)"\""${NC}"\n" >&2
        export cpu="$(nproc)"
    fi

    # Check maximum number of parallel processes
    if [ "$maxProc" -gt "$cpu" ]; then 
        printf ""${BLUE}"Number of parallel processes entered ("$maxProc") excedes total number of cores available ("$(nproc)")"${NC}"\n" >&2
        printf ""${BLUE}"Number of parallel processes has been set to \""$(nproc)"\""${NC}"\n" >&2
        export maxProc="$cpu"
    fi

    # Check if memory is more that 85% of memory available
    if [ "$mem" -gt "$maxMem" ]; then
        printf ""${BLUE}"Memory requested ("$mem") excedes 85% of available memory ("$maxMem")"${NC}"\n" >&2
        printf ""${BLUE}"Number of cores has been set to \""$maxMem"\""${NC}"\n" >&2
        export mem="$maxMem"
    fi
    export memJava="-Xmx"$mem"g"
    memCdHit=$((mem*1000))

    # Check if kmer values are good
    # Comma-separated list of k-mer sizes to be used
    # All values must be odd, less than 128 and listed in ascending order
    IFS=',' read -r -a kmer_array <<< "$kmer"
    j=0
    for i in "${kmer_array[@]}"; do
        if [ $((i%2)) -eq 0 ]; then  # even
            printf ""${BLUE}"kmer value(s) must be odd"${NC}"\n" >&2
            displayHelp && exit 1
        elif [ "$i" -gt 127 ]; then
            printf ""${BLUE}"kmer value(s) must not exceed 127"${NC}"\n" >&2
            displayHelp && exit 1
        elif [ "${#kmer_array[@]}" -gt 1 ] && [ "$i" -lt "$j" ]; then
            printf ""${BLUE}"kmer values must be listed in ascending order"${NC}"\n" >&2
            displayHelp && exit 1
        fi
        j="$i"
    done

    # Clustering settings
    IFS=',' read -ra clust_array <<< "$clust"
    seq_id="${clust_array[0]}"  # -c
    len_cutoff="${clust_array[1]}"  # -s
    word_size="${clust_array[2]}"  # -n

    # http://weizhongli-lab.org/lab-wiki/doku.php?id=cd-hit-user-guide
    # -c  sequence identity threshold, default 0.9
    #     this is the default cd-hit's "global sequence identity" calculated as:
    #     number of identical amino acids in alignment
    #     divided by the full length of the shorter sequence
    # -s  length difference cutoff, default 0.0
    #     if set to 0.9, the shorter sequences need to be
    #     at least 90% length of the representative of the cluster

    # For DNAs:

        # Word size 10-11 is for thresholds 0.95 ~ 1.0
        # Word size 8,9 is for thresholds 0.90 ~ 0.95
        # Word size 7 is for thresholds 0.88 ~ 0.9
        # Word size 6 is for thresholds 0.85 ~ 0.88
        # Word size 5 is for thresholds 0.80 ~ 0.85
        # Word size 4 is for thresholds 0.75 ~ 0.8

    # Folder where the PST was cloned from GitHub
    export pst_path="$( cd "$(dirname "$0")" ; pwd -P )"  # script path
}

function create_folder_structure()
{
    fastq=""${baseDir}"/fastq"
    export logs=""${baseDir}"/logs"
    export trimmed=""${baseDir}"/trimmed"
    export merged=""${baseDir}"/merged"
    export assembly=""${baseDir}"/assembly"
    export phaster=""${baseDir}"/phaster"
    qiime=""${baseDir}"/qiime"

    # Create folders if do not exist
    # "||" if test is false
    # "&&" if test is true
    [ -d "$baseDir" ] || mkdir -p "$baseDir"
    [ -d "$fastq" ] || mkdir -p "$fastq"
    [ -d "$logs" ] || mkdir -p "$logs"
    [ -d "$trimmed" ] || mkdir -p "$trimmed"
    [ -d "$merged" ] || mkdir -p "$merged"
    [ -d "$assembly" ] || mkdir -p "$assembly"
    [ -d "$phaster" ] || mkdir -p "$phaster"
    [ -d "$qiime" ] || mkdir -p "$qiime"
}

function initiate_log()
{
    echo -e "$(date)\n" | tee "${logs}"/log.txt
    echo -e "User: $(whoami)" | tee -a "${logs}"/log.txt
    echo -e "Processors: "$cpu"" | tee -a "${logs}"/log.txt
    echo -e "Memory: "$mem"GB" | tee -a "${logs}"/log.txt

    # Pipeline version
    echo -e "\n"$0" v"${version}"\n" | tee -a "${logs}"/log.txt
}

function check_dependencies()
{
    # Java
    if hash java 2>/dev/null; then 
        java -version 2>&1 1>/dev/null | grep "java version" | tr -d '"' | tee -a "${logs}"/log.txt
    else
        echo >&2 "java was not found. Aborting." | tee -a "${logs}"/log.txt
        exit 1
    fi

    # GNU parallel
    if hash parallel 2>/dev/null; then 
        parallel --version | head -n 1 | tee -a "${logs}"/log.txt
    else
        echo >&2 "GNU parallel was not found. Aborting." | tee -a "${logs}"/log.txt
        exit 1
    fi

    # Python3
    if hash python3 &>/dev/null; then
        python3 --version | tee -a "${logs}"/log.txt
    else
        echo >&2 "python3 was not found. Aborting." | tee -a "${logs}"/log.txt
        exit 1
    fi

    # Pigz
    if hash pigz 2>/dev/null; then 
        pigz --version | tee -a "${logs}"/log.txt
    else
        echo >&2 "pigz was not found. Aborting." | tee -a "${logs}"/log.txt
        exit 1
    fi

    # BBDuk
    if hash bbduk.sh 2>/dev/null; then 
        bbduk.sh -v 2>&1 1>/dev/null | grep "version" | tee -a "${logs}"/log.txt
    else
        echo >&2 "bbduk.sh was not found. Aborting." | tee -a "${logs}"/log.txt
        exit 1
    fi

    # BBmap
    if hash bbmerge.sh 2>/dev/null; then 
        bbmerge.sh -v 2>&1 1>/dev/null | grep "version" | tee -a "${logs}"/log.txt
    else
        echo >&2 "bbmerge.sh was not found. Aborting." | tee -a "${logs}"/log.txt
        exit 1
    fi

    # SPAdes
    if hash spades.py 2>/dev/null; then
        spades.py -v 2>&1 1>/dev/null | tee -a "${logs}"/log.txt
    else
        echo >&2 "spades.py was not found. Aborting." | tee -a "${logs}"/log.txt
        exit 1
    fi

    # CD-HIT-EST
    if hash cd-hit-est 2>/dev/null; then
        cd-hit-est -h | head -n 1 | tr -d "=" | sed 's/^[ \t]*//;s/[ \t]*$//' | tee -a "${logs}"/log.txt
    else
        echo >&2 "cd-hit-est was not found. Aborting." | tee -a "${logs}"/log.txt
        exit 1
    fi

    # Biom
    if hash biom 2>/dev/null; then
        biom --version | tee -a "${logs}"/log.txt
    else
        echo >&2 "biom was not found. Aborting." | tee -a "${logs}"/log.txt
        exit 1
    fi

    # QIMME
    source activate qiime1
    if [ $? -eq 0 ]; then  # environment was found and activated without error
        v=$(print_qiime_config.py -t | grep "QIIME script version" | cut -d $'\t' -f 2)
        echo "QIIME v"${v}""
        source deactivate
    else
        echo >&2 "QIIME was not found. Aborting." | tee -a "${logs}"/log.txt
        exit 1
    fi
}

function trim()
{
    r1="$1"
    r2=$(echo "$r1" | sed 's/_R1/_R2/')   # 2014-SEQ-0729_S5_L001_R1_001.fastq.gz

    #find bbduk.sh installation path
    bb_path=$(dirname $(which bbduk.sh))

    bbduk.sh "$memJava" \
        threads=$((cpu/maxProc)) \
        in1="$r1" \
        in2="$r2" \
        ref="${bb_path}"/resources/adapters.fa \
        ktrim=r k=23 mink=11 hdist=1 tbo tpe \
        qtrim=lr trimq=10 \
        minlen=64 \
        out1="${trimmed}"/"${sample}"_Trimmed_R1.fastq.gz \
        out2="${trimmed}"/"${sample}"_Trimmed_R2.fastq.gz \
        pigz=t \
        unpigz=t \
        2> >(tee "${logs}"/trimming/"${sample}".txt)
}

function merge()
{
    t1="${trimmed}"/"${sample}"_Trimmed_R1.fastq.gz
    t2=$(echo "$t1" | sed 's/_R1/_R2/')

    bbmerge.sh "$memJava" \
        threads=$((cpu/maxProc)) \
        in1="$t1" \
        in2="$t2" \
        out="${merged}"/"${sample}"_merged.fastq.gz \
        outu1="${merged}"/"${sample}"_unmerged_R1.fastq.gz \
        outu2="${merged}"/"${sample}"_unmerged_R2.fastq.gz \
        pigz=t \
        unpigz=t \
        2> >(tee -a "${logs}"/merging/"${sample}".txt)
}

function assemble()
{
    m="${merged}"/"${sample}"_merged.fastq.gz
    u1="${merged}"/"${sample}"_unmerged_R1.fastq.gz
    u2="${merged}"/"${sample}"_unmerged_R2.fastq.gz

    # Create a separate output directory for each sample
    spadesOut=""${assembly}"/"${sample}""
    [ -d "$spadesOut" ] || mkdir -p "$spadesOut"

    spades.py \
        -t $((cpu/maxProc)) \
        -m $((mem/maxProc)) \
        -k "$kmer" \
        --careful \
        --s1 "$m" \
        --pe1-1 "$u1" \
        --pe1-2 "$u2" \
        -o "$spadesOut" &> /dev/null

    # Remane scaffolds.fasta
    mv "${spadesOut}"/scaffolds.fasta "${spadesOut}"/"${sample}"_assembly.fasta

    # Do some cleanup after assembly (remove temporary files)
    find "$spadesOut" | awk 'NR > 1' | grep -v "${spadesOut}"/"${sample}"_assembly.fasta | xargs rm -rf
}

function TrimMergeAssemble()
{
    sample=$(basename "$1" | cut -d '_' -f 1)

    trim "$1"
    merge
    # Free up some disk space; remove trimmed reads
    rm "${trimmed}"/"${sample}"*
    assemble
    # Free up some disk space; remove merged reads
    rm  "${merged}"/"${sample}"*
}

# Make functions available to parallel
# -f is to export functions
export -f trim
export -f merge
export -f assemble
export -f TrimMergeAssemble

function assemble_from_fastq()
{
    # Populate fastq folder with symbolic links pointing to paired-end fastq files
    find -L "$reads" -type f -name "*.fastq.gz" \
        | parallel "[ -e "${fastq}"/{/} ] || ln -s {} "${fastq}"/{/}"
        # {/} means basename using GNU parallel synthax 
        # https://www.gnu.org/software/parallel/man.html

    #Create report output directory
    [ -d "${logs}"/trimming ] || mkdir -p "${logs}"/trimming
    [ -d "${logs}"/merging ] || mkdir -p "${logs}"/merging

    # Run trimming on multiple samples in parallel
    find -L "$fastq" -type f -name "*.fastq.gz" -name "*_R1*" \
        | parallel  --env TrimMergeAssemble \
                    --env cpu \
                    --env maxProc \
                    --env memJava \
                    --env trimmed \
                    --env merged \
                    --env assembly \
                    --env mem \
                    --env kmer \
                    --env logs \
                    --jobs "$maxProc" \
                    "TrimMergeAssemble {}"
}

function assemblyTrimm()
    {
        sample=$(cut -d '_' -f 1 <<< $(basename "$1"))

        # http://phaster.ca/instructions
        if [ $(cat "$1" | grep -Ec "^>") -gt 1 ]; then  # If more than one contig
            # Remove contigs smaller than 2000 bp from assembly
            perl "${pst_path}"/removesmallscontigs.pl \
                2000 \
                "$1" \
                > "${phaster}"/assemblies/"${sample}"_trimmed2000.fasta
        elif [ $(cat "$1" | grep -Ec "^>") -eq 1 ]; then  # If only one contig
            # Check if contig is at least 1500 bp
            seqlen=$(cat "$1" | awk '!/^>/ {l+=length($0)} END {print l}')
            if [ "$seqlen" -lt 1500 ]; then
                echo "Assembly is one contig, but smaller than 1500bp! Skipping."
            else
                ln -s "$1" "${phaster}"/assemblies/"${sample}".fasta
            fi
        else
            echo "No assembly for "$sample""  # Should not get here!
            exit 1
        fi
    }

# Make function available to parallel
export -f assemblyTrimm  # -f is to export functions

# It seems that submitting the assemblies with wget is faster than using python
function phasterSubmit ()
{
    sample=$(basename "$1" | cut -d '_' -f 1)

    # {"job_id":"ZZ_7aed0446a6","status":"You're next!..."}
    wget --post-file="$i" \
        http://phaster.ca/phaster_api?contigs=1 \
        -O "${phaster}"/"${sample}".json \
        -o "${phaster}"/"${sample}"_wget.log
}

function submit_to_PHASTER_server()
{
    [ -d "${phaster}"/assemblies ] || mkdir -p "${phaster}"/assemblies

    # Run trimming on multiple assemblies in parallel
    find "$assembly" -type f -name "*_assembly.fasta" \
        | parallel  --bar \
                    --env assemblyTrimm \
                    --env phaster \
                    --env pst_path \
                    'assemblyTrimm {}'

    # Submit to phaster sequencially
    counter=0
    total=$(find "${phaster}"/assemblies -type f -name "*.fasta" | wc -l)
    for i in $(find "${phaster}"/assemblies -type f -name "*.fasta"); do
        let counter+=1
        sample=$(cut -d '_' -f 1 <<< $(basename "$1"))
        echo -ne "Submitting "${sample}" ("${counter}"/"${total}")"\\r
        phasterSubmit "$i"
    done
}

function check_PHASTER_for_results()
{
   # Check if results are ready. If so, fetch them, else wait and try again
    python3 "${pst_path}"/checkPhasterServer.py \
        --check \
        -i "${phaster}"/assemblies \
        -o "$phaster" 
}


function extract_fasta()
{
    sample=$(cut -d '_' -f 1 <<< $(basename "$1"))

    # Only get the fasta file out of the zip
    unzip -p \
        -j "${phaster}"/"${sample}"_phaster.zip \
        "phage_regions.fna" \
        > "${phaster}"/"${sample}"_phages.fasta

    # Add sample name and entry number to fasta header
    sed -i "s/^>/>"${sample}"_/" "${phaster}"/"${sample}"_phages.fasta
}

export -f extract_fasta

function check_n_download_PHASTER_results()
{
    check_PHASTER_for_results

    find "$phaster" -type f -name "*_phaster.zip" |
    parallel    --bar \
                --env extract_fasta \
                --env phaster \
                'extract_fasta {}'
}

function cluster_sequences()
{
    # Concatenate all phage sequences in a single file
    for i in $(find "$phaster" -type f -name "*_phages.fasta"); do
        cat "$i" >> "${phaster}"/phages_all.fasta
    done

    # cluster sequences
    cd-hit-est \
        -i "${phaster}"/phages_all.fasta \
        -o "${phaster}"/phages_clustered_c"${seq_id}"_s"${len_cutoff}".fasta \
        -c "$seq_id" \
        -s "$len_cutoff" \
        -n "$word_size" \
        -T "$cpu" \
        -M "$memCdHit" \
        -d 0
}

function create_biom_table()
{
    # Create a list of samples from the phages_all.fasta file
    # Better if tha just listing the *_phages.fasta files, in case multiple 
    # phages_all.fasta are merged
    sampleList=($(cat "${phaster}"/phages_all.fasta | \
                grep -e "^>" | cut -d "_" -f 1 | \
                tr -d ">" | sort | uniq))
    echo "${sampleList[@]}" | tr " " "\n" > "${phaster}"/sampleList.txt

    # Convert CD-HIT-EST ".clstr" output file to OTU table
    # Usage: perl cdHitClstr2table.pl -s sampleList -c cd-hit.clstr -o outputTable.tsv
    perl "${pst_path}"/cdHitClstr2table.pl \
        -s "${phaster}"/sampleList.txt \
        -i "${phaster}"/phages_clustered_c"${seq_id}"_s"${len_cutoff}".fasta.clstr \
        -o "${phaster}"/phages_clustered_c"${seq_id}"_s"${len_cutoff}".tsv
}

function update_phage_database()
{
    # Check for latest prophage database from PHASTER website
    if [ -e "${pst_path}"/prophage_virus.db ]; then
        # cd "$pst_path"
        wget \
            --timestamp \
            http://phaster.ca/downloads/prophage_virus.db \
            -q \
            --show-progress \
            -P "$pst_path"
    else
        wget \
            http://phaster.ca/downloads/prophage_virus.db \
            --show-progress \
            -P "$pst_path"
    fi

    makeblastdb \
        -in "${pst_path}"/prophage_virus.db \
        -dbtype prot
}

function blast_clusted_sequences()
{
    # Identify the phage of the repesentative sequences from CD-HIT-EST
    # ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.fna.tar.gz
    blastn \
        -query "${phaster}"/phages_clustered_c"${seq_id}"_s"${len_cutoff}".fasta \
        -db "${pst_path}"/prophage_virus.db \
        -out "${phaster}"/clusterID_c"${seq_id}"_s"${len_cutoff}".blastn.tsv \
        -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
        -num_threads "$cpu" \
        -evalue "1e-10" \
        -culling_limit 5

    # Best hit only
    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" \
        >  "${phaster}"/clusterID_c"${seq_id}"_s"${len_cutoff}".blastn.tsv.tmp

    cat  "${phaster}"/clusterID_c"${seq_id}"_s"${len_cutoff}".blastn.tsv \
        | sort -t $'\t' -k1,1 -k3,3r \
        | sort -t $'\t' -uk1,1 \
        >>  "${phaster}"/clusterID_c"${seq_id}"_s"${len_cutoff}".blastn.tsv.tmp

    mv  "${phaster}"/clusterID_c"${seq_id}"_s"${len_cutoff}".blastn.tsv.tmp \
         "${phaster}"/clusterID_c"${seq_id}"_s"${len_cutoff}".bestHit.blastn.tsv

    # Add header to blast output
    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" \
        > "${phaster}"/tmp.txt
    cat "${phaster}"/clusterID_c"${seq_id}"_s"${len_cutoff}".blastn.tsv \
        >> "${phaster}"/tmp.txt
    mv "${phaster}"/tmp.txt "${phaster}"/clusterID_c"${seq_id}"_s"${len_cutoff}".blastn.tsv
}

function identify_phages()
{
    update_phage_database
    blast_clusted_sequences
}

function run_qiime()
{
    # Create biom v1 (json) file
    biom convert \
        -i "${phaster}"/phages_clustered_c"${seq_id}"_s"${len_cutoff}".tsv \
        -o "${qiime}"/phages_clustered_c"${seq_id}"_s"${len_cutoff}".biom \
        --table-type="OTU table" \
        --to-json

    # Activate python virtual environment for QIIME
    source activate qiime1

    # Beta Diversity (non-phylogenetic)
    beta_diversity.py \
        -i "${qiime}"/phages_clustered_c"${seq_id}"_s"${len_cutoff}".biom \
        -m euclidean \
        -o "${qiime}"/beta_div_non-phylo/

    # N-J tree from distance matrix from beta-diversity
    neighbor_joining.py \
        -i "${qiime}"/beta_div_non-phylo/euclidean_phages_clustered_c"${seq_id}"_s"${len_cutoff}".txt \
        -o "${qiime}"/nj_euclidean_phages_clustered_c"${seq_id}"_s"${len_cutoff}".tree

    #Deactivate the python virtual environment
    source deactivate
}

function get_sample_clustering()
{
    cluster_sequences
    create_biom_table
    # identify_phages
    run_qiime
}


############
#          #
#   MAIN   #
#          #
############


check_options

case "$subcommand" in
    assemble)
        create_folder_structure
        initiate_log
        check_dependencies
        assemble_from_fastq
        ;;
    submit)
        create_folder_structure
        submit_to_PHASTER_server
        ;;
    check)
        create_folder_structure
        check_n_download_PHASTER_results
        ;;
    cluster)
        create_folder_structure
        get_sample_clustering
        ;;
    all)
        create_folder_structure
        initiate_log
        check_dependencies
        assemble_from_fastq
        submit_to_PHASTER_server
        check_n_download_PHASTER_results
        get_sample_clustering
        ;;
    \?)
        printf ""${BLUE}"Invalid subcommand: "$subcommand""${NC}"" >&2
        displayHelp
        exit 1
        ;;
esac
