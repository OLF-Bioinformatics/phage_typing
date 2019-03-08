#!/bin/bash

###################################
#                                 #
#    Prophage Sequence Typing     #
#                                 #
###################################


author='duceppemo'
version='0.1.1'


##############
#            #
#    Note    #
#            #
##############


# The pipleline will take as sample name everything that is before the first "_"
# in the file name. Thus, it requires that the "_" is present in your file name,
# which is usually the case with Illumina reads.


######################
#                    #
#    User Defined    #
#                    #
######################


# Folder where the PST was cloned from GitHub
export pst_path=""${HOME}"/prog/phage_typing"

# Analysis folder
baseDir=""${HOME}"/analyses/my_analysis"

# Folder where all raw Illumina paired-end reads are located
# The pipeline will look for ".fastq.gz" files
reads="/media/30tb_raid10/data/my_data"

# Phage blast database
phage_db="/media/30tb_raid10/db/blast/prophages/viruses_Dec_2018.fasta"

# Maximum number of cores used per sample for parallel processing
#A highier value reduces the memory footprint.
export maxProc=8

#k-mer size for SPAdes assembler (must be odd number(s))
#Should be smaller that minimum trimmed read length
export kmer="21,33,55,77,99,127"

# Clustering settings
seq_id=0.99  # -c
len_cutoff=0.99  # -s
word_size=11  # -n

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


######################
#                    #
#     Resources      #
#                    #
######################


# Computer performance
export cpu=$(nproc) #total number of cores
export mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
export memJava="-Xmx"$mem"g"
memCdHit=$((mem*1000))


#######################
#                     #
#   Data Stucture     #
#                     #
#######################


# Folder structure
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


#######################
#                     #
#   Initiating log    #
#                     #
#######################


#Date
echo -e "$(date)\n" | tee "${logs}"/log.txt
echo -e "User: $(whoami)" | tee -a "${logs}"/log.txt
echo -e "Processors: "$cpu"" | tee -a "${logs}"/log.txt
echo -e "Memory: "$mem"G" | tee -a "${logs}"/log.txt

#pipeline version
echo -e "\nphage_diversity_HT.sh version "${version}"\n" | tee -a "${logs}"/log.txt  # $0

#check if depenencies are installed
#if so, log version

#java
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


########################
#                      #
#     Fastq files      #
#                      #
########################


# Populate fastq folder with symbolic links pointing to paired-end fastq files
find -L "$reads" -type f -name "*.fastq.gz" \
    | parallel --bar "ln -s {} "${fastq}"/{/}"
    # {/} means basename using GNU parallel synthax 
    # https://www.gnu.org/software/parallel/man.html


########################################
#                                      #
#     Trimming/Merging/Assembling      #
#                                      #
########################################


function TrimMergeAssemble()
{
    # Trim
    r1="$1"
    r2=$(echo "$r1" | sed 's/_R1/_R2/')   # 2014-SEQ-0729_S5_L001_R1_001.fastq.gz
    sample=$(basename "$r1" | cut -d '_' -f 1) 

    bbduk.sh "$memJava" \
        threads=$((cpu/maxProc)) \
        in1="$r1" \
        in2="$r2" \
        ref=nextera.fa.gz \
        ktrim=r k=23 mink=11 hdist=1 tbo tpe \
        qtrim=lr trimq=10 \
        minlen=64 \
        out1="${trimmed}"/"${sample}"_Trimmed_1P.fastq.gz \
        out2="${trimmed}"/"${sample}"_Trimmed_2P.fastq.gz \
        pigz=t \
        unpigz=t \
        2> >(tee "${logs}"/trimming/"${sample}".txt)


    # Merge
    t1="${trimmed}"/"${sample}"_Trimmed_1P.fastq.gz
    t2=$(echo "$t1" | sed 's/_1P/_2P/')

    bbmerge.sh "$memJava" \
        threads=$((cpu/maxProc)) \
        in1="$t1" \
        in2="$t2" \
        out="${merged}"/"${sample}"_merged.fastq.gz \
        outu1="${merged}"/"${sample}"_unmerged_1P.fastq.gz \
        outu2="${merged}"/"${sample}"_unmerged_2P.fastq.gz \
        pigz=t \
        unpigz=t \
        2> >(tee -a "${logs}"/merging/"${sample}".txt)

    # Free up some disk space; remove trimmed reads
    rm "${trimmed}"/"${sample}"*


    # Assemble
    # Create a separate output directory for each sample
    spadesOut=""${assembly}"/"${sample}""
    [ -d "$spadesOut" ] || mkdir -p "$spadesOut"

    spades.py \
        -t $((cpu/maxProc)) \
        -m $((mem/maxProc)) \
        -k "$kmer" \
        --careful \
        --s1 "${merged}"/"${sample}"_merged.fastq.gz \
        --pe1-1 "${merged}"/"${sample}"_unmerged_1P.fastq.gz \
        --pe1-2 "${merged}"/"${sample}"_unmerged_2P.fastq.gz \
        -o "$spadesOut" &> /dev/null

    # Remane scaffolds.fasta
    mv "${spadesOut}"/scaffolds.fasta "${spadesOut}"/"${sample}"_assembly.fasta

    # Do some cleanup after assembly (remove temporary files)
    find "$spadesOut" | awk 'NR > 1' | grep -v "${spadesOut}"/"${sample}"_assembly.fasta | xargs rm -rf

    # Free up some disk space; remove merged reads
    rm  "${merged}"/"${sample}"*
}

# Make function available to parallel
export -f TrimMergeAssemble  # -f is to export functions

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


#########################
#                       #
#   Assembly trimming   #
#                       #
#########################


[ -d "${phaster}"/assemblies ] || mkdir -p "${phaster}"/assemblies

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

# Run trimming on multiple assemblies in parallel
find "$assembly" -type f -name "*_assembly.fasta" \
    | parallel  --bar \
                --env assemblyTrimm \
                --env phaster \
                --env pst_path \
                'assemblyTrimm {}'



###############
#             #
#   PHASTER   #
#             #
###############


# It seems that submitting the assemlies with wget is faster than using python
function phasterSubmit ()
{
    sample=$(basename "$1" | cut -d '_' -f 1)

    # {"job_id":"ZZ_7aed0446a6","status":"You're next!..."}
    wget --post-file="$i" \
        http://phaster.ca/phaster_api?contigs=1 \
        -O "${phaster}"/"${sample}".json \
        -o "${phaster}"/"${sample}"_wget.log
}

# Submit to phaster sequencially
counter=0
total=$(find "${phaster}"/assemblies -type f -name "*.fasta" | wc -l)
for i in $(find "${phaster}"/assemblies -type f -name "*.fasta"); do
    let counter+=1
    sample=$(cut -d '_' -f 1 <<< $(basename "$1"))
    echo -ne "Submitting "${sample}" ("${counter}"/"${total}")"\\r
    phasterSubmit "$i"
done


# Submit assemblies to PHASTER server and fetch results when ready
python3 "${pst_path}"/checkPhasterServer.py \
    --check \
    -i "${phaster}"/assemblies \
    -o "$phaster"

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

find "$phaster" -type f -name "*_phaster.zip" |
parallel    --bar \
            --env extract_fasta \
            --env phaster \
            'extract_fasta {}'


#######################
#                     #
#   Phage diversity   #
#                     #
#######################


for i in $(find "$phaster" -type f -name "*_phages.fasta"); do
    cat "$i" >> "${phaster}"/phages_all.fasta
done

# Cluster similar phages together with CD-HIT-EST
cd-hit-est \
    -i "${phaster}"/phages_all.fasta \
    -o "${phaster}"/phages_clustered_c"${seq_id}"_s"${len_cutoff}".fasta \
    -c "$seq_id" \
    -s "$len_cutoff" \
    -n "$word_size" \
    -T "$cpu" \
    -M "$memCdHit" \
    -d 0


#######################
#                     #
#   Phage diversity   #
#                     #
#######################


# Create a list of samples
sampleList=($(ls "${phaster}"/*_phages.fasta | sed -e 's/_phages.fasta//' -e "s%$phaster/%%g"))
echo "${sampleList[@]}" | tr " " "\n" > "${phaster}"/sampleList.txt

# Convert CD-HIT-EST ".clstr" output file to OTU table
# Usage: perl cdHitClstr2table.pl -s sampleList -c cd-hit.clstr -o outputTable.tsv
perl "${pst_path}"/cdHitClstr2table.pl \
    -s "${phaster}"/sampleList.txt \
    -i "${phaster}"/phages_clustered_c"${seq_id}"_s"${len_cutoff}".fasta.clstr \
    -o "${phaster}"/phages_clustered_c"${seq_id}"_s"${len_cutoff}".tsv


#############
#           #
#   BLAST   #
#           #
#############


# Identify the phage of the repesentative sequences from CD-HIT-EST
# ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.fna.tar.gz
blastn \
    -query "${phaster}"/phages_clustered_c99_s90.fasta \
    -db "$phage_db" \
    -out "${phaster}"/clusterID_c99_s90.blastn.tsv \
    -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
    -num_threads "$cpu" \
    -evalue "1e-10" \
    -culling_limit 5

# Best hit only
blastn \
    -query "${phaster}"/phages_clustered_c"${seq_id}"_s"${len_cutoff}".fasta \
    -db "$phage_db" \
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
mv "${phaster}"/tmp.txt "${phaster}"/clusterID_c9"${seq_id}"_s"${len_cutoff}".blastn.tsv


#############
#           #
#   QIIME   #
#           #
#############


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

# UPGM tree from distance matrix from beta-diversity
upgma_cluster.py \
    -i "${qiime}"/beta_div_non-phylo/euclidean_phages_clustered_c"${seq_id}"_s"${len_cutoff}".txt \
    -o "${qiime}"/euclidean_phages_clustered_c"${seq_id}"_s"${len_cutoff}".tree

# N-J tree from distance matrix from beta-diversity
neighbor_joining.py \
    -i "${qiime}"/beta_div_non-phylo/euclidean_phages_clustered_c"${seq_id}"_s"${len_cutoff}".txt \
    -o "${qiime}"/nj_euclidean_phages_clustered_c"${seq_id}"_s"${len_cutoff}".tree

#Deactivate the python virtual environment
source deactivate
