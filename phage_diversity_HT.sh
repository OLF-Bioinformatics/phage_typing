#!/bin/bash


######################
#                    #
#    User Defined    #
#                    #
######################


#Analysis folder
baseDir=""${HOME}"/analyses/salmonella_walid_extra3"

#reads
reads="/media/6tb_raid10/data/salmonella_walid_extra3"

#program location
export prog=""${HOME}"/prog"


#Maximum number of cores used per sample for parallel processing
#A highier value reduces the memory footprint.
export maxProc=3

#k-mer size for SPAdes assembler (must be odd number(s))
#Should be smaller that minimum trimmed read length
export kmer="21,33,55,77,99,127"


######################
#                    #
#     Resources      #
#                    #
######################


#computer performance
export cpu=$(nproc) #total number of cores
export mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
export memJava="-Xmx"$mem"g"
memCdHit=$((mem*1000))


#######################
#                     #
#   Data Stucture     #
#                     #
#######################


#Folder structure
fastq=""${baseDir}"/fastq"
export logs=""${baseDir}"/logs"
export trimmed=""${baseDir}"/trimmed"
export merged=""${baseDir}"/merged"
export assembly=""${baseDir}"/assembly"
export phaster=""${baseDir}"/phaster"
qiime=""${baseDir}"/qiime"

#create folders if do not exist
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
echo -e "\nphage_diversity_HT.sh version 0.1.0\n" | tee -a "${logs}"/log.txt  # $0

#check if depenencies are installed
#if so, log version

#java
if hash java 2>/dev/null; then 
    java -version 2>&1 1>/dev/null | grep "java version" | tr -d '"' | tee -a "${logs}"/log.txt
else
    echo >&2 "java was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#BBDuk
if hash bbduk.sh 2>/dev/null; then 
    bbduk.sh -v 2>&1 1>/dev/null | grep "version" | tee -a "${logs}"/log.txt
else
    echo >&2 "bbduk.sh was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#BBmap
if hash bbmerge.sh 2>/dev/null; then 
    bbmerge.sh -v 2>&1 1>/dev/null | grep "version" | tee -a "${logs}"/log.txt
else
    echo >&2 "bbmerge.sh was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#SPAdes
if hash spades.py 2>/dev/null; then
    spades.py -v 2>&1 1>/dev/null | tee -a "${logs}"/log.txt
else
    echo >&2 "spades.py was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#CD-HIT-EST
if hash cd-hit-est 2>/dev/null; then
    cd-hit-est -h | head -n 1 | tr -d "=" | sed 's/^[ \t]*//;s/[ \t]*$//' | tee -a "${logs}"/log.txt
else
    echo >&2 "cd-hit-est was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi


########################
#                      #
#     Fastq files      #
#                      #
########################


#populate read folder with symbolic links pointing to paired-end fastq files
find -L "$reads" -type f -name "*.fastq.gz" \
    | parallel --bar "ln -s {} "${fastq}"/{/}"  # {/} means basename in parallel -> https://www.gnu.org/software/parallel/man.html


########################################
#                                      #
#     Trimming/Merging/Assembling      #
#                                      #
########################################


function TrimMergeAssemble()
{
    #trim
    r1="$1"
    r2=$(echo "$r1" | sed 's/_R1/_R2/')   # 2014-SEQ-0729_S5_L001_R1_001.fastq.gz
    sample=$(basename "$r1" | cut -d '_' -f 1) 

    bbduk.sh "$memJava" \
        threads=$((cpu/maxProc)) \
        in1="$r1" \
        in2="$r2" \
        ref="${prog}"/bbmap/resources/nextera.fa.gz \
        ktrim=r k=23 mink=11 hdist=1 tbo tpe \
        qtrim=lr trimq=10 \
        minlen=64 \
        out1="${trimmed}"/"${sample}"_Trimmed_1P.fastq.gz \
        out2="${trimmed}"/"${sample}"_Trimmed_2P.fastq.gz \
        pigz=t \
        unpigz=t \
        2> >(tee "${logs}"/trimming/"${sample}".txt)


    #merge
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

    #Free up some disk space; remove trimmed reads
    rm  "${trimmed}"/"${sample}"*


    #assemble
    #create a separate output directory for each sample
    spadesOut=""${assembly}"/"${sample}""
    [ -d "$spadesOut" ] || mkdir -p "$spadesOut"

    #Splitting the assembly in two separate processes decreases the total run time of assembly
    #error correction only
    spades.py \
        --only-error-correction \
        -t $((cpu/maxProc)) \
        -m $((mem/maxProc)) \
        -k "$kmer" \
        --careful \
        --s1 "${merged}"/"${sample}"_merged.fastq.gz \
        --pe1-1 "${merged}"/"${sample}"_unmerged_1P.fastq.gz \
        --pe1-2 "${merged}"/"${sample}"_unmerged_2P.fastq.gz \
        -o "${spadesOut}"/ec &> /dev/null

    #corrected reads
    cor_m=$(find "${spadesOut}"/ec/corrected -type f -name "*_merged*cor.fastq.gz")
    cor_1p=$(find "${spadesOut}"/ec/corrected -type f -name "*_unmerged_1P*cor.fastq.gz")
    cor_2p=$(find "${spadesOut}"/ec/corrected -type f -name "*_unmerged_2P*cor.fastq.gz")

    #Assembly only
    spades.py \
        --only-assembler \
        -t $((cpu/maxProc)) \
        -m $((mem/maxProc)) \
        -k "$kmer" \
        --careful \
        --s1 "$cor_m" \
        --pe1-1 "$cor_1p" \
        --pe1-2 "$cor_2p" \
        -o "$spadesOut" &> /dev/null

    #remane scaffolds.fasta
    mv "${spadesOut}"/scaffolds.fasta "${spadesOut}"/"${sample}"_assembly.fasta

    #do some cleanup after assembly (remove temporary files)
    find "$spadesOut" | awk 'NR > 1' | grep -v "${spadesOut}"/"${sample}"_assembly.fasta | xargs rm -rf

    #Free up some disk space; remove merged reads
    rm  "${merged}"/"${sample}"*
}

#make function available to parallel
export -f TrimMergeAssemble  # -f is to export functions

#Create report output directory
[ -d "${logs}"/trimming ] || mkdir -p "${logs}"/trimming
[ -d "${logs}"/merging ] || mkdir -p "${logs}"/merging

#run trimming on multiple samples in parallel
find -L "$fastq" -type f -name "*.fastq.gz" -name "*_R1*" \
    | parallel  --env TrimMergeAssemble \
                --env cpu \
                --env maxProc \
                --env memJava \
                --env prog \
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


function assemblyTrimm()
{
    name=$(basename "$1")
    sample=$(cut -d '_' -f 1 <<< "$name")

    # http://phaster.ca/instructions
    if [ $(cat "$1" | grep -Ec "^>") -gt 1 ]; then  # if more than one contig
        #remove contigs smaller than 2000 bp from assembly
        perl "${prog}"/phage_typing/removesmallscontigs.pl \
            2000 \
            "$1" \
            > "${1%.fasta}"_trimmed2000.fasta
    elif [ $(cat "$1" | grep -Ec "^>") -eq 1 ]; then  # if only one contig
        #remove contigs smaller than 2000 bp from assembly
        perl "${prog}"/phage_typing/removesmallscontigs.pl \
            1500 \
            "$1" \
            > "${1%.fasta}"_trimmed1500.fasta
    else
        echo "No assembly for "$sample""  # Should not get here!
        exit 1
    fi

    #replace old file
    # mv "${1}".tmp "${1}"
}

#make function available to parallel
export -f assemblyTrimm  # -f is to export functions

#run trimming on multiple assemblies in parallel
find "$assembly" -type f -name "*_assembly.fasta" \
    | parallel --env assemblyTrimm --env prog 'assemblyTrimm {}'



###############
#             #
#   PHASTER   #
#             #
###############


function phasterSubmit()
{
    sample=$(basename "$1" | cut -d '_' -f 1)

    wget --post-file="$1" \
        http://phaster.ca/phaster_api?contigs=1 \
        -O "${phaster}"/"${sample}"_query.json \
        -o "${phaster}"/"${sample}"_wget.log
}

# Submit samples to PHASTER
total=$(find "$assembly" -type f -name "*_assembly_trimmed*.fasta" | wc -l)
counter=0
for i in $(find "$assembly" -type f -name "*_assembly_trimmed*.fasta"); do
    let counter=counter+1
    echo -ne "Progress: "${counter}"/"${total}"\r"

    phasterSubmit "$i"
done

function phasterResults()
{
    name=$(basename "$1")
    sample=$(cut -d '_' -f 1 <<< "$name")

    #Retrieve job ID from json file
    # {"job_id":"ZZ_7aed0446a6","status":"You're next!..."}
    jobID=$(cat "${phaster}"/"${sample}"_query.json | cut -d ',' -f 1 | cut -d ":" -f 2 | tr -d '"')
    # echo ""${sample}": "$jobID""  # debug

    #Check if all submission were successful
    while [ -z "$jobID" ]; do  # if no jobID (unsuccessful submission)
        phasterSubmit "$1"  # resubmit the sample
        jobID=$(cat "${phaster}"/"${sample}"_query.json | cut -d ',' -f 1 | cut -d ":" -f 2 | tr -d '"')
    done

    #get status
    wget http://phaster.ca/phaster_api?acc="$jobID" -O "${phaster}"/"${sample}"_status.json
    status=$(cat "${phaster}"/"${sample}"_status.json | cut -d ',' -f 2 | cut -d ":" -f 2 | tr -d '"')

    # echo "PHASTER analysis of "$sample" is "$status""

    #check if PHASTER job is finished running
    while [ "$status" != "Complete" ]; do
        waitTime="30s" # sleep 

        #check status every X time
        echo "Job status of "$sample" is "$status". Checking status back in "$waitTime"." 
        sleep "$waitTime"

        #get status
        wget http://phaster.ca/phaster_api?acc="$jobID" -O "${phaster}"/"${sample}"_status.json

        #check job status
        status=$(cat "${phaster}"/"${sample}"_status.json | cut -d ',' -f 2 | cut -d ":" -f 2 | tr -d '"')
    done

    echo "PHASTER analysis of "$sample" is "$status""

    #get the PHASTER output file
    phasterZip=$(cat "${phaster}"/"${sample}"_status.json | cut -d ',' -f 4 | cut -d ":" -f 2 | tr -d '"')

    # Check if was already downloaded
    if [ ! -s "${phaster}"/"${sample}"_phaster.zip ]; then
        wget "$phasterZip" -O "${phaster}"/"${sample}"_phaster.zip

        #Only get the fasta file out of the zip
        unzip -p \
            -j "${phaster}"/"${sample}"_phaster.zip \
            "phage_regions.fna" \
            > "${phaster}"/"${sample}"_phages.fasta

        #Add sample name and entry number to fasta header
        sed -i "s/^>/>"${sample}"_/" "${phaster}"/"${sample}"_phages.fasta
    fi
}

# #make function available to parallel
# export -f phasterResults  # -f is to export functions

# Parallel downlaoding results in an overload of the CPUs and no more download afer ~123 samples.
# find "$assembly" -type f -name "*_assembly.fasta" \
#     | parallel  --bar \
#                 --delay 0.3 \
#                 --env phasterResults \
#                 --env phaster \
#                 --jobs "$cpu" \
#                 'phasterResults {}'


total=$(find "$assembly" -type f -name "*_assembly.fasta" | wc -l)
counter=0
for i in $(find "$assembly" -type f -name "*_assembly.fasta"); do
    let counter=counter+1
    echo -ne "Progress: "${counter}"/"${total}"\r"
    phasterResults "$i"
done



#######################
#                     #
#   Phage diversity   #
#                     #
#######################


for i in $(find "$phaster" -type f -name "*_phages.fasta"); do
    cat "$i" >> "${phaster}"/phages_all.fasta
done


# http://weizhongli-lab.org/lab-wiki/doku.php?id=cd-hit-user-guide

# For DNAs:

    # Word size 10-11 is for thresholds 0.95 ~ 1.0
    # Word size 8,9 is for thresholds 0.90 ~ 0.95
    # Word size 7 is for thresholds 0.88 ~ 0.9
    # Word size 6 is for thresholds 0.85 ~ 0.88
    # Word size 5 is for thresholds 0.80 ~ 0.85
    # Word size 4 is for thresholds 0.75 ~ 0.8

#cluster similar phages together with CD-HIT-EST
cd-hit-est \
    -i "${phaster}"/phages_all.fasta \
    -o "${phaster}"/phages_clustered.fasta \
    -c 0.9 \
    -s 0.9 \
    -T $(nproc) \
    -M "$memCdHit" \
    -n 8 \
    -d 0


#######################
#                     #
#   Phage diversity   #
#                     #
#######################


# Create a list of samples
sampleList=($(ls "${phaster}"/*_phages.fasta | sed -e 's/_phages.fasta//' -e "s%$phaster/%%g"))
echo "${sampleList[@]}" | tr " " "\n" > "${phaster}"/sampleList.txt

#Convert CD-HIT-EST ".clstr" output file to OTU table
# Usage: perl cdHitClstr2table.pl -s sampleList -c cd-hit.clstr -o outputTable.tsv
perl "${prog}"/phage_typing/cdHitClstr2table.pl \
    -s "${phaster}"/sampleList.txt \
    -i "${phaster}"/phages_clustered.fasta.clstr \
    -o "${phaster}"/phages_clustered.tsv


#############
#           #
#   BLAST   #
#           #
#############


#Identify the phage og the repesentative sequences from CD-HIT-EST
blastx \
    -query "${phaster}"/phages_clustered.fasta \
    -db /media/3tb_hdd/db/blast/prophages/prophage_virus.db \
    -out "${phaster}"/clusterID.blastx \
    -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
    -num_threads $(nproc) \
    -max_target_seqs 1


#############
#           #
#   QIIME   #
#           #
#############


#activate python virtual environment for QIIME
source activate qiime1

#Create biom fole from OTU table
make_otu_table.py \
    -i "${phaster}"/phages_clustered.tsv \
    -o "${qiime}"/phages_clustered.biom


#Deactivate the python virtual environment
source deactivate

