#!/bin/bash


######################
#                    #
#    User Defined    #
#                    #
######################


#Analysis folder
baseDir=""${HOME}"/analyses/assembly_test"

#reads
reads="/media/3tb_hdd/data/salmonella_panel"

#program location
export prog=""${HOME}"/prog"
export scripts=""${HOME}"/scripts"

#Kraken DB
export db="/media/3tb_hdd/db/kraken/refseq_BV_old"

#Maximum number of cores used per sample for parallel processing
#A highier value reduces the memory footprint.
export maxProc="4"

#k-mer size for SPAdes assembler (must be odd number(s))
#Should be smaller that minimum trimmed read length
export kmer="21,33,55,77,99,127"

#Minimun contig length (for PHASTER)
# export minContigLen=2000


#######################
#                     #
#   Data Stucture     #
#                     #
#######################


#Folder structure
fastq=""${baseDir}"/fastq"
export logs=""${baseDir}"/logs"
qc=""${baseDir}"/QC"
export krakenOut=""${qc}"/kraken/raw"
export trimmed=""${baseDir}"/trimmed"
export merged=""${baseDir}"/merged"
export assembly=""${baseDir}"/assembly"
export phaster=""${baseDir}"/phaster"

#create folders if do not exist
# "||" if test is false
# "&&" if test is true
[ -d "$baseDir" ] || mkdir -p "$baseDir"
[ -d "$fastq" ] || mkdir -p "$fastq"
[ -d "$logs" ] || mkdir -p "$logs"
[ -d "$qc" ] || mkdir -p "$qc"
[ -d "$krakenOut" ] || mkdir -p "$krakenOut"
[ -d "$trimmed" ] || mkdir -p "$trimmed"
[ -d "$merged" ] || mkdir -p "$merged"
[ -d "$assembly" ] || mkdir -p "$assembly"
[ -d "$phaster" ] || mkdir -p "$phaster"


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
#   Initiating log    #
#                     #
#######################


#Date
echo -e "$(date)\n" | tee "${logs}"/log.txt
echo -e "User: $(whoami)" | tee -a "${logs}"/log.txt
echo -e "Processors: "$cpu"" | tee -a "${logs}"/log.txt
echo -e "Memory: "$mem"G" | tee -a "${logs}"/log.txt

#software version
echo -e "\nphage_diversity.sh version 0.1.0\n" | tee -a "${logs}"/log.txt  # $0

# kraken -v | grep "Kraken" | tee -a "${logs}"/log.txt
# java -version 2>&1 1>/dev/null | grep "java version" | tee -a "${logs}"/log.txt
# fastqc -v | tee -a "${logs}"/log.txt
# bbduk.sh -v 2>&1 1>/dev/null | grep "version" | tee -a "${logs}"/log.txt
# bbmerge.sh -v 2>&1 1>/dev/null | grep "version" | tee -a "${logs}"/log.txt
# spades.py -v | tee -a "${logs}"/log.txt
# quast.py -v | tee -a "${logs}"/log.txt
# cd-hit-est -h | head -n 1 | tr -d "=" | sed 's/^[ \t]*//;s/[ \t]*$//' | tee -a "${logs}"/log.txt


#check if depenencies are installed
#if so, log version

#Kraken
if hash kraken 2>/dev/null; then  # if installed
    kraken -v | grep "Kraken" | tee -a "${logs}"/log.txt
else
    echo >&2 "kraken was not found. Aborting."
    exit 1
fi

#java
if hash java 2>/dev/null; then 
    java -version 2>&1 1>/dev/null | grep "java version" | tr -d '"' | tee -a "${logs}"/log.txt
else
    echo >&2 "java was not found. Aborting."
    exit 1
fi

#FastQC
if hash fastqc 2>/dev/null; then 
    fastqc -v | tee -a "${logs}"/log.txt
else
    echo >&2 "fastQC was not found. Aborting."
    exit 1
fi

#BBDuk
if hash bbduk.sh 2>/dev/null; then 
    bbduk.sh -v 2>&1 1>/dev/null | grep "version" | tee -a "${logs}"/log.txt
else
    echo >&2 "bbduk.sh was not found. Aborting."
    exit 1
fi

#BBmap
if hash bbmerge.sh 2>/dev/null; then 
    bbmerge.sh -v 2>&1 1>/dev/null | grep "version" | tee -a "${logs}"/log.txt
else
    echo >&2 "bbmerge.sh was not found. Aborting."
    exit 1
fi

#SPAdes
if hash spades.py 2>/dev/null; then
    spades.py -v | tee -a "${logs}"/log.txt
else
    echo >&2 "spades.py was not found. Aborting."
    exit 1
fi

#QUAST
if hash spades.py 2>/dev/null; then
    quast.py -v | tee -a "${logs}"/log.txt
else
    echo >&2 "quast.py was not found. Aborting."
    exit 1
fi

#CD-HIT-EST
if hash cd-hit-est 2>/dev/null; then
    cd-hit-est -h | head -n 1 | tr -d "=" | sed 's/^[ \t]*//;s/[ \t]*$//' | tee -a "${logs}"/log.txt
else
    echo >&2 "cd-hit-est was not found. Aborting."
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

# for b in $(find "$reads" -type f -name "*.fastq.gz") #look in folder recursively.
# do
#   name=$(basename "$b") 
#   ln -s "$b" "${fastq}"/"$name"
# done


####################
#                  #
#   FastQc - Raw   #
#                  #
####################


[ -e  "${qc}"/fastqc/list.txt ] && rm  "${qc}"/fastqc/list.txt
for i in $(find -L "$fastq" -type f -name "*.fastq.gz"); do 
    echo "$i" >> "${qc}"/fastqc/list.txt
done

ARRAY=($( cat  "${qc}"/fastqc/list.txt ))
list=$(echo "${ARRAY[@]}")

[ -d  "${qc}"/fastqc/raw ] || mkdir -p  "${qc}"/fastqc/raw

fastqc \
    --o  "${qc}"/fastqc/raw \
    --noextract \
    --threads "$cpu" \
    $list  # don't put in quotes


####################
#                  #
#   Kraken - Raw   #
#                  #
####################


#put kraken db into memory
cat "${db}"/database.kdb > /dev/null  # About 6min to run


# Can't run Kraken with parallel because it's using the database stored into memory

# function krack()
# {
#     r1="$1"
#     r2=$(echo "$r1" | sed 's/_R1/_R2/')
#     name=$(basename "$r1") #discard the path
#     sample=$(cut <<< "${name}" -d '_' -f 1) #only keep the first part of the file name, before the first "_"

#     #Output folder
#     out=""${krakenOut}"/"${sample}""
#     [ -d "$out" ] || mkdir -p "$out"

#     kraken --db "$db" \
#         --output "${out}"/"${sample}".kraken \
#         --threads $((cpu/maxProc)) \
#         --gzip-compressed \
#         --check-names \
#         --fastq-input \
#         --paired "$r1" "$r2" \
#         &> >(tee "${out}"/"${sample}".kraken.log)

#     #Prepare result for display with Krona
#     cat "${out}"/"${sample}".kraken \
#         | cut -f 2-3 \
#         | ktImportTaxonomy /dev/stdin -o "${out}"/"${sample}".html
# }


# #make function available to parallel
# export -f krack  # -f is to export functions

# #Run Kraken in parallel
# find -L "$fastq" -type f -name "*.fastq.gz" -name "*_R1*" \
#     | parallel --env db --env maxProc --env krack --env krakenOut --jobs "$maxProc" "krack {}"


start=$(date +%s)

#Find fastq reads and get sample
for k in $(find -L "$fastq" -type f -name "*_R1*" -name "*.fastq.gz"); do #find forward (R1) read files recursively
    r1="$k"
    r2=$(echo "$r1" | sed 's/_R1/_R2/')
    name=$(basename "$r1") #discard the path
    sample=$(cut <<< "${name}" -d '_' -f 1) #only keep the first part of the file name, before the first "_"

    #Output folder
    out=""${krakenOut}"/"${sample}""
    [ -d "$out" ] || mkdir -p "$out"

    #run Kraken
    kraken --db "$db" \
        --output "${out}"/"${sample}".kraken \
        --threads "$cpu" \
        --gzip-compressed \
        --check-names \
        --fastq-input \
        --paired "$r1" "$r2" \
        &> >(tee "${out}"/"${sample}".kraken.log)

    #Prepare result for display with Krona
    cat "${out}"/"${sample}".kraken \
        | cut -f 2-3 \
        | ktImportTaxonomy /dev/stdin -o "${out}"/"${sample}".html

    #visualize the resutls in Firefow browser
    # firefox file://"${out}"/"${sample}".html &
done

#Elapsed time
end=$(date +%s)
elapsed=$(($end - $start))
printf "Kraken analysis finished in %dh:%dm:%ds\n" \
    $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60)) | tee -a "${logs}"/log.txt


#####################
#                   #
#     Trimming      #
#                   #
#####################


start=$(date +%s)

function trimm()
{
    #sequence nomenclature:
    # 2014-SEQ-0729_S5_L001_R1_001.fastq.gz
    r1="$1"
    r2=$(echo "$r1" | sed 's/_R1/_R2/')
    name=$(basename "$r1") 
    sample=$(cut -d '_' -f 1 <<< "$name")

    bbduk.sh "$memJava" \
        threads=$((cpu/maxProc)) \
        in1="$r1" \
        in2="$r2" \
        ref="${prog}"/bbmap/resources/nextera.fa.gz \
        ktrim=r k=23 mink=11 hdist=1 tbo tpe \
        qtrim=lr trimq=5 \
        qtrim=w trimq=15 \
        minlen=64 \
        out1="${trimmed}"/"${sample}"_Trimmed_1P.fastq.gz \
        out2="${trimmed}"/"${sample}"_Trimmed_2P.fastq.gz \
        pigz=t \
        unpigz=t \
        2> >(tee "${logs}"/trimming/"${sample}".txt)
}

#make function available to parallel
export -f trimm  # -f is to export functions

#Create report output directory
[ -d "${logs}"/trimming ] || mkdir -p "${logs}"/trimming

#run trimming on multiple samples in parallel
find -L "$fastq" -type f -name "*.fastq.gz" -name "*_R1*" \
    | parallel --env trimm --env cpu --env maxProc --env memJava --env prog --env trimmed --env logs --jobs "$maxProc" "trimm {}"

#Elapsed time
end=$(date +%s)
elapsed=$(($end - $start))
printf "Trimming finished in %dh:%dm:%ds\n" \
    $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60)) | tee -a "${logs}"/log.txt


########################
#                      #
#   FastQc - Trimmed   #
#                      #
########################


[ -e "${qc}"/list.txt ] && rm "${qc}"/fastqc/list.txt
for i in $(find "$trimmed" -type f -name "*.fastq.gz"); do 
    echo "$i" >> "${qc}"/fastqc/list.txt
done

ARRAY=($( cat "${qc}"/fastqc/list.txt ))
list=$(echo "${ARRAY[@]}")

[ -d "${qc}"/fastqc/trimmed ] || mkdir -p "${qc}"/fastqc/trimmed

fastqc \
    --o "${qc}"/fastqc/trimmed \
    --noextract \
    --threads "$cpu" \
    $list  # don't put in quotes


####################
#                  #
#     Merging      #
#                  #
####################


#paired-end read merging
start=$(date +%s)

function merge()
{
    #sequence nomenclature:
    # 2014-SEQ-0729_S5_L001_R1_001.fastq.gz
    r1="$1"
    r2=$(echo "$r1" | sed 's/_1P/_2P/')
    name=$(basename "$r1")
    sample=$(cut -d '_' -f 1 <<< "$name")

    bbmerge.sh "$memJava" \
        threads=$((cpu/maxProc)) \
        in1="$r1" \
        in2="$r2" \
        out="${merged}"/"${sample}"_merged.fastq.gz \
        outu1="${merged}"/"${sample}"_unmerged_1P.fastq.gz \
        outu2="${merged}"/"${sample}"_unmerged_2P.fastq.gz \
        pigz=t \
        unpigz=t \
        2> >(tee -a "${logs}"/merging/"${sample}".txt)
}

#make function available to parallel
export -f merge  # -f is to export functions

#Create report output directory
[ -d "${logs}"/merging ] || mkdir -p "${logs}"/merging

#run paired-end merging on multiple samples in parallel
find -L "$trimmed" -type f -name "*.fastq.gz" -name "*_1P*" \
    | parallel --env merge --env maxProc --env cpu --env memJava --env merged --env logs --jobs "$maxProc" "merge {}"

#Elapsed time
end=$(date +%s)
elapsed=$(($end - $start))
printf "Merging finished in %dh:%dm:%ds\n" \
    $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60)) | tee -a "${logs}"/log.txt

#remove trimmed files
rm "${trimmed}"/*


#######################
#                     #
#   FastQc - Merged   #
#                     #
#######################


# Prepare list of files
[ -e "${qc}"/fastqc/list.txt ] && rm "${qc}"/fastqc/list.txt
for i in $(find "$merged" -type f -name "*.fastq.gz"); do 
    echo "$i" >> "${qc}"/fastqc/list.txt
done
ARRAY=($( cat "${qc}"/fastqc/list.txt ))
list=$(echo "${ARRAY[@]}")

#Create report output directory
[ -d "${qc}"/fastqc/merged ] || mkdir -p "${qc}"/fastqc/merged

#run FastQC
fastqc \
    --o "${qc}"/fastqc/merged \
    --noextract \
    --threads "$cpu" \
    $list  # don't put in quotes

#cleanup
rm "${qc}"/fastqc/list.txt


#####################
#                   #
#     Assembly      #
#                   #
#####################


# "${merged}"/"${sample}"_merged.fastq.gz \
# "${merged}"/"${sample}"_unmerged_1P.fastq.gz \
# "${merged}"/"${sample}"_unmerged_2P.fastq.gz \

#runtime
start=$(date +%s)

function assemble()
{
    path=$(dirname "$1")
    name=$(basename "$1")
    sample=$(cut <<< "$name" -d '_' -f 1)

    r1=""${path}"/"${sample}"_unmerged_1P.fastq.gz"
    r2=$(echo "$r1" | sed 's/_1P/_2P/')

    #create a separate output directory for each sample
    spadesOut=""${assembly}"/"$sample""
    [ -d "$spadesOut" ] || mkdir -p "$spadesOut"

    #Splitting the assembly in two separate processes decreases the total run time of assembly
    #error correction only
    spades.py \
        --only-error-correction \
        -t $((cpu/maxProc)) \
        -m $((mem/maxProc)) \
        -k "$kmer" \
        --careful \
        --s1 "$1" \
        --pe1-1 "$r1" \
        --pe1-2 "$r2" \
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
}


#make function available to parallel
export -f assemble  # -f is to export functions

#run assembly on multiple samples in parallel
find "$merged" -type f -name "*_merged.fastq.gz" \
    | parallel --bar --env maxProc --env assemble --env assembly --env mem --env kmer --env cpu --env scripts --jobs "$maxProc" 'assemble {}'

#Elapsed time
end=$(date +%s)
elapsed=$(($end - $start))
printf "Assembly finished in %dh:%dm:%ds\n" \
    $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60)) | tee -a "${logs}"/log.txt

# Remove merged files
rm "${merged}"/*


#############
#           #
#   Quast   #
#           #
#############


[ -e "${assembly}"/list.txt ] && rm "${assembly}"/list.txt
for i in $(find "$assembly" -type f -name "*_assembly.fasta"); do 
    echo "$i" >> "${assembly}"/list.txt
done

ARRAY=($(cat "${assembly}"/list.txt))
list=$(echo "${ARRAY[@]}")

quast.py \
    -s \
    -o "${assembly}"/quast \
    -t "$cpu" \
    $list  # don't put in quotes


# TODO
# Make quast report on individual assembly


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
        perl "${scripts}"/removesmallscontigs.pl \
            2000 \
            "$1" \
            > "${1}".tmp
    elif [ $(cat "$1" | grep -Ec "^>") -eq 1 ]; then  # if only one contig
        #remove contigs smaller than 2000 bp from assembly
        perl "${scripts}"/removesmallscontigs.pl \
            1500 \
            "$1" \
            > "${1}".tmp
    else
        echo "No assembly for "$sample""  # Should not get here!
        exit 1
    fi

    #replace old file
    mv "${1}".tmp "${1}"
}

#make function available to parallel
export -f assemblyTrimm  # -f is to export functions

#run trimming on multiple assemblies in parallel
find "$assembly" -type f -name "*_assembly.fasta" \
    | parallel --env assemblyTrimm --env scripts 'assemblyTrimm {}'


###############
#             #
#   PHASTER   #
#             #
###############


#Batch submit samples to PHASTER server
#Allow little deley between each sample
function phasterSubmit()
{
    name=$(basename "$1")
    sample=$(cut -d '_' -f 1 <<< "$name")

    # {"job_id":"ZZ_7aed0446a6","status":"You're next!..."}
    wget --post-file="$1" \
        http://phaster.ca/phaster_api?contigs=1 \
        -O "${phaster}"/"${sample}"_query.json \
        -o "${phaster}"/"${sample}"_wget.log
}

#make function available to parallel
export -f phasterSubmit  # -f is to export functions

find "$assembly" -type f -name "*_assembly.fasta" \
    | parallel --delay 30 --env phasterSubmit --env phaster 'phasterSubmit {}'


function phasterResults()
{
    name=$(basename "$1")
    sample=$(cut -d '_' -f 1 <<< "$name")

    #Retrieve job ID from json file
    jobID=$(cat "${phaster}"/"${sample}"_query.json | cut -d ',' -f 1 | cut -d ":" -f 2 | tr -d '"')
    # echo ""${sample}": "$jobID""  # debug

    #Check if all submission were successful
    while [ -z "$jobID" ]; do  # if no jobID (unsuccessful submission)
        phasterSubmit "$1"  # resubmit the sample
        jobID=$(cat "${phaster}"/"${sample}"_query.json | cut -d ',' -f 1 | cut -d ":" -f 2 | tr -d '"')
    done

    #check if PHASTER job is finished running
    status="Submitted"

    # echo "PHASTER analysis of "$sample" is "$status""
    while [ "$status" != "Complete" ]; do
        #check status every 2 minutes
        echo "Job is "$status". Checking status back in 2 minutes." 
        sleep 2m  # sleep 2 minutes

        #get status
        wget http://phaster.ca/phaster_api?acc="$jobID" -O "${phaster}"/"${sample}"_status.json

        #check job status
        status=$(cat "${phaster}"/"${sample}"_status.json | cut -d ',' -f 2 | cut -d ":" -f 2 | tr -d '"')
    done

    echo "PHASTER analysis of "$sample" is "$status""
}
    #get the PHASTER output file
    phasterZip=$(cat "${phaster}"/"${sample}"_status.json | cut -d ',' -f 4 | cut -d ":" -f 2 | tr -d '"')
    wget "$phasterZip" -O "${phaster}"/"${sample}"_phaster.zip

    #Only get the fasta file out of the zip
    unzip -p \
        -j "${phaster}"/"${sample}"_phaster.zip \
        "phage_regions.fna" \
        > "${phaster}"/"${sample}"_phages.fasta

    #Add sample name and entry number to fasta header
    sed -i "s/^>/>"${sample}"_/" "${phaster}"/"${sample}"_phages.fasta
}

#make function available to parallel
export -f phasterResults  # -f is to export functions

find "$assembly" -type f -name "*_assembly.fasta" \
    | parallel --delay 0.3 --env phasterResults --env phaster 'phasterResults {}'









# for i in $(find "$assembly" -type f -name "*_assembly.fasta"); do
#     name=$(basename "$i")
#     sample=$(cut -d '_' -f 1 <<< "$name")

#     #send query
#     wget --post-file="$i" \
#         http://phaster.ca/phaster_api?contigs=1 \
#         -O "${phaster}"/"${sample}"_query.json  # {"job_id":"ZZ_7aed0446a6","status":"You're next!..."}
# done

#     #Retrieve job ID from json file
#     jobID=$(cat "${phaster}"/"${sample}"_query.json | cut -d ',' -f 1 | cut -d ":" -f 2 | tr -d '"')

#     #check if PHASTER job is finished running
#     status="Submitted"
#     while [ "$status" != "Complete" ]; do
#         #check status every 2 minutes
#         echo "Job is "$status". Checking status back in 2 minutes" 
#         sleep 2m  # sleep 2 minutes

#         #get status
#         wget http://phaster.ca/phaster_api?acc="$jobID" -O "${phaster}"/"${sample}"_status.json

#         #check job status
#         # {"job_id":"ZZ_7a1852db39","status":"1072 submissions ahead of yours..."}
#         status=$(cat "${phaster}"/"${sample}"_status.json | cut -d ',' -f 2 | cut -d ":" -f 2 | tr -d '"' | tr -d "}")
#     done

#     echo "PHASTER analysis of "$sample" is "$status""

#     #get the PHASTER output file
#     phasterZip=$(cat "${phaster}"/"${sample}"_status.json | cut -d ',' -f 4 | cut -d ":" -f 2 | tr -d '"')
#     wget "$phasterZip" -O "${phaster}"/"${sample}"_phaster.zip

#     #Only get the fasta file out of the zip
#     unzip -p \
#         -j "${phaster}"/"${sample}"_phaster.zip \
#         "phage_regions.fna" \
#         > "${phaster}"/"${sample}"_phages.fasta

#     #Add sample name and entry number to fasta header
#     sed -i "s/^>/>"${sample}"_/" "${phaster}"/"${sample}"_phages.fasta
# done


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
echo "${sampleList[@]}" | tr " " "\n" >> "${phaster}"/sampleList.txt

#Convert CD-HIT-EST ".clstr" output file to OTU table
# Usage: perl cdHitClstr2table.pl <sampleList.txt> <cd-hit.clstr> <outputTable.tsv>
perl "${scripts}"/cdHitClstr2table.pl \
    "${phaster}"/sampleList.txt \
    "${phaster}"/phages_clustered.fasta.clstr \
    "${phaster}"/phages_clustered.tsv



