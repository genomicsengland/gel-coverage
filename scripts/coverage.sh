#!/bin/bash

# Specifies the interpreting shell for the job
#$ -S /bin/bash

# Run the job from the current working directory
#$ -cwd

# The queue to which the job is to be submitted
#$ -q all.q

# Run at most 4 jobs per node
# (use 2 of the 8 large_file_scan consumables)
#$ -l large_file_scan=2

export PATH=/genomes/software/apps/gel-coverage/scripts:${PATH}

module load python-tiering/2.7.10
module load java/jdk1.8.0_45

# Check number of arguments
if [[ $# -ne 1 ]]
  then
  echo "Exactly one argument is required, the path to a text file where each line lists the path corresponding to an array job, e.g. one line might be /genomes/by_date/2015-12-16/0000234403/LP2000950-DNA_A10"
  exit 1
fi

sge_array_file=$1
pth=$(awk "NR==$SGE_TASK_ID" ${sge_array_file})
lp=$(echo $pth | awk -F"/" '{print $6}')
bam=${pth}/Assembly/${lp}.bam

pth2=$(echo $pth | sed 's:/genomes/by_date/:/genomes/analysis/by_date/:')
direc=${pth2}/coverage/
bigwig=${direc}/${lp}.bw

# Exit if the bigwig file already exists
if [ -f "${bigwig}" ]
then
  exit 0
fi

# Create the directory structure if necessary
if [ ! -d "${direc}" ]
then
  mkdir -p $direc
fi

# 1. Make chromosome size file
python /genomes/software/apps/gel-coverage/bam2wig/get_chr_sizes.py \
  --bam $bam \
  --output ${direc}/${lp}.chr

# 2. Make bigwig file
java \
  -jar /genomes/software/apps/gel-coverage/bam2wig/gel-coverage-jar-with-dependencies.jar \
  -bam ${bam} \
  -stdout | \
    /genomes/software/src/ucsc/wigToBigWig stdin ${direc}/${lp}.chr $bigwig

# 3. Make coverage summary
python /genomes/software/apps/gel-coverage/scripts/coverage_summary.py \
  --bw $bigwig \
  --xlim 101 \
  --outdir $direc \
  --outprefix $lp \
  --genome_n /genomes/software/apps/gel-coverage/resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.bed
