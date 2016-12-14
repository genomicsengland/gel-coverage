bam=$1
bigwig=$2
jar=$3
config=$4

module load java/jdk1.8.0_45
java \
  -jar ${jar} \
  --bam ${bam} \
  --wig - \
  --config ${config} | \
/genomes/software/src/ucsc/wigToBigWig stdin ${bam}.chr ${bigwig}
