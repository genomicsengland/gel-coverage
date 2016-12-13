bam=$1
bigwig=$2


java \
  -jar /genomes/software/apps/gel-coverage/bam2wig/gel-coverage-jar-with-dependencies.jar \
  -bam ${bam} \
  -output - | \
    /genomes/software/src/ucsc/wigToBigWig stdin $bam.chr $bigwig