bam=$1
bigwig=$2


/genomes/software/apps/jdk1.8.0_45/bin/java \
  -jar /home/pferreiro/src/gel-coverage//bam2wig/gel-coverage-jar-with-dependencies.jar \
  --bam ${bam} \
  --output - \
  --config /home/pferreiro/src/gel-coverage/resources/bigwig_analyser.config | \
/genomes/software/src/ucsc/wigToBigWig stdin $bam.chr $bigwig