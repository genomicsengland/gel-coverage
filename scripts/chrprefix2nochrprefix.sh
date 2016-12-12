
BIGWIG=$1

/genomes/software/src/ucsc/bigWigToWig $BIGWIG stdout | sed -e 's/chr//g' | /genomes/software/src/ucsc/bigWigToWig stdin $BIGWIG.wochr.bw