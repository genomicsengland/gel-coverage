
# bam2wig

The objective is to create a WIG file with the depth of coverage data from a given BAM, after performing some basic
filtering on the reads.

## Read filtering

Some of the reads in the BAM are filtered out following certain criteria like base quality, mapping quality and
duplicate reads.
These parameters are provided from a configuration file with the following structure:
```
[bam2wig]
base_quality : 30
mapping_quality : 10
filter_duplicates : true
```

## Build

Java version: JDK 1.8.0
Dependencies are managed by Maven.
Run at gel-coverage/bam2wig:
```
mvn compile
mvn package
```

This will create the file `gel-coverage-jar-with-dependencies.jar` under the folder target.

## Running

To create the coverage WIG file from a BAM file run:
```
java -jar bam2wig-jar-with-dependencies.jar --bam $BAM --wig $WIG --config $CONFIG
```

If you want to create a BIGWIG file you can pipe the output as follows:
```
java \
  -jar bam2wig-jar-with-dependencies.jar \
  --bam $BAM \
  --wig -
  --config $CONFIG| \
    /genomes/software/src/ucsc/wigToBigWig stdin $BAM.chr $BIGWIG
```

where $BAM.chr is a file with two tab-separated columns with the chrommosome name and the length. This file is created
either as $BAM.chr or by replacing the wig extension by chr from $WIG.


