
# bam2wig

The objective is to create a WIG file containing the coverage information from a gievn BAM. 

## Read filtering
Some of the reads in the BAM are filtered out following certain criteria like base quality, mapping quality and duplicate reads. The parameters applied by default are:

* Minimum base quality = 30
* Minimum mapping quality = 10
* Filter out duplicate reads = True

## Build
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
java -jar gel-coverage-jar-with-dependencies.jar -bam $BAM -output $WIG
```

If you want to create a BIGWIG file you can pipe the output as follows:
```
java \
  -jar gel-coverage-jar-with-dependencies.jar \
  -bam $BAM \
  -stdout | \
    /genomes/software/src/ucsc/wigToBigWig stdin $CHR $BIGWIG
```

where $CHR is a file with two tab-separated columns with the chrommosome name and the length.


