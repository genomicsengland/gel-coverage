bam2wig
=======

This is a java program that creates a coverage WIG file from a BAM file, filtering some of the reads based on base quality, mapping quality and duplicates.

Compile and build:

.. code-block:: bash

   mvn compile
   mvn package


Mode of use:

.. code-block:: bash

   module load java/jdk1.8.0_45
   java -jar /genomes/software/apps/gel-coverage/target/bam2wig-jar-with-dependencies.jar --bam <BAM_FILE> --wig <WIG_FILE> --config <CONFIG>


Where the config file has the following content:

.. code-block:: bash

   [bam2wig]
   base_quality : 30
   mapping_quality : 10
   filter_duplicates : true

.. note::

   You will get a wig file and a file with *.chr extension.

Due to the size of the output it is quite convenient produce a bigWig. In our pipeline should be run in the following way:

.. code-block:: bash

   module load java/jdk1.8.0_45
   module load python
   java -jar /genomes/software/apps/gel-coverage/bam2wig-jar-with-dependencies.jar --bam <BAM_FILE> --wig <WIG_FILE>
   wigToBigWig <WIG_FILE> <CHR_FILE> <BIGWIG_FILE>





