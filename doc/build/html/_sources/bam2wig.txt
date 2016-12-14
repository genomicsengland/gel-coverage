Bam2Wig
=======

This is a java program that calculate the coverage from a BAM file. Mode of use:

.. code-block:: bash

   module load java/jdk1.8.0_45
   java -jar /genomes/software/apps/gel-coverage/gel-coverage-jar-with-dependencies.jar -bam <BAM_FILE> -output <OUTPUT_PREFIX>

.. note::

   You will get <OUTPUT_PREFIX>.wig <OUTPUT_PREFIX>.chr

Due to the size of the output it is quite convenient produce a bigWig. In our pipeline should be run in the following way:

.. code-block:: bash
   module load java/jdk1.8.0_45
   module load python
   java -jar /genomes/software/apps/gel-coverage/gel-coverage-jar-with-dependencies.jar -bam <BAM_FILE>  -output <OUTPUT_PREFIX>
   wigToBigWig <OUTPUT_PREFIX>.wig <OUTPUT_PREFIX>.chr <OUTPUT_PREFIX>.bw





