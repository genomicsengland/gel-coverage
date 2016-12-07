#!/bin/bash
#$ -N gel-coverage.test6
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test6.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2016-02-15/CH00349553/LP2000873-DNA_B02/coverage/LP2000873-DNA_B02.bw --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config  --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test6.json
