#!/bin/bash
#$ -N gel-coverage.test7
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test7.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2016-09-27/HX01166491/CancerLP3000079-DNA_G03_NormalLP3000067-DNA_A07/coverage/LP3000079-DNA_G03.bw --panel "Charcot-Marie-Tooth disease" --panel-version 1.1 --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test7.json
