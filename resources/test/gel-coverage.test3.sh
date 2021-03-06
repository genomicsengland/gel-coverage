#!/bin/bash
#$ -N test3
#$ -q all.q
#$ -S /bin/bash
#$ -o /genomes/scratch/pferreiro/gelcoverage/gel-coverage.test3.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw --config /home/pferreiro/src/gel-coverage/resources/bigwig_analyser.config  --output /genomes/scratch/pferreiro/gelcoverage/gel-coverage.test3.json
