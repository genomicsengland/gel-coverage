#!/bin/bash
#$ -N test11-cvg
#$ -q all.q
#$ -S /bin/bash
#$ -o /genomes/scratch/pferreiro/gelcoverage/gel-coverage.test11.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2016-11-18/RAREP40001/LP2000274-DNA_B11/coverage/LP2000274-DNA_B11.wochr.bw --panel "Charcot-Marie-Tooth disease" --panel-version 1.1 --config /home/pferreiro/src/gel-coverage/resources/bigwig_analyser.GRCh38.config --output /genomes/scratch/pferreiro/gelcoverage/gel-coverage.test11.json
