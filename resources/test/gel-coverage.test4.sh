#!/bin/bash
#$ -N gel-coverage.test4
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test4.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/exon_coverage_summary.py --bw /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw --panel "Familial colon cancer" --panel-version 1.3 --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test4.json
