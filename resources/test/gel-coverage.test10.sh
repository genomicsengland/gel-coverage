#!/bin/bash
#$ -N gel-coverage.test10
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test10.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config  --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test10.json --wg-regions /home/pferreiro/src/gel-coverage/resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.bed --disable-exon-stats
