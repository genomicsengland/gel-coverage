#!/bin/bash
#$ -N gel-coverage.test5
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test5.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2016-09-27/HX01166465/CancerLP3000079-DNA_H03_NormalLP3000067-DNA_E05/coverage/LP3000079-DNA_H03.bw --gene-list ABL1,EVI1,MYC,APC,IL2,TNFAIP3,ABL2,EWSR1,MYCL1,ARHGEF12,JAK2,TP53,AKT1,FEV,MYCN,ATM,MAP2K4,TSC1,AKT2,FGFR1 --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test5.json
