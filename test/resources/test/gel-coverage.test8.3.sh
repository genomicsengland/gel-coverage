#!/bin/bash
#$ -N test8.3-cvg
#$ -q all.q
#$ -S /bin/bash
#$ -o /genomes/scratch/pferreiro/gelcoverage/gel-coverage.test8.3.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2015-01-08/RAREP01497/LP2000274-DNA_B11/coverage/LP2000274-DNA_B11.bw --gene-list DIAPH1,KCNQ4,GJB3,GJB2,MYH14,DFNA5,WFS1,TECTA,COCH,EYA4,MYO7A,COL11A2,POU4F3,MYH9,ACTG1,MYO6,SIX1,SLC17A8,GRHL2,TMC1,DSPP,P2RX2,CCDC50,MYO1A,MIR96,TJP2 --config /home/pferreiro/src/gel-coverage/resources/bigwig_analyser.config --output /genomes/scratch/pferreiro/gelcoverage/gel-coverage.test8.3.json --wg-regions /home/pferreiro/src/gel-coverage/resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.prefix.bed
