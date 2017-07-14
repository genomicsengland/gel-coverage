#!/bin/bash
#$ -N gel-coverage.test3.verification
#$ -q all.q
#$ -S /bin/bash
#$ -o /genomes/scratch/pferreiro/gelcoverage/gel-coverage.test3.verification.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/json_verifier.py --json /genomes/scratch/pferreiro/gelcoverage/gel-coverage.test3.json
