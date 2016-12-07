# Testing GEL-coverage

## Test plan

The GEL_coverage application covers three scenarios:
1. Coverage analysis of all genes in a panel from PanelApp
2. Coverage analysis of genes in a provided gene list
3. Coverage analysis of all genes in the genome

Consider extreme cases and error control:
TBD

Test cases 1 and 2 are not time consuming and they are used for unit testing.
Test cases 3,4,5,6,7 and 8 are using whole genome datasets and they are not intended for unit testing.

* **Test case 1**: Coverage analysis of PanelApp panel in a reduced dataset (no whole genome stats)
* **Test case 2**: Coverage analysis of a gene list in a reduced dataset (no whole genome stats)
* **Test case 3**: Coverage analysis of whole exome in a cancer whole genome dataset
* **Test case 4**: Coverage analysis of PanelApp panel in a cancer whole genome dataset
* **Test case 5**: Coverage analysis of a gene list in a cancer whole genome dataset
* **Test case 6**: Coverage analysis of whole exome in a rare disease whole genome dataset
* **Test case 7**: Coverage analysis of PanelApp panel in a rare disease whole genome dataset
* **Test case 8**: Coverage analysis of a gene list in a rare disease whole genome dataset

Parameters are tested on unit tests.

## Test data

### Small datasets

We need to create small datasets for testing. We will create a subset BAM for very specific intervals that contains the genes that we will be analysing as our tests are defined.

Reference genome: GRCh37
Intervals were obtained adding 100k bp up and downstream of gene interval as reported by Ensembl.

Intervals of interest for test case 1. Genes included in panel "Epileptic encephalopathy" v0.2:
```
SCN2A: 2: 165,995,882-166,349,242
SPTAN1: 9:131,214,837-131,495,941
PLCB1: 20:8,012,824-9,049,003
SLC25A22: 11:690,475-898,333
SCN8A: 12:51,884,050-52,306,648
STXBP1: 9:130,341,649-130,458,394
PNKP: 19:50,263,139-50,481,608
```

Intervals of interest for test case 2. Gene list (BRCA1, BRCA2, CFTR and IGHE):
```
CFTR: 7:117,005,838-117,456,025
BRCA1: 17:41,096,312-41,422,262
BRCA2: 13:32,789,611-33,074,403
IGHE: 14:105,964,224-106,168,065
```

Create the BAMs as follows:
```
samtools view -b /genomes/by_date/2016-11-25/HX01579108/LP3000160-DNA_A02/Assembly/LP3000160-DNA_A02.bam SCN2A: 2:165,995,882-166,349,242 9:131,214,837-131,495,941 20:8,012,824-9,049,003 11:690,475-898,333 12:51,884,050-52,306,648 9:130,341,649-130,458,394 19:50,263,139-50,481,608  > ~/test1.bam
samtools index test1.bam
samtools view -b /genomes/by_date/2016-11-25/HX01579108/LP3000160-DNA_A02/Assembly/LP3000160-DNA_A02.bam CFTR: 7:117,005,838-117,456,025 17:41,096,312-41,422,262 13:32,789,611-33,074,403 14:105,964,224-106,168,065 > ~/test2.bam
samtools index test1.bam
```

We will need a file describing chromosome lengths:
```
/genomes/software/apps/python2.7/bin/python /genomes/software/apps/gel-coverage/bam2wig/get_chr_sizes.py --bam ~/data/test1.bam --output ~/data/test1.chr
/genomes/software/apps/python2.7/bin/python /genomes/software/apps/gel-coverage/bam2wig/get_chr_sizes.py --bam ~/data/test2.bam --output ~/data/test2.chr
```

Create the bigwigs as follows:
```
/genomes/software/apps/jdk1.8.0_45/bin/java -jar ~/src/gel-coverage/bam2wig/gel-coverage-jar-with-dependencies.jar -bam ~/data/test1.bam -stdout | /genomes/software/src/ucsc/wigToBigWig stdin ~/data/test1.chr ~/data/test1.bw
/genomes/software/apps/jdk1.8.0_45/bin/java -jar ~/src/gel-coverage/bam2wig/gel-coverage-jar-with-dependencies.jar -bam ~/data/test2.bam -stdout | /genomes/software/src/ucsc/wigToBigWig stdin ~/data/test1.chr ~/data/test2.bw
```

**NOTE:** tests on reduced datasets need to run using the flag --disable-wg-stats so running time is limited

### Big datasets

#### Cancer datasets
The folllowing files contains the coverage information of cancer whole genome samples:
```
/genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw
LP3000079-DNA_G03
LP3000079-DNA_H03
LP3000079-DNA_A04
LP3000079-DNA_B04
LP3000079-DNA_D04
LP3000079-DNA_C05
LP3000079-DNA_E05
LP3000079-DNA_G05
LP3000079-DNA_B06
```

#### Rare diseases datasets
The folllowing files contains the coverage information of rare disease whole genome samples.

Reference GRCh37:
```
/genomes/analysis/by_date/2016-02-15/CH00349553/LP2000873-DNA_B02/coverage/LP2000873-DNA_B02.bw
/genomes/analysis/by_date/2015-01-08/RAREP01497/LP2000274-DNA_B11/coverage/LP2000274-DNA_B11.bw
/genomes/analysis/by_date/2014-12-18/RAREP01384/LP2000274-DNA_D01/coverage/LP2000274-DNA_D01.bw
/genomes/analysis/by_date/2015-03-20/RAREP01974/LP2000731-DNA_H10/coverage/LP2000731-DNA_H10.bw
```

Expected results for percentage of bases >= 15x for the previous samples on GRCh37:
```
93.20164
95.96680
98.25317
99.27863
```

Reference GRCh38:
```
# we don't have bigwigs for this... generating them...
/genomes/by_date/2016-11-18/RAREP40001/LP2000274-DNA_B11
/genomes/by_date/2016-11-23/RAREP40046/LP2000861-DNA_G06
```

Expected results for percentage of bases >= 15x for the previous samples on GRCh38:
```
94.3
96.2
```

## Test descriptions

All tests are run with the following configuration unless explicitly stated otherwise:
```
[cellbase]
species: hsapiens
version: latest
assembly: GRCh37
host: 10.5.8.201:8080/cellbase-4.5.0-rc

[panelapp]
host: bioinfo.extge.co.uk/crowdsourcing/WebServices
# Comma separated list of gene confidence association to disease for genes in panels (possible values: HighEvidence, MediumEvidence, LowEvidence)
gene_confidence : HighEvidence

[transcript_filtering]
# Comma separated list determining those transcript flags that will be considered for analysis. Any of this
flags: basic
# Comma separated list determining those transcript biotypes that will be considered for analysis. Any of this
biotypes: IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,IG_V_gene,protein_coding,nonsense_mediated_decay,non_stop_decay,TR_C_gene,TR_D_gene,TR_J_gene,TR_V_gene

[exon]
exon_padding: 15
```

### Test 1

Run the coverage analysis in a PanelApp panel on a reduced dataset.

* Panel name: `Epileptic encephalopathy`
* Panel version: `0.2`
* Gene confidence: `HighEvidence` (this is set in the config file)
* Disable whole genome stats: True

Run:
```
/genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/exon_coverage_summary.py --bw TBD --output test1.json --panel "Epileptic encephalopathy" --panel-version 0.2 --disable-wg-stats
```

### Test 2

Run the coverage analysis in a gene list on a reduced dataset.

* Gene list: `BRCA1, BRCA2, CFTR, IGHE`
* Disable whole genome stats: True

Run:
```
/genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/exon_coverage_summary.py --bw TBD --output test1.json --gene-list BRCA1, BRCA2, CFTR, IGHE  --disable-wg-stats
```

### Test 3

Run the whole exome coverage analysis on a whole genome sample.

* Bigwig: /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw

Run the following job:
```
qsub /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test3.sh
```

where gel-coverage.test3.sh:
```
#!/bin/bash
#$ -N gel-coverage.test3
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test3.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/exon_coverage_summary.py --bw /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config  --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test3.json
```

### Test 4

Run a panel coverage analysis  on a whole genome sample.

* Bigwig: /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw
* Panel name: Familial colon cancer
* Panel version: 1.3

Run the following job:
```
qsub /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test4.sh
```

where gel-coverage.test4.sh:
```
#!/bin/bash
#$ -N gel-coverage.test4
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test4.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/exon_coverage_summary.py --bw /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw --panel "Familial colon cancer" --panel-version 1.3 --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test4.json
```

### Test 5

TBD


