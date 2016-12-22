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

Reduced datasets:
* **Test case 1**: Coverage analysis of PanelApp panel in a reduced dataset (no whole genome stats)
* **Test case 2**: Coverage analysis of a gene list in a reduced dataset (no whole genome stats)

Whole genome datasets:
* **Test case 3**: Coverage analysis of whole exome in a cancer whole genome dataset
* **Test case 4**: Coverage analysis of PanelApp panel in a cancer whole genome dataset
* **Test case 5**: Coverage analysis of a gene list in a cancer whole genome dataset
* **Test case 6**: Coverage analysis of whole exome in a rare disease whole genome dataset
* **Test case 7**: Coverage analysis of PanelApp panel in a rare disease whole genome dataset
* **Test case 8**: Coverage analysis of a gene list in a rare disease whole genome dataset
* **Test case 9**: Coverage analysis of panel with a BED of NonN regions 
* **Test case 10**: Coverage analysis of whole exome with exon stats disabled
* **Test case 11**: Coverage analysis of whole genome aligned to GRCh38

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
/genomes/analysis/by_date/2016-09-27/HX01166491/CancerLP3000079-DNA_G03_NormalLP3000067-DNA_A07/coverage/LP3000079-DNA_G03.bw
/genomes/analysis/by_date/2016-09-27/HX01166465/CancerLP3000079-DNA_H03_NormalLP3000067-DNA_E05/coverage/LP3000079-DNA_H03.bw
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

# change chromosome idetifiers to non-prefixed
samtools view -H /genomes/by_date/2016-11-18/RAREP40001/LP2000274-DNA_B11/Assembly/LP2000274-DNA_B11.bam | sed  -e 's/SN:chrM/SN:MT/' -e 's/SN:chr/SN:/'  | samtools reheader - /genomes/by_date/2016-11-18/RAREP40001/LP2000274-DNA_B11/Assembly/LP2000274-DNA_B11.bam > /genomes/analysis/by_date/2016-11-18/RAREP40001/LP2000274-DNA_B11/coverage/LP2000274-DNA_B11.wochr.bam

# Extract chromosomes and their lengths from BAM header (adds chr prefix...)
/genomes/software/apps/python2.7/bin/python /genomes/software/apps/gel-coverage/bam2wig/get_chr_sizes.py --bam /genomes/analysis/by_date/2016-11-18/RAREP40001/LP2000274-DNA_B11/coverage/LP2000274-DNA_B11.wochr.bam --output /genomes/analysis/by_date/2016-11-18/RAREP40001/LP2000274-DNA_B11/coverage/LP2000274-DNA_B11.wochr.chr

# Creates the bigwig
/genomes/software/apps/jdk1.8.0_45/bin/java -jar ~/src/gel-coverage/bam2wig/gel-coverage-jar-with-dependencies.jar -bam /genomes/analysis/by_date/2016-11-18/RAREP40001/LP2000274-DNA_B11/coverage/LP2000274-DNA_B11.wochr.bam -stdout | /genomes/software/src/ucsc/wigToBigWig stdin /genomes/analysis/by_date/2016-11-18/RAREP40001/LP2000274-DNA_B11/coverage/LP2000274-DNA_B11.wochr.chr /genomes/analysis/by_date/2016-11-18/RAREP40001/LP2000274-DNA_B11/coverage/LP2000274-DNA_B11.wochr.bw
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

Run the whole exome coverage analysis on a cancer whole genome sample.

* Bigwig: /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw

Beware that the __output JSON is 4 GB__.

Status: **PASSED**
Execution time: 4h28m

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
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config  --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test3.json
```

After the job has finished you should find the following files:
```
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test3.log
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test3.json
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test3.json.bed
```

### Test 4

Run a panel coverage analysis  on a cancer whole genome sample.

* Bigwig: /genomes/analysis/by_date/2016-09-27/HX01166491/CancerLP3000079-DNA_G03_NormalLP3000067-DNA_A07/coverage/LP3000079-DNA_G03.bw
* Panel name: Familial colon cancer
* Panel version: 1.3

Status: **PASSED**
Execution time: 57m

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
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2016-09-27/HX01166491/CancerLP3000079-DNA_G03_NormalLP3000067-DNA_A07/coverage/LP3000079-DNA_G03.bw --panel "Familial colon cancer" --panel-version 1.3 --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test4.json
```

After the job has finished you should find the following files:
```
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test4.log
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test4.json
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test4.json.bed
```

Run the automated verifications on the output JSON as follows:
```
/genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/json_verifier.py --json /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test4.json
```

### Test 5

Run a gene list coverage analysis  on a cancer whole genome sample.

* Bigwig: /genomes/analysis/by_date/2016-09-27/HX01166465/CancerLP3000079-DNA_H03_NormalLP3000067-DNA_E05/coverage/LP3000079-DNA_H03.bw
* Gene list: ABL1,EVI1,MYC,APC,IL2,TNFAIP3,ABL2,EWSR1,MYCL1,ARHGEF12,JAK2,TP53,AKT1,FEV,MYCN,ATM,MAP2K4,TSC1,AKT2,FGFR1 

Status: **PASSED**
Execution time: 58m

Run the following job:
```
qsub /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test5.sh
```

where gel-coverage.test5.sh:
```
#!/bin/bash
#$ -N gel-coverage.test5
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test5.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2016-09-27/HX01166465/CancerLP3000079-DNA_H03_NormalLP3000067-DNA_E05/coverage/LP3000079-DNA_H03.bw --gene-list ABL1,EVI1,MYC,APC,IL2,TNFAIP3,ABL2,EWSR1,MYCL1,ARHGEF12,JAK2,TP53,AKT1,FEV,MYCN,ATM,MAP2K4,TSC1,AKT2,FGFR1 --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test5.json
```

After the job has finished you should find the following files:
```
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test5.log
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test5.json
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test5.json.bed
```

Run the automated verifications on the output JSON as follows:
```
/genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/json_verifier.py --json /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test5.json
```

### Test 6

Run the whole exome coverage analysis on a rare disease whole genome sample.

* Bigwig: /genomes/analysis/by_date/2016-02-15/CH00349553/LP2000873-DNA_B02/coverage/LP2000873-DNA_B02.bw

Beware that the __output JSON is 4.2 GB__.

Status: **PASSED**
Execution time: 4h24m

Run the following job:
```
qsub /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test6.sh
```

where gel-coverage.test6.sh:
```
#!/bin/bash
#$ -N gel-coverage.test6
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test6.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2016-02-15/CH00349553/LP2000873-DNA_B02/coverage/LP2000873-DNA_B02.bw --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config  --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test6.json
```

After the job has finished you should find the following files:
```
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test6.log
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test6.json
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test6.json.bed
```

### Test 7

Run a panel coverage analysis  on a rare disease whole genome sample.

* Bigwig: /genomes/analysis/by_date/2015-01-08/RAREP01497/LP2000274-DNA_B11/coverage/LP2000274-DNA_B11.bw
* Panel name: "Charcot-Marie-Tooth disease"
* Panel version: 1.1

Status: **FAILED (%>=15x needs further verification)** (uncovered genes chromosome tag)
Execution time: 56m

Run the following job:
```
qsub /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test7.sh
```

where gel-coverage.test7.sh:
```
#!/bin/bash
#$ -N gel-coverage.test7
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test7.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2016-09-27/HX01166491/CancerLP3000079-DNA_G03_NormalLP3000067-DNA_A07/coverage/LP3000079-DNA_G03.bw --panel "Charcot-Marie-Tooth disease" --panel-version 1.1 --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test7.json
```

After the job has finished you should find the following files:
```
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test7.log
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test7.json
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test7.json.bed
```

Run the automated verifications on the output JSON as follows:
```
/genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/json_verifier.py --json /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test7.json
```

Check that the percentage of bases >=15x at whole genome is of 95.96680. **It is of 89.534%**

### Test 8

Run a gene list coverage analysis on a rare disease whole genome sample.

* Bigwig: /genomes/analysis/by_date/2014-12-18/RAREP01384/LP2000274-DNA_D01/coverage/LP2000274-DNA_D01.bw
* Gene list: DIAPH1,KCNQ4,GJB3,GJB2,MYH14,DFNA5,WFS1,TECTA,COCH,EYA4,MYO7A,COL11A2,POU4F3,MYH9,ACTG1,MYO6,SIX1,SLC17A8,GRHL2,TMC1,DSPP,P2RX2,CCDC50,MYO1A,MIR96,TJP2
(Clinical Manifestations and Molecular Genetics of Known Genes Causing Autosomal Dominant Nonsyndromic Hearing Impairment from GeneReviews)

Status: **%>=15x needs further verification**
Execution time: 56m 

Run the following job:
```
qsub /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test8.sh
```

where gel-coverage.test4.sh:
```
#!/bin/bash
#$ -N gel-coverage.test8
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test8.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2014-12-18/RAREP01384/LP2000274-DNA_D01/coverage/LP2000274-DNA_D01.bw --gene-list DIAPH1,KCNQ4,GJB3,GJB2,MYH14,DFNA5,WFS1,TECTA,COCH,EYA4,MYO7A,COL11A2,POU4F3,MYH9,ACTG1,MYO6,SIX1,SLC17A8,GRHL2,TMC1,DSPP,P2RX2,CCDC50,MYO1A,MIR96,TJP2 --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test8.json
```

After the job has finished you should find the following files:
```
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test8.log
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test8.json
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test8.json.bed
```

Run the automated verifications on the output JSON as follows:
```
/genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/json_verifier.py --json /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test5.json
```

Check that the percentage of bases >=15x at whole genome is of 98.25317. **It is of 89.504%**

### Test 8.1

Same as 8 but computing the whole genome statistics only on non N regions.

* Bigwig: /genomes/analysis/by_date/2014-12-18/RAREP01384/LP2000274-DNA_D01/coverage/LP2000274-DNA_D01.bw
* Whole genome region: /home/pferreiro/src/gel-coverage/resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.bed

**Observed %>=15x**: 96.925%
**Expected %>=15x**: 98.25317%

### Test 8.2

Same as 8.1

* Bigwig: /genomes/analysis/by_date/2016-02-15/CH00349553/LP2000873-DNA_B02/coverage/LP2000873-DNA_B02.bw
* Whole genome region: /home/pferreiro/src/gel-coverage/resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.bed

**Observed percentage >= 15x**: 91.028%
**Expected percentage >= 15x**: 93.20164%

### Test 8.3

Same as 8.1

* Bigwig: /genomes/analysis/by_date/2015-01-08/RAREP01497/LP2000274-DNA_B11/coverage/LP2000274-DNA_B11.bw
* Whole genome region: /home/pferreiro/src/gel-coverage/resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.bed

**Observed percentage >= 15x**: 92.529%
**Expected percentage >= 15x**: 95.96680

### Test 8.4

Same as 8.1 

* Bigwig: /genomes/analysis/by_date/2015-03-20/RAREP01974/LP2000731-DNA_H10/coverage/LP2000731-DNA_H10.bw
* Whole genome region: /home/pferreiro/src/gel-coverage/resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.bed

**Observed percentage >= 15x**: 94.75%
**Expected percentage >= 15x**: 99.27863


### Test 9

Run a panel coverage analysis providing a bed file with the nonN regions for the whole genome statistics.

* Bigwig: /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw
* Whole genome region: /home/pferreiro/src/gel-coverage/resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.bed
* Panel name: Familial colon cancer
* Panel version: 1.3

Status: **PASSED**
Execution time: 51m

Run the following job:
```
qsub /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test9.sh
```

where gel-coverage.test9.sh:
```
#!/bin/bash
#$ -N gel-coverage.test9
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test9.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config  --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test9.json --wg-regions /home/pferreiro/src/gel-coverage/resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.bed --panel "Familial colon cancer" --panel-version 1.3
```

After the job has finished you should find the following files:
```
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test9.log
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test9.json
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test9.json.bed
```

Run the automated verifications on the output JSON as follows (you may want to enqueue this in SGE):
```
/genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/json_verifier.py --json /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test9.json
```

### Test 10

Run the whole exome coverage analysis providing a bed file with the nonN regions for the whole genome statistics and disabling the exon level statatistics to have a smaller JSON.

* Bigwig: /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw
* Whole genome region: /home/pferreiro/src/gel-coverage/resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.bed
* Exon stats disabled

Beware that the __output JSON is 262 MB__.

Status: **FAILED** (invalid chromosome tag in JSON)
Execution time: 4h28m

Run the following job:
```
qsub /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test10.sh
```

where gel-coverage.test10.sh:
```
#!/bin/bash
#$ -N gel-coverage.test10
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test10.log
#$ -j y
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/bigwig_analyser.py --bw /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw --config /home/pferreiro/src/gel-coverage/resources/exon_coverage_summary.config  --output /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test9.json --wg-regions /home/pferreiro/src/gel-coverage/resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.bed --disable-exon-stats
```

After the job has finished you should find the following files:
```
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test10.log
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test10.json
/home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test10.json.bed
```

Run the automated verifications on the output JSON as follows (you may want to enqueue this in SGE):
```
/genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/json_verifier.py --json /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test10.json
```

### Test 11

Run a panel coverage analysis  on a rare disease whole genome sample aligned to the GRCh38.

* Bigwig: /genomes/analysis/by_date/2016-11-18/RAREP40001/LP2000274-DNA_B11/coverage/LP2000274-DNA_B11.wochr.bw
* Panel name: "Charcot-Marie-Tooth disease"
* Panel version: 1.1

Status: **PASSED**
Execution time: 52 m

Run the following job:
```
qsub /home/pferreiro/src/gel-coverage/resources/test/gel-coverage.test11.sh
```

## Other test results

### Uneveness
Testing the uneveness results against the former implementation results.

```
Sample	Expected	Observed	Diff
LP2000904-DNA_E01	6.19	6.25	0.06
LP2000904-DNA_B01	6.57	6.647	0.077
LP2000906-DNA_C01	13.55	13.879	0.329
LP2000907-DNA_A01	23.2	23.615	0.415
LP2000906-DNA_B02	14.81	15.167	0.357
LP2000906-DNA_D02	14.02	14.383	0.363
LP2000907-DNA_D01	46.48	47.18	0.7
LP2000907-DNA_A02	48.6	49.182	0.582
LP2000907-DNA_G01	20.56	20.754	0.194
LP2000907-DNA_D02	46.4	47.166	0.766
LP2000904-DNA_C01	6.86	6.942	0.082
LP2000904-DNA_H01	6.54	6.622	0.082
LP2000907-DNA_B02	25.0	25.283	0.283
LP2000907-DNA_H01	50.7	51.236	0.536
LP2000906-DNA_F01	12.91	13.201	0.291
```

### Cosmic variants coverage
Testing the percentage of cosmic variants regions covered below 30x against the former implementation results (GATK's DepthOfCoverage).

```
Sample	Expected	Observed	100-(observed*100)	Diff
LP2000906-DNA_A01	1.15	0.98911	1.089	0.061
LP2000906-DNA_A02	1.04	0.98985	1.015	0.025
LP2000906-DNA_B02	1.26	0.98776	1.224	0.036
LP2000906-DNA_C01	1.14	0.98904	1.096	0.044
LP2000906-DNA_C02	1.05	0.98986	1.014	0.036
LP2000906-DNA_D01	1.05	0.98997	1.003	0.047
LP2000906-DNA_D02	0.95	0.99106	0.894	0.056
LP2000906-DNA_F01	0.94	0.9911	0.89	0.05
LP2000906-DNA_G01	7.62	0.92245	7.755	-0.135
LP2000906-DNA_H01	1.18	0.98832	1.168	0.012
LP2000907-DNA_A01	1.11	0.98956	1.044	0.066
LP2000907-DNA_A02	6.17	0.93863	6.137	0.033
LP2000907-DNA_B01	0.98	0.99066	0.934	0.046
LP2000907-DNA_B02	1.55	0.98485	1.515	0.035
LP2000907-DNA_D01	3.96	0.95982	2.515	1.035
LP2000907-DNA_D02	4.16	0.95929	4.071	0.089
LP2000907-DNA_E01	1.16	0.9886	1.14	0.02
LP2000907-DNA_E02	14.63	0.85423	14.577	0.053
LP2000907-DNA_F02	1.41	0.98634	1.366	0.044
LP2000907-DNA_G01	1.60	0.98402	1.598	0.002
LP2000907-DNA_G02	1.16	0.98886	1.114	0.046
LP2000907-DNA_H01	4.64	0.95467	4.533	0.107
```

### Percentage  >= 15x
Testing the percentage coverage at a depth equal or greater than 15x on some rare diseases samples against the former coverage.

```
Samples	Expected %>=15x	Observed %>=15x	Diff
LP2000873-DNA_B02	93.20164	91.028	-2.17364
LP2000274-DNA_B11	95.96680	92.529	-3.4378
LP2000274-DNA_D01	98.25317	96.925	-1.32817
LP2000731-DNA_H10	99.27863	94.75	-4.52863
```
