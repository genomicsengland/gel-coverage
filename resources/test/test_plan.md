# Testing GEL-coverage

## Test plan

The GEL_coverage application covers three scenarios:
1. Coverage analysis of all genes in a panel from PanelApp
2. Coverage analysis of genes in a provided gene list
3. Coverage analysis of all genes in the genome

Consider extreme cases and error control:
TBD

### Test case 1
Coverage analysis of PanelApp panel in a reduced dataset.

### Test case 2
Coverage analysis of a gene list in a reduced dataset.

### Test case 3
Coverage analysis of whole exome in a whole genome dataset.

### Test case 4
Coverage analysis of PanelApp panel in a whole genome dataset.

### Test case 5
Coverage analysis of a gene list in a whole genome dataset.

## Test data

### Small datasets

We need to create small datasets for testing. We will create a subset BAM for very specific intervals that contains the genes that we will be analysing as our tests are defined.

Reference genome: GRCh37
Intervals were obtained with a 10x zoom out from the UCSC Genome Browser in order to include the gene interest together with others that won't be analysed.

Intervals of interest for test case 1. Genes included in panel "Epileptic encephalopathy" v0.2:
```
SCN2A: 2:165,239,488-165,392,032
SPTAN1: 9:128,187,585-128,998,634
PLCB1: 20:4,745,389-12,271,778
SLC25A22: 11:755,109-833,698
SCN8A: 12:50,593,905-52,810,194
STXBP1: 9:127,250,330-128,054,669
PNKP: 19:49,831,130-49,897,959
```

Intervals of interest for test case 2. Gene list (BRCA1, BRCA2, CFTR and IGHE):
```
CFTR: 7:116,630,799-118,517,828
BRCA1: 17:43,045,629-43,125,483
BRCA2: 13:31,933,905-32,781,834
IGHE: 14:105,591,895-105,609,774
```

Create the BAMs as follows:
```
samtools view -b /genomes/by_date/2016-11-25/HX01579108/LP3000160-DNA_A02/Assembly/LP3000160-DNA_A02.bam 2:165,239,488-165,392,032 9:49,831,130-49,897,959 9:127,250,330-128,054,669 9:128,187,585-128,998,634 11:755,109-833,698 12:50,593,905-52,810,194 20:4,745,389-12,271,778  > ~/test1.bam
samtools index test1.bam
samtools view -b /genomes/by_date/2016-11-25/HX01579108/LP3000160-DNA_A02/Assembly/LP3000160-DNA_A02.bam 7:116,630,799-118,517,828 17:43,045,629-43,125,483 13:31,933,905-32,781,834 14:105,591,895-105,609,774 > ~/test2.bam
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

### Big datasets

The folllowing file contains the coverage information of a whole genome cancer sample:
```
/genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw
```

## Test description

### Test 1

Run the coverage analysis in a PanelApp panel on a reduced dataset.

* Panel name: `Epileptic encephalopathy`
* Panel version: `0.2`
* Gene confidence: `HighEvidence`

Run:
```
/genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/exon_coverage_summary.py --bw TBD --output test1.json --panel "Epileptic encephalopathy" --panel-version 0.2
```

### Test 2

Run the coverage analysis in a gene list on a reduced dataset.

* Gene list: `BRCA1, BRCA2, CFTR, IGHE`

Run:
```
/genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/exon_coverage_summary.py --bw TBD --output test1.json --gene-list BRCA1, BRCA2, CFTR, IGHE
```

### Test 3

Run the whole exome coverage analysis on a whole genome sample.

* Bigwig: /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw

Run the following in SGE:
```
#!/bin/bash
#$ -cwd -V
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/data
#$ -e /home/pferreiro/data
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/exon_coverage_summary.py --bw /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw --output LP3000079-DNA_F03.wg.json > LP3000079-DNA_F03.wg.time
```

### Test 4

Run a panel coverage analysis  on a whole genome sample.

* Bigwig: /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw
* Panel name: Familial colon cancer
* Panel version: 1.3

Run the following in SGE:
```
#!/bin/bash
#$ -cwd -V
#$ -q all.q
#$ -S /bin/bash
#$ -o /home/pferreiro/data
#$ -e /home/pferreiro/data
time /genomes/software/apps/python2.7-coverage_tests/bin/python /home/pferreiro/src/gel-coverage/scripts/exon_coverage_summary.py --bw /genomes/analysis/by_date/2016-09-27/HX01166477/CancerLP3000079-DNA_F03_NormalLP3000067-DNA_C12/coverage/LP3000079-DNA_F03.bw --panel "Familial colon cancer" --panel-version 1.3 --output LP3000079-DNA_F03.panel.json > LP3000079-DNA_F03.panel.time
```

### Test 5

TBD


