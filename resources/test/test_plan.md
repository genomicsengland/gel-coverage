# Testing GEL-coverage

## Test plan

The GEL_coverage application covers three scenarios and as such we have mainly three test cases:
1. Coverage analysis of all genes in a panel from PanelApp
2. Coverage analysis of genes in a provided gene list
3. Coverage analysis of all genes in the genome

Consider extreme cases and error control:
TBD

### Test case 1
Coverage analysis of all gene in a panel from PanelApp.
Panel name: `Epileptic encephalopathy`
Panel version: `0.2`
Gene confidence: `HighEvidence`

### Test case 2
Coverage analysis of a gene list.
Gene list: `BRCA1, BRCA2, CFTR, IGHE`

### Test case 3
Coverage analysis of whole exome. This will be a long running test....

## Test data

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