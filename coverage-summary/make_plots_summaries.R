#!/usr/bin/env Rscript
require(gridExtra)
require(reshape)
require(ggplot2)
require(grid)
require(plyr)
library(RColorBrewer)
require(extrafont)
library(gelCoverageR)
library("optparse")

option_list = list(
  make_option(c("-w", "--wgs"), type="character", default=NULL,
              help="files to process", metavar="character"),
  make_option(c("-e", "--exon"), type="character", default=NULL,
              help="files to process", metavar="character"),
  make_option(c("-g", "--gc"), type="character", default=NULL,
              help="files to process", metavar="character"),
  make_option(c("-l", "--label"), type="character", default=NULL,
              help="labels for plotting", metavar="labels"),
  make_option(c("-c", "--coverage"), type="character", default=NULL,
              help="xlim for plots", metavar="labels")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

wgs=opt$wgs
exon=opt$exon
gc=opt$gc
label=opt$label
covs= as.integer(opt$coverage)

print(wgs)
print(exon)
print(gc)
print(label)
print(covs)

coverage_summary(wgs,label,covs,"wg")

coverage_summary(exon,label,covs,"exon")

gc_exon_boxplots(gc)

cumulative_coverage(label,wgs,exon,covs)


