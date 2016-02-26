#!/usr/bin/Rscript
require(gridExtra)
require(reshape)
require(ggplot2)
require(grid)
require(plyr)
library(RColorBrewer)
require(extrafont)

###############get options##################

library("optparse")

option_list = list(
  make_option(c("-w", "--wgenome"), type="character", default=NULL,
              help="whole genome coverage summary file", metavar="character"),
  make_option(c("-e", "--exon"), type="character", default=NULL,
              help="exon coverage summary file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$wgenome)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
if (is.null(opt$exon)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

wgfile = opt$wgenome
exonfile = opt$exon

##############set theme############################
my.rcumsum <- function(x){
  y=100-cumsum(x)
  y1=c(100,y)
  y1[1:length(y1)-1]
}

theme_gel_proper <- function(base_size = 6, base_family="Calibri") {
  theme(
    text = element_text(size=base_size,family=base_family),
    axis.line =         element_blank(),
    axis.text.x =       element_text(size=base_size-2,lineheight = 0.9, vjust = 1),
    axis.text.y =       element_text(size=base_size-2,lineheight = 0.9, hjust = 1),
    axis.ticks =        element_line(colour = "black", size = 0.2),
    axis.title.x =      element_text(size = base_size, vjust = 1),
    axis.title.y =      element_text(size = base_size, angle = 90, vjust = 0.5),
    axis.ticks.length = unit(0.3, "lines"),
    axis.ticks.margin = unit(0.5, "lines"),
    
    legend.background = element_rect(colour=NA),
    legend.key =        element_blank(),
    legend.key.size =   unit(0.8, "lines"),
    legend.title =      element_text(size = base_size,face = "bold", hjust = 0),
    legend.position =   "right",
    
    panel.background =  element_rect(fill = "white", colour = NA),
    panel.border = element_rect(color="black",size=0.5,fill="NA"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.margin =      unit(0.5, "lines"),
    
    strip.text.x =      element_text(),
    strip.text.y =      element_text(angle = -90),
    strip.background = element_rect(colour="black",size=0.6, fill="#FFFFFF"),
    
    plot.background =   element_rect(colour = NA),
    plot.title =        element_text(size = base_size * 1.2),
    plot.margin =       unit(c(1, 1, 0.5, 0.5), "lines")
    
    
    
    
  )
}

gel_colours=c("#0ead84","#44546b","#addce9","#27b7cc","#90c684","#d3922d")
####################do stuff##############################

wgfile="/Volumes/SCRATCH/temp/countstest.wg.coverage.counts.txt.exon.coverage.counts.txt"
exonfile="/Volumes/SCRATCH/temp/countstest.wg.coverage.counts.txt.exon.coverage.counts.txt"
covs=101
wg.g<-read.table(wgfile,header=T,sep="\t",check.names = FALSE)
exon.g<-read.table(exonfile,header=T,sep="\t")
xmin<-min(dim(exon.g)[[1]])
exon.g.total <- cbind(exon.g$Total[1:xmin])
exon.g.prop <- as.data.frame(prop.table(exon.g.total/100,2))
colnames(exon.g.prop)<-c("ex.g")
rownames(exon.g.prop)<-exon.g$Coverage[1:xmin]
exon.g.prop$germline<-cumsum(exon.g.prop[,1])

### combine genomic and exonic coverage ####
we.gd.total<-cbind(wg.g$Total[1:covs],exon.g$Total[1:covs])
we.gd.prop<-as.data.frame(100*prop.table(we.gd.total/100,2))
colnames(we.gd.prop)<-c("wg.g","ex.g")
we.gd.prop$wgg.cum<-my.rcumsum(we.gd.prop$wg.g)
we.gd.prop$exg.cum<-my.rcumsum(we.gd.prop$ex.g)
we.gd.prop$cov<- seq(0,(covs-1),1)
png(paste(wgfile,"cumulative_plot","png",sep="."),width=400,height=400,type="cairo")
plot(wgg.cum ~ cov,data=we.gd.prop, xlab="Coverage",ylab="Percent of Bases (%)",
     main=paste(wgfile," cumulative coverage",sep=""),col="green",pch=22)
points(cum ~ cov,data=we.gd.prop,col="green",pch=19)
legend("topright",c("Genomic Normal","Genomic Disease","Exonic Normal","Exonic Disease"),
       pch=c(22,22,19,19),col=c("green","red"),bty="n",cex=0.9)
dev.off()

names(we.gd.prop)=c("Whole Genome","Exome","Whole Genome Cum.","Exon Cum.","cov")
m.we.gd.prop = melt(we.gd.prop[,c(3,4,5)],id="cov")

ggplot(m.we.gd.prop,aes(cov,value,color=variable))+
  theme_gel_proper()+
  geom_point()+
  ggtitle("Cumulative Coverage")+
  scale_y_continuous("Percent of Total")+
  scale_color_manual("Region",values=gel_colours)+
  scale_x_continuous("Coverage")

filename=paste(wgfile,"cumulative coverage.png",sep=".")
#filename=sub(".coverage.counts.txt", "", filename, ignore.case =FALSE, fixed=FALSE)
ggsave(filename, width = 9, height = 5, units = "cm",dpi=300)

names(we.gd.prop)=c("Whole Genome","Exome","Whole Genome Cum.","Exon Cum.","cov")
m.we.gd.prop = melt(we.gd.prop[,c(1,2,5)],id="cov")
ggplot(m.we.gd.prop,aes(cov,value,fill=variable))+
  geom_bar(stat="identity",position="dodge")+
  ggtitle("All Chromosomes Coverage Distribution")+
  scale_y_continuous("Percent Total")+
  scale_fill_manual("Region",values=gel_colours)+
  scale_x_continuous("Coverage")+
  theme_gel_proper()+
  coord_cartesian(xlim = c(0,101), ylim = c(0,5))
  
filename=paste(wgfile,"all.coverage.distribution.all.chrs.png",sep=".")
#filename=sub(".coverage.counts.txt", "", filename, ignore.case =FALSE, fixed=FALSE)
ggsave(filename, width = 9, height = 5, units = "cm",dpi=300)

png(paste(wgfile,"coverage.distribution.all.chrs","png",sep="."),height=1200,width=1600,type="cairo")
par(mfrow=c(2,1))
barplot(as.matrix(t(we.gd.prop[,c(1,2)])),names.arg=wg.g$Coverage[1:covs],ylim=c(0,5),
        main="Normal Coverage",xlab="Coverage",ylab="Percent of Total (%)",
        axis.lty=1,cex.names=0.7,legend.text=c("Genomic","Exonic"),beside=T,col=c("green",rgb(0,30/255,0)))

barplot(as.matrix(t(we.gd.prop[,c(2,4)])),names.arg=exon.g$Coverage[1:covs],ylim=c(0,5),
        main="Disease Coverage",xlab="Coverage",ylab="Percent of Total (%)",
        axis.lty=1,cex.names=0.7,legend.text=c("Genomic","Exonic"),beside=T,col=c("red",rgb(30/255,0,0)))
dev.off()