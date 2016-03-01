#!/usr/local/bin/Rscript
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
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="file to process", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

file = opt$file

##############set theme############################


theme_gel_proper <- function(base_size = 10, base_family="Calibri") {
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
    legend.key.size =   unit(1.2, "lines"),
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


####################do stuff##############################

dat<-read.table(file,header=T,sep="\t")
print(head(dat))
dat$V5=dat$cov/median(dat$cov)
gc25<-quantile(dat[,6],0.25)
gc75<-quantile(dat[,6],0.75)

y.cov=4*median(dat$cov,na.rm=T)

m25<-dat[dat$gc<=gc25 & dat$exon<=30,]
m75<-dat[dat$gc>=gc75 & dat$exon<=30,]
m50<-dat[dat$gc>gc25 & dat$gc<gc75  & dat$exon<=30,]
png(paste(file,"coverage.boxplots.png",sep="."),height=250,width=1000,type="cairo")
par(mfrow=c(1,3))

############do average coverage plots#################

xlabel="Average Coverage"

p1=ggplot(m25,aes(exon,cov,group=exon))+
  geom_boxplot(fill="olivedrab",outlier.size = 0.5,size = 0.3)+
  theme_gel_proper()+
  scale_x_continuous("Exons")+
  scale_y_continuous(xlabel)+
  ggtitle("Low GC") +
  coord_cartesian(ylim = c(0, y.cov))+
  geom_hline(yintercept=mean(dat$cov))

p2=ggplot(m50,aes(exon,cov,group=exon))+
  geom_boxplot(fill="goldenrod1",outlier.size = 0.5,size = 0.3)+
  theme_gel_proper()+
  scale_x_continuous("Exons")+
  scale_y_continuous(xlabel)+
  ggtitle("Moderate GC") +
  coord_cartesian(ylim = c(0, y.cov))+
  geom_hline(yintercept=mean(dat$cov))

p3=ggplot(m75,aes(exon,cov,group=exon))+
  geom_boxplot(fill="firebrick",outlier.size = 0.5,size = 0.3)+
  theme_gel_proper()+
  scale_x_continuous("Exons")+
  scale_y_continuous(xlabel)+
  ggtitle("High GC") +
  coord_cartesian(ylim = c(0, y.cov))+
  geom_hline(yintercept=mean(dat$cov))


#grid.arrange(p1, p2, p3, nrow=1,ncol=3)
g <- arrangeGrob(p1, p2, p3, nrow=1,ncol=3) #generates g
filename=paste(file,"average.coverage_gc.boxplots.png",sep=".")
filename=sub(".exon.coverage.means.with.GC.txt", "", filename, ignore.case =FALSE, fixed=FALSE)
ggsave(filename,g,width = 20, height = 10, units = "cm",dpi=300)

#####################do relative coverage plots######################

xlabel="Relative Average Coverage"

p1=ggplot(m25,aes(exon,V5,group=exon))+
  geom_boxplot(fill="olivedrab",outlier.size = 0.5,size = 0.3)+
  theme_gel_proper()+
  scale_x_continuous("Exons")+
  scale_y_continuous(xlabel)+
  ggtitle("Low GC") +
  coord_cartesian(ylim = c(0, 4))+
  geom_hline(yintercept=mean(dat$V5))


p2=ggplot(m50,aes(exon,V5,group=exon))+
  geom_boxplot(fill="goldenrod1",outlier.size = 0.5,size = 0.3)+
  theme_gel_proper()+
  scale_x_continuous("Exons")+
  scale_y_continuous(xlabel)+
  ggtitle("Moderate GC") +
  coord_cartesian(ylim = c(0, 4))+
  geom_hline(yintercept=mean(dat$V5))

p3=ggplot(m75,aes(exon,V5,group=exon))+
  geom_boxplot(fill="firebrick",outlier.size = 0.5,size = 0.3)+
  theme_gel_proper()+
  scale_x_continuous("Exons")+
  scale_y_continuous(xlabel)+
  ggtitle("High GC") +
  coord_cartesian(ylim = c(0, 4))+
  geom_hline(yintercept=mean(dat$V5))

#grid.arrange(p1, p2, p3, nrow=1,ncol=3)
g <- arrangeGrob(p1, p2, p3, nrow=1,ncol=3) #generates g
filename=paste(file,"relative.coverage_gc.boxplots.png",sep=".")
filename=sub(".exon.coverage.means.with.GC.txt", "", filename, ignore.case =FALSE, fixed=FALSE)
ggsave(filename,g,width = 20, height = 10, units = "cm",dpi=300)

