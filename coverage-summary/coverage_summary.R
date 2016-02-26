#parse the arguments

library("optparse")

option_list = list(
  make_option(c("-f", "--files"), type="character", default=NULL, 
              help="files to process", metavar="character"),
  make_option(c("-l", "--labels"), type="character", default=NULL, 
              help="labels for plotting", metavar="labels"),
  make_option(c("-c", "--coverage"), type="character", default=NULL, 
              help="labels for plotting", metavar="labels"),
  make_option(c("-s", "--scope"), type="character", default=NULL, 
              help="labels for plotting", metavar="labels")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$files)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

files=opt$files
labels=opt$labels
covs=opt$coverage
scope=opt$scope
print(labels)
print(files)
#files="/Volumes/SCRATCH/temp/countstest.wg.coverage.counts.txt,/Volumes/SCRATCH/temp/countstest2.wg.coverage.counts.txt"
#labels="LP20001-DNA_A01,LP20001-DNA_A02"
##++++++++++ define functions ++++++++++++++++++++++
myprod<-function(x){as.numeric(x[1])*as.numeric(x[2])}
mysq<-function(x){((as.numeric(x[1])-as.numeric(x[3]))^2)*as.numeric(x[2])}
mymsn <- function(x,p,name) {  
  cof<-qnorm(1-p/2)
  subtotal <- apply(x,1,myprod)
  total <- sum(subtotal) #sum(k(i))
  n <- sum(as.numeric(x[,2]))
  x$m<-total/n
  m<-x$m[1]
  this.sq<-sqrt(sum(apply(x,1,mysq))/n)
  out<-data.frame(sample=name,
                  n=round(n,digits=0),
                  mean=round(m,digits=1),                  
                  sd=round(this.sq,digits=1),
                  lower=round(m-cof*this.sq,digits=1),
                  upper=round(m+cof*this.sq,digits=1),
                  pct25=round(qnorm(0.25,mean=m,sd=this.sq),digits=1),
                  median=round(qnorm(0.5,mean=m,sd=this.sq),digits=1),
                  pct75=round(qnorm(0.75,mean=m,sd=this.sq),digits=1))
  out
}
my.rcumsum <- function(x){
  y=100-cumsum(x)
  y1=c(100,y)
  y1[1:length(y1)-1]
}
chrs<-c(1:22,"X","Y")


#################my new code
require(reshape)
require(ggplot2)
require(grid)
require(plyr)
library(RColorBrewer)
require(extrafont)

theme_gel_proper <- function(base_size = 12, base_family="Calibri") {
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

gel_colours=c("#0ead84","#44546b","#addce9","#27b7cc","#90c684","#d3922d")

files=strsplit(files,",")[[1]]
print(files)
labels=strsplit(labels,",")[[1]]
print(files)

for(i in 1:length(files)) {

  file=files[i]
  #print(file)
  wg<-read.table(as.character(file),header=T,sep="\t",check.names = FALSE)
  
  #write out coverage summary
  summary=mymsn(wg[2:covs,c("Coverage","Total")],0.05,labels[i])
  filename=paste(file,scope,"coverage.summary.table","txt",sep=".")
  filename=sub(".coverage.counts.txt", "", filename, ignore.case =FALSE, fixed=FALSE)
  write.table(summary,filename,quote=F,row.names=F,col.names=T,sep="\t")
  
  prop=as.data.frame(prop.table(as.matrix(wg), 2) )
  prop$Coverage=as.numeric(row.names(prop))-1
  m.wg = melt(prop,id=c("Coverage","Total"))
  m.wg$sample=labels[i]
  print(head(m.wg))
  ggplot(m.wg,aes(Coverage,as.numeric(value)*100))+
    geom_bar(stat="identity",color=gel_colours[1])+
    ggtitle(labels[i])+
    facet_wrap(~variable,scales="free")+
    scale_fill_manual("Sample",values=gel_colours)+
    scale_y_continuous("Percent Total")+
    scale_x_continuous("Coverage")+
    theme_gel_proper()+
    coord_cartesian(xlim = c(0,101), ylim = c(0,5))
  filename=paste(file,scope,"coverage.distribution.chr-by-chr.png",sep=".")
  print(filename)

  filename=sub(".coverage.counts.txt", "", filename, ignore.case =FALSE, fixed=FALSE)
  
  ggsave(filename, width = 30, height = 20, units = "cm",dpi=600)
  if (i > 1){
    final= rbind(old,m.wg)
    old=final
  }else{
    old=m.wg
  }
}
if (length(files) > 1){
  ggplot(final,aes(Coverage,as.numeric(value)*100,fill=sample))+
    geom_bar(stat="identity",position="dodge")+
    facet_wrap(~variable,scales="free")+
    scale_y_continuous("Percent Total")+
    scale_fill_manual("Sample",values=gel_colours)+
    scale_x_continuous("Coverage")+
    theme_gel_proper()+
    coord_cartesian(xlim = c(0,101), ylim = c(0,5))
  filename=paste("all",scope,".coverage.distribution.chr-by-chr.png",sep=".")
  filename=sub(".coverage.counts.txt", "", filename, ignore.case =FALSE, fixed=FALSE)
  ggsave(filename, width = 30, height = 20, units = "cm",dpi=600)
}