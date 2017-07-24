my.rcumsum <- function(x){
    y=100-cumsum(x)
    y1=c(100,y)
    y1[1:length(y1)-1]
}

#' A function to produce the summary table
#'
#'
#' @param base_size base font size, default 12
#' @param base_family base font family, default Calibri
#' @keywords summary
#' @export
#' @examples
#' summary_table()
summary_table <- function(freq_table,well_id,column){

    mean=wtd.mean(freq_table[,"Coverage"], freq_table[,column])
    pct25=wtd.quantile(freq_table[,"Coverage"], freq_table[,column], probs = 0.25)
    median=wtd.quantile(freq_table[,"Coverage"], freq_table[,column], probs = 0.5)
    pct75=wtd.quantile(freq_table[,"Coverage"], freq_table[,column], probs = 0.75)
    var=wtd.var(freq_table[,"Coverage"], freq_table[,column])
    sd=sqrt(var)
    n=sum(freq_table[,column])

    out=data.frame(scope=column,
                    wellId=well_id,
                    n=round(n,digits=0),
                    mean=round(mean,digits=1),
                    sd=round(sd,digits=1),
                    pct25=round(pct25,digits=1),
                    median=round(median,digits=1),
                    pct75=round(pct75,digits=1))
    out

}

chrs<-c(1:22,"X","Y")

#' A function to gel a plot
#'
#'
#' @param base_size base font size, default 12
#' @param base_family base font family, default Calibri
#' @keywords theme
#' @export
#' @examples
#' theme_gel_proper()
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
    #axis.ticks.margin = unit(0.5, "lines"), DEPRECIATED

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


#' A function to summarise coverage
#'
#'
#' @param files a single file
#' @param labels a single label
#' @param covs ylim
#' @param scope wg or exome
#' @keywords summary
#' @export
#' @examples
#' coverage_summary()
coverage_summary <- function(file,label,covs,scope){

    #################make plots for all files###########################

    wg<-read.table(as.character(file),header=T,sep="\t",check.names = FALSE)

    allchr=summary_table(wg,label,"allchrs")
    autosomes=summary_table(wg,label,"autosomes")
    x=summary_table(wg,label,"X")
    y=summary_table(wg,label,"Y")

    summary_new=rbind(allchr,autosomes,x,y)

    filename_new=paste(file,"coverage.summary.table.txt",sep=".")
    filename_new=sub(".coverage.counts.txt", "", filename_new, ignore.case =FALSE, fixed=FALSE)
    write.table(summary_new,filename_new,quote=F,row.names=F,col.names=T,sep="\t")

    prop=as.data.frame(prop.table(as.matrix(wg), 2) )
    prop$Coverage=as.numeric(row.names(prop))-1
    print(head(prop))
    m.wg = melt(prop,id=c("Coverage","allchrs","autosomes"))
    m.wg$sample=label
    ymax=max(subset(m.wg,Coverage!=0)$value*100)
    print(head(m.wg))
    print(covs)
    print("plotting...")
    ggplot(m.wg,aes(Coverage,as.numeric(value)*100))+
        geom_bar(stat="identity",color=gel_colours[1])+
        ggtitle(label)+
        facet_wrap(~variable,scales="free")+
        scale_y_continuous("Percent Total",limits=c(0,ymax))+
        scale_fill_manual("Sample",values=gel_colours)+
        scale_x_continuous("Coverage",limits=c(0,covs))+
        theme_gel_proper()
    filename=paste(file,"coverage.distribution.chr-by-chr.png",sep=".")
    filename=sub(".coverage.counts.txt", "", filename, ignore.case =FALSE, fixed=FALSE)
    print(filename)
    ggsave(filename, width = 30, height = 20, units = "cm",dpi=600)

}


#' Plots multiple samples on same graphs
#'
#'
#' @param files a single files
#' @param labels a single label
#' @param covs ylim
#' @param scope wg or exome
#' @keywords summary
#' @export
#' @examples
#' multiple_sample_plots()
multiple_sample_plots <- function(files,labels,covs,scope){
  outfile_prefix=gsub(",", "_",labels)
  print(outfile_prefix)
  files=strsplit(files,",")[[1]]
  labels=strsplit(labels,",")[[1]]

  for(i in 1:length(files)) {

    file=files[i]
    wg<-read.table(as.character(file),header=T,sep="\t",check.names = FALSE)
    prop=as.data.frame(prop.table(as.matrix(wg), 2) )
    prop$Coverage=as.numeric(row.names(prop))-1
    m.wg = melt(prop,id=c("Coverage","allchrs"))
    m.wg$sample=labels[i]
    ymax=max(subset(m.wg,Coverage!=0)$value*100)

    if (i > 1){
      final= rbind(old,m.wg)
      old=final
    }else{
      old=m.wg
    }

  }


  ymax=max(subset(m.wg,Coverage!=0)$value*100)
  ggplot(final,aes(Coverage,as.numeric(value)*100,fill=sample))+
    geom_bar(stat="identity",position="dodge")+
    facet_wrap(~variable,scales="free")+
    scale_y_continuous("Percent Total",limits=c(0,ymax))+
    scale_fill_manual("Sample",values=gel_colours)+
    scale_x_continuous("Coverage",limits=c(0,covs))+
    theme_gel_proper()
  filename=paste(outfile_prefix,scope,"coverage.distribution.chr-by-chr.png",sep=".")
  filename=sub(".coverage.counts.txt", "", filename, ignore.case =FALSE, fixed=FALSE)
  print(filename)
  ggsave(filename, width = 30, height = 20, units = "cm",dpi=600)

}

#' A function to draw boxplots of coverage split by exon in transcript and gc content
#'
#'
#' @param file a file of exon coverage and gc cotent
#' @keywords boxplots
#' @export
#' @examples
#' gc_exon_boxplots()
gc_exon_boxplots <- function(file){

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


}

#' A function to plot cumulative coverage
#'
#'
#' @param label a single label i.e. well_id
#' @param wgfile the whole genome counts file
#' @param exonfile the exon counts file
#' @param covs ylim
#' @keywords cumulative
#' @export
#' @examples
#' cumulative_coverage()
cumulative_coverage <- function(label,wgfile,exonfile,covs){

  ####################do stuff##############################

  wg.g<-read.table(wgfile,header=T,sep="\t",check.names = FALSE)
  exon.g<-read.table(exonfile,header=T,sep="\t",check.names = FALSE)
  xmin<-min(dim(exon.g)[[1]])
  exon.g.total <- cbind(exon.g$allchrs[1:xmin])

  exon.g.prop <- as.data.frame(prop.table(exon.g.total/100,2))
  print("hello")
  colnames(exon.g.prop)<-c("ex.g")
  rownames(exon.g.prop)<-exon.g$Coverage[1:xmin]
  exon.g.prop$germline<-cumsum(exon.g.prop[,1])

  ### combine genomic and exonic coverage ####
  wg=wg.g$allchrs[1:covs]
  exon=exon.g$allchrs[1:covs]
  exon[is.na(exon)] <- 0
  wg[is.na(wg)] <- 0
  we.gd.total<-cbind(wg,exon)
  we.gd.prop<-as.data.frame(100*prop.table(we.gd.total/100,2))
  colnames(we.gd.prop)<-c("wg.g","ex.g")
  we.gd.prop$wgg.cum<-my.rcumsum(we.gd.prop$wg.g)
  we.gd.prop$exg.cum<-my.rcumsum(we.gd.prop$ex.g)
  we.gd.prop$cov<- seq(0,(covs-1),1)

  names(we.gd.prop)=c("Whole Genome","Exome","Whole Genome Cum.","Exon Cum.","cov")
  m.we.gd.prop = melt(we.gd.prop[,c(3,4,5)],id="cov")
  print(head(m.we.gd.prop))

  ggplot(m.we.gd.prop,aes(cov,value,color=variable))+
    theme_gel_proper()+
    geom_point()+
    ggtitle("Cumulative Coverage")+
    scale_y_continuous("Percent of Total")+
    scale_color_manual("Region",values=gel_colours)+
    scale_x_continuous("Coverage")+
    coord_cartesian(xlim = c(0,covs), ylim = c(0,100))

  filename=paste(wgfile,"cumulative_coverage.png",sep=".")
  #filename=sub(".coverage.counts.txt", "", filename, ignore.case =FALSE, fixed=FALSE)
  ggsave(filename, width = 20, height = 10, units = "cm",dpi=150)

  names(we.gd.prop)=c("Whole Genome","Exome","Whole Genome Cum.","Exon Cum.","cov")
  m.we.gd.prop = melt(we.gd.prop[,c(1,2,5)],id="cov")
  print(head(m.we.gd.prop))
  ymax=max(subset(m.we.gd.prop,cov!=0)$value)
  ggplot(m.we.gd.prop,aes(cov,value,fill=variable))+
    geom_bar(stat="identity",position="dodge")+
    ggtitle("All Chromosomes Coverage Distribution")+
    scale_y_continuous("Percent Total")+
    scale_fill_manual("Region",values=gel_colours)+
    scale_x_continuous("Coverage")+
    theme_gel_proper()+
    coord_cartesian(xlim = c(0,covs), ylim = c(0,ymax))

  filename=paste(wgfile,"all.coverage.distribution.all.chrs.png",sep=".")
  #filename=sub(".coverage.counts.txt", "", filename, ignore.case =FALSE, fixed=FALSE)
  ggsave(filename, width = 20, height = 10, units = "cm",dpi=150)


}


#' A function to plot a gene with coverage and variants
#'
#'
#' @param gene the gene you want to print
#' @param bw the bw file
#' @param vcf the vcf file
#' @keywords plot
#' @export
#' @examples
#' plot_gene_coverage()
plot_gene_coverage <- function(gene,bw,vcf){
    data(genesymbol, package = "biovizBase")
    wh <- genesymbol[c(gene)]
    wh <- range(wh, ignore.strand = TRUE)
    p.tx <- autoplot(Homo.sapiens, which  = wh,columns = c("TXNAME"), names.expr = "TXNAME", label=FALSE)
    p.tx = p.tx+
      theme_gel_proper()

    #vcf plot
    vcf <- readVcf(vcf, "hg19")
    vr <- as(vcf, "VRanges")
    p.vr <- autoplot(vr, which = wh,geom="rect",facet=FALSE,arrow=FALSE)
    p.vr = p.vr+
      ggtitle("test")+
      theme_gel_proper()+
      theme(legend.position="None",
            strip.text=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())+
      scale_y_continuous(limits = c(1.6,1.8))+
      scale_fill_manual(values=c("firebrick2","firebrick2","firebrick2","firebrick2","firebrick2"))+
      scale_color_manual(values=c("firebrick2","firebrick2","firebrick2","firebrick2","firebrick2"))

    #bw plot
    bw.data=summary(BigWigFile(bw),size=1000,which=wh)
    p.bw <- autoplot(bw.data,geom="bar",fill="red",color="red")
    p.bw <- p.bw+
      scale_y_continuous("Coverage")+
      theme_gel_proper()+
      geom_hline(yintercept=30)

    tracks(p.vr,p.bw,p.tx,heights = c(1, 3, 4),title=gene)
}