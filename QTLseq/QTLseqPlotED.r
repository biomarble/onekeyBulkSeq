args<-commandArgs(trailingOnly = T)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))

linefile=args[1]
dotfile=args[2]
cutoff=args[3]
title=args[4]
outprefix=args[5]

# 
# rm(list=ls())
# linefile="../../BSAoutput/2.QTLseq.windowed.ED.txt"
# dotfile="../../BSAoutput/1.QTLseq.raw.ED.txt"
# cutoff="NULL"
# title=""
# outprefix="test"

if(cutoff!="NULL"){
	cutoff=as.numeric(cutoff)
}else{
	cutoff=NULL
}

data=read.table(linefile,header=F,stringsAsFactors = T)
Usenames=c('chr','pos','ED')
colnames(data)=Usenames

minchrLength=2000000
CheckChr=tapply(data$pos, data$chr, max)>minchrLength
data=data[data$chr %in% names(CheckChr)[CheckChr],Usenames]

dots=read.table(dotfile,header=F,stringsAsFactors = T,col.names=c('chr','pos','oriED','ED'))
dots=dots[dots$chr %in% names(CheckChr)[CheckChr],Usenames]
nchr=length(unique(data$chr))

colors=rep(c('#00A087FF','#008B45FF'),ceiling(nchr/2))[1:nchr]



if(is.null(cutoff)){
    cutoff=median(data$ED)+3*sd(data$ED)
}



######manhattan
g=ggplot(data,aes(x=pos/1000000,y=ED))+
    geom_point(data=dots,aes(color=chr),size=0.7,alpha=0.6)+
    geom_line(color='black',size=0.8)+
    geom_hline(aes(yintercept = cutoff) ,color='red',linetype='dashed',size=0.8)

g=g+facet_grid(~chr,scales="free",switch='both')+
  scale_color_manual(values=colors)+
  scale_x_continuous(breaks=pretty(data$pos/1000000,n=10))+
  theme_classic()+
  labs(x=paste0("chromosome position(Mb)  cutoff:",sprintf('%.4f',cutoff)),y="Euclidean Distance",title=title)+
  theme(
    legend.position = 'none',
    panel.background = element_blank(),  
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    strip.placement = "outside",
    panel.grid.major.y  = element_line(colour = "grey80",size=.25,linetype ="dotted" ),
    panel.grid.major.x  = element_line(colour = "grey80",size=.25,linetype ="dotted" ),
    panel.grid.minor.y  = element_blank(),
    panel.grid.minor.x  = element_blank(),
    plot.title=element_text(size=15,hjust=.5),
    axis.text.x = element_text(size = 8,colour="black"),
    axis.text.y = element_text(size = 10,colour="black")
  )

pdf(paste0(outprefix,'.total.pdf'),w=nchr*5,h=6)
print(g)
dev.off()

#####################single chromosome manhattan
singlemanhattan <-function(linedata,dotdata,a,prefix,cutoff){
  line.forplot=linedata[linedata$chr==a,]
  dots.forplot=dotdata[dotdata$chr==a,]
  
  g=ggplot(line.forplot,aes(x=pos/1000000,y=ED))+
      geom_point(data=dots.forplot,color='blue',size=0.7,alpha=0.8)+
      geom_line(color='black',size=0.8)+
      geom_hline(aes(yintercept = cutoff) ,color='red',linetype='dashed')
  
  g=g+facet_grid(~chr,scales="free",switch='both')+
    scale_x_continuous(breaks = pretty(line.forplot$pos/1000000,n=20))+
    theme_classic()+
    labs(x=paste0("chromosome position(Mb) cutoff:",sprintf('%.4f',cutoff)),y="Euclidean Distance",title=paste(title,a))+
    theme(
      legend.position = 'none',
      panel.background = element_blank(),  
      strip.background = element_blank(),
      strip.text = element_text(size = 12),
      strip.placement = "outside",
      panel.grid.major.y  = element_line(colour = "grey80",size=.25,linetype ="dotted" ),
      panel.grid.major.x  = element_line(colour = "grey80",size=.25,linetype ="dotted" ),
      panel.grid.minor.y  = element_blank(),
      panel.grid.minor.x  = element_blank(),
      plot.title=element_text(size=15,hjust=.5),
      axis.text.x = element_text(size = 8,colour="black"),
      axis.text.y = element_text(size = 10,colour="black")
    )

  pdf(paste0(prefix,".",a,".pdf"),w=8,h=6)
  print(g)
  dev.off()
}

for(i in unique(data$chr)){
  cat('plotting ',i,'\r')
  singlemanhattan(data,dots,i,outprefix,cutoff)
}
