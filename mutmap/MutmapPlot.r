args<-commandArgs(trailingOnly = T)
suppressPackageStartupMessages(library(ggplot2))

linefile=args[1]
dotfile=args[2]
cutoff=args[3]
title=args[4]
outprefix=args[5]

if(cutoff!="NULL"){
	cutoff=as.numeric(cutoff)
}else{
	cutoff=NULL
}

data=read.table(linefile,header=F,stringsAsFactors = T)
colnames(data)=c('chr','pos','index')

minchrLength=2000000
CheckChr=tapply(data$pos, data$chr, max)>minchrLength
data=data[data$chr %in% names(CheckChr)[CheckChr],]

dots=NULL
if (dotfile!="NULL" ){ 
  dots=read.table(dotfile,header=F,stringsAsFactors = T)
  colnames(dots)=c('chr','pos','index')
  dots=dots[dots$chr %in% names(CheckChr)[CheckChr],]
}
nchr=length(unique(data$chr))
colors=rainbow(n=nchr)
######manhattan
g=ggplot(data,aes(x=pos/1000000,y=index))
if (! is.null(dots) ){ 
  g=g+geom_point(data=dots,aes(color=as.factor(dots$chr)),size=0.8,alpha=0.6)
}
g=g+geom_line(color='black',size=0.8)+
  facet_grid(.~chr,scales="free",switch='both')+
  scale_fill_manual(values=colors)+
  scale_x_continuous(breaks=pretty(data$pos/1000000,n=10))+
  theme_bw()+
  labs(x="chromosome position(M)",y=element_blank(),title=title)+
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
if(! is.null(cutoff)){
  g=g+geom_hline(aes(yintercept = cutoff) ,color='red',linetype='dashed')
}

#png(paste0(outprefix,'.total.png'),w=6000,h=2000,res=300,units="px")
pdf(paste0(outprefix,'.total.pdf'),w=nchr*5,h=6)
print(g)
dev.off()

#####################single chromosome manhattan
singlemanhattan <-function(linedata,dotdata,a,prefix,cutoff){
  line.forplot=linedata[linedata$chr==a,]
  
  if (! is.null(dotdata) ){ 
    dots.forplot=dotdata[dotdata$chr==a,]
  }
  g=ggplot(line.forplot,aes(x=pos/1000000,y=index))
  
  if (! is.null(dots) ){ 
    g=g+geom_point(data=dots.forplot,color='blue',alpha=0.6)
  }
  
  g=g+geom_line(size=0.8)+
    scale_x_continuous(breaks = pretty(line.forplot$pos/1000000,n=20))+
    theme_bw()+
    labs(x="chromosome position(M)",y=element_blank(),title=a)+
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
  if(! is.null(cutoff)){
    g=g+geom_hline(aes(yintercept = cutoff) ,color='red',linetype='dashed')
  }
  
  #png(paste0(prefix,".",a,".png"),w=4000,h=2000,res=300,units="px")
  pdf(paste0(prefix,".",a,".pdf"),w=8,h=6)
  print(g)
  dev.off()
}

for(i in unique(data$chr)){
  cat('plotting ',i,'\r')
  singlemanhattan(data,dots,i,outprefix,cutoff)
}
