args<-commandArgs(trailingOnly = T)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))

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

#rm(list=ls())
#linefile="../../BSAoutput/2.QTLseq.windowed.index.txt"
#dotfile="../../BSAoutput/1.QTLseq.raw.index.txt"
#title=""
#outprefix="../../BSAoutput/test"
#cutoff="NULL"


Usenames=c('chr','pos','SNPindexMut','SNPindexWT',"deltaSNPindex")
UseLables=c(expression('SNP-index '[Mut]),expression('SNP-index '[WT]),expression(paste(Delta,'SNP-index')))

dots=read.table(dotfile,header=F,stringsAsFactors = T)
colnames(dots)=Usenames
minchrLength=2000000
CheckChr=tapply(dots$pos, dots$chr, max)>minchrLength


dots=dots[dots$chr %in% names(CheckChr)[CheckChr],Usenames]
usedot=melt(dots,id.vars=c('chr','pos'))
usedot$variable <- factor(usedot$variable,
                        levels=c('SNPindexMut','SNPindexWT','deltaSNPindex'),
                        labels=UseLables
)

data=read.table(linefile,header=F,stringsAsFactors = T)

if(!is.null(cutoff)){
    cutoff=as.numeric(cutoff)
    data=data[data[,1] %in% names(CheckChr)[CheckChr],1:5]
    colnames(data)=Usenames
    cutData=data[,c('chr','pos')]
    cutData$cutoff=cutoff
    cutData$variable=factor('deltaSNPindex',
                            levels=c('SNPindexMut','SNPindexWT','deltaSNPindex'),
                            labels=UseLables
    )
    
}else{
    data=data[data[,1] %in% names(CheckChr)[CheckChr],c(1:5,10,11,12,13)]
    colnames(data)=c(Usenames,'L95','H95','L99','H99')
    simData=data[,c('chr','pos','L95','H95','L99','H99')]
    simData$variable=factor('deltaSNPindex',
                            levels=c('SNPindexMut','SNPindexWT','deltaSNPindex'),
                            labels=UseLables
                            )
    data=data[,Usenames]
}

nchr=length(unique(data$chr))



colors=rep(c('#00A087FF','#008B45FF'),ceiling(nchr/2))[1:nchr]


useline=melt(data,id.vars=c('chr','pos'))

useline$variable <- factor(useline$variable,
                       levels=c('SNPindexMut','SNPindexWT','deltaSNPindex'),
                       labels=UseLables
)


limDat=data.frame(chr=rep(names(CheckChr)[1],6),
                  pos=c(1,1,1,2,2,2),
                  variable=factor(rep(Usenames[3:5],2),
                                  levels=c('SNPindexMut','SNPindexWT','deltaSNPindex'),
                                  labels=UseLables
                                  ),
                  value=c(0,0,-1,1,1,1)
                  )

######manhattan
g=ggplot(useline,aes(x=pos/1000000,y=value))+
    geom_blank(data=limDat,aes(x=pos/1000000,y=value))+
    geom_point(data=usedot,aes(color=chr),size=0.7,alpha=0.6)+
    geom_line(color='black',size=0.8)

if(is.null(cutoff)){
    g=g+geom_line(data=simData,aes(x=pos/1000000,y=L95),color='#5F559BFF',size=0.6)+
        geom_line(data=simData,aes(x=pos/1000000,y=H95),color='#5F559BFF',size=0.6)+
        geom_line(data=simData,aes(x=pos/1000000,y=L99),color='#EE0000FF',size=0.6)+
        geom_line(data=simData,aes(x=pos/1000000,y=H99),color='#EE0000FF',size=0.6)
}else{
    g=g+geom_hline(data=cutData,aes(yintercept = cutoff), color='red',linetype='dashed',size=0.8)
}

g=g+facet_grid(variable~chr,scales="free",switch='both',labeller = label_parsed)+
  scale_color_manual(values=colors)+
  scale_x_continuous(breaks=pretty(data$pos/1000000,n=10))+
  theme_classic()+
  labs(x="chromosome position(Mb)",y=element_blank(),title=title)+
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
  limDat=data.frame(chr=rep(a,6),
                    pos=c(1,1,1,2,2,2),
                    variable=factor(rep(Usenames[3:5],2),
                                    levels=c('SNPindexMut','SNPindexWT','deltaSNPindex'),
                                    labels=UseLables
                    ),
                    value=c(0,0,-1,1,1,1)
  )
  g=ggplot(line.forplot,aes(x=pos/1000000,y=value))+
      geom_blank(data=limDat,aes(x=pos/1000000,y=value))+
      geom_point(data=dots.forplot,color='blue',size=0.7,alpha=0.8)+
      geom_line(color='black',size=0.8)
  
  if(is.null(cutoff)){
      sim=simData[simData$chr==a,]
      g=g+geom_line(data=sim,aes(x=pos/1000000,y=L95),color='#5F559BFF',size=0.6)+
          geom_line(data=sim,aes(x=pos/1000000,y=H95),color='#5F559BFF',size=0.6)+
          geom_line(data=sim,aes(x=pos/1000000,y=L99),color='#EE0000FF',size=0.6)+
          geom_line(data=sim,aes(x=pos/1000000,y=H99),color='#EE0000FF',size=0.6)
  }else{
      cut=cutData[cutData$chr==a,]
      g=g+geom_hline(data=cut,aes(yintercept = cutoff) ,color='red',linetype='dashed',size=0.8)
  }
    g=g+facet_grid(variable~chr,scales="free",switch='both',labeller = label_parsed)+
    scale_x_continuous(breaks = pretty(line.forplot$pos/1000000,n=20))+
    theme_classic()+
    labs(x="chromosome position(Mb)",y=element_blank(),title=a)+
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
  singlemanhattan(useline,usedot,i,outprefix,cutoff)
}
