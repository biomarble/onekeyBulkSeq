args<-commandArgs(trailingOnly = T)
options(scipen = 200)
dotfile=args[1]
winwidth=as.numeric(args[2])
winshift=as.numeric(args[3])
mincount=as.numeric(args[4])
outfile=args[5]
useSim=as.numeric(args[6])

data=read.table(dotfile,header=F,stringsAsFactors = F)
myncol <- ncol(data)
if(useSim==1){
  colnames(data)=c('chr','pos','index1','index2','delta','L95','H95','L99','H99')
}else{
  colnames(data)=c('chr','pos','index1','index2','delta')
}
chrname=unique(data$chr)
chrnum=length(chrname)

if(file.exists(outfile)){
  cat('removing existing file', outfile,"\n")
  file.remove(outfile)
}
if(useSim==1){
head='#Chr\tPos\tSNPindexMut\tSNPindexWT\t△SNPindex\twindSize\tCount\tStart\tEnd\tL95\tH95\tL99\tH99'
}else{
head='#Chr\tPos\tSNPindexMut\tSNPindexWT\t△SNPindex\twindSize\tCount\tStart\tEnd'
}
write(head,file=outfile)

for (nn in 1:chrnum){
  cat('sliding: chromosome ',chrname[nn],'\r')
  chrdata <- data[data$chr == chrname[nn],]
  if(nrow(chrdata) < 1) next
  chrsize=max(chrdata$pos)
  steps <- ceiling(chrsize / winshift)
  for (ss in 1:steps){
    pos0 <- winshift * (ss - 1)
    pos1 <- pos0 - winwidth / 2
    pos2 <- pos0 + winwidth / 2
    if(pos1 < 0) pos1 <- 0
    if(pos2 > chrsize) pos2 <- chrsize
    actwidth <- pos2 - pos1 + 1
    chrposdata <- chrdata[(chrdata$pos>=pos1 & chrdata$pos<=pos2), ]
    chrposdata_mean1 <- apply(chrposdata['index1'], 2, mean)
    chrposdata_mean2 <- apply(chrposdata['index2'], 2, mean)
    chrposdata_mean <- apply(chrposdata['delta'], 2, mean)
    if(useSim==1){
   	L95 <- apply(chrposdata['L95'], 2, mean)
   	H95 <- apply(chrposdata['H95'], 2, mean)
   	L99 <- apply(chrposdata['L99'], 2, mean)
   	H99 <- apply(chrposdata['H99'], 2, mean)
        chrposdata_cnt  <- nrow(chrposdata)
        newdata <- c(chrname[nn], pos0,chrposdata_mean1,chrposdata_mean2, chrposdata_mean,actwidth, chrposdata_cnt,pos1,pos2,L95,H95,L99,H99)
    }else{
        chrposdata_cnt  <- nrow(chrposdata)
        newdata <- c(chrname[nn], pos0,chrposdata_mean1,chrposdata_mean2, chrposdata_mean,actwidth, chrposdata_cnt,pos1,pos2)
    }
    
    if(mincount > 0){
      if(chrposdata_cnt < mincount){
        chrposdata_null <- NaN
        if(useSim==1){
            newdata <- c(chrname[nn], pos0, chrposdata_null,chrposdata_null, chrposdata_null,actwidth, chrposdata_cnt,pos1,pos2,L95,H95,L99,H99)
	}else{
            newdata <- c(chrname[nn], pos0, chrposdata_null,chrposdata_null, chrposdata_null,actwidth, chrposdata_cnt,pos1,pos2)
	}
      }else{
       write(t(as.matrix(newdata)), file = outfile, ncolumns = myncol + 10, append = TRUE, sep = "\t")
      }
    }
  }
}
