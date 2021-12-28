args<-commandArgs(trailingOnly = T)
options(scipen = 200)
data=read.table(args[1],header=F,stringsAsFactors = F)
colnames(data)=c('chr','pos','index')
chrname=unique(data$chr)
chrnum=length(chrname)
winwidth=as.numeric(args[2])
winshift=as.numeric(args[3])
mincount=as.numeric(args[4])
outfile=args[5]
myncol <- ncol(data)

if(file.exists(outfile)){
  cat('removing existing file', outfile,"\n")
  file.remove(outfile)
}
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
    chrposdata_mean <- apply(chrposdata[3], 2, mean)
    chrposdata_cnt  <- nrow(chrposdata)
    newdata <- c(chrname[nn], pos0, chrposdata_mean,actwidth,  chrposdata_cnt,pos1,pos2)
#    write(newdata, file = outfile, append = T, ncolumns = myncol + 10,sep = "\t")
    
    if(mincount > 0){
      if(chrposdata_cnt < mincount){
        chrposdata_null <- NaN
        newdata <- c(chrname[nn], pos0, actwidth, chrposdata_null, chrposdata_cnt,pos1,pos2)
      }else{
       write(t(as.matrix(newdata)), file = outfile, ncolumns = myncol + 10, append = TRUE, sep = "\t")
      }
    }
  }
}
