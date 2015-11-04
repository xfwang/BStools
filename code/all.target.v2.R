# This R-script is used to generate the mean and median histograms at the target-region level
# as well as generate a target summary table

args<-commandArgs(trailingOnly = TRUE)

 chrData <- read.table(file=paste(args[1],"/",args[3],".",args[2],".nonCG.covered.out",sep=''), header=F)
  colnames(chrData)<-c("chr","start","end","tot.nonCGc","methylation")

  m.chrData <- read.table(file=paste(args[1],"/",args[3],".",args[2],".nonCG.match",sep=''), header=F)
  colnames(m.chrData)<-c("chr","start","end","tot.nonCGc","methylation")

  splc<-split(chrData, paste(chrData$chr, chrData$start, chrData$end))
  m.splc<-split(m.chrData, paste(m.chrData$chr, m.chrData$start, m.chrData$end))

    
  df.summ<-as.data.frame(t(sapply(splc, function(x) summary(1-x$methylation))))
  m.summ.chrData<-as.data.frame(t(sapply(m.splc, function(z) list(chr=z$chr[1], start=z$start[1], end=z$end[1], NCGC=NROW(z$end), NCGCwC=sum(z[['tot.nonCGc']]>=1), percent=sum(z[['tot.nonCGc']]>=1)/NROW(z$end) )     )))
  summ<-cbind(m.summ.chrData, df.summ)



out<-as.data.frame(lapply(summ,function(x) factor(unlist(x))))
write.table(out, file=paste(args[1],"/",args[3],".",args[2],".target.summary.table.txt",sep=''), quote=F, row.names=F, sep="\t", col.names=T)


postscript(file=paste(args[1],"/",args[3],".",args[2],".mean.median.ps",sep=''), paper="letter", horizontal=T)
par(mfrow=c(1,2))

hist(summ$Median, prob=T, xlab="Bisulfite rate", main=paste(args[2]," Median Bisulfite Rate",sep=''), breaks=seq(0,1, by=0.05), density=10, angle=45, cex.main=2, cex.lab=1.5, cex.axis=1.5)
hist(summ$Mean, prob=T, xlab="Bisulfite rate", main=paste(args[2]," Mean Bisulfite Rate",sep=''), breaks=seq(0,1, by=0.05), density=10, angle=45,cex.main=2, cex.lab=1.5, cex.axis=1.5)
dev.off()




