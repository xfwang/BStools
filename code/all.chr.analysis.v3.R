# This R-script is used to generate the bisulfite rate histogram at the chromosome-level
# as well as a small summary table for whole genome and the individual chromosomes

args<-commandArgs(trailingOnly = TRUE)

chrData <- read.table(file=paste(args[1],"/",args[2],".nonCG.covered.base",sep=''), header=F)
colnames(chrData)<-c("chr","methylation")

chrTR <- length(count.fields(file=paste(args[1],"/",args[2],".nonCG.base",sep='')))

summ<-c(as.numeric(chrTR),as.integer(dim(chrData)[1]), as.numeric((dim(chrData)[1])/chrTR), as.numeric(summary(1-chrData[,2])))
summ<-as.data.frame(t(summ))
splc<-split(chrData, paste(chrData$chr))
#df.chrData<-as.data.frame(t(sapply(splc, function(x) list(chr=x$chr[1], TNCGC=chrTR, TNCGCwC=NROW(x$methylation), percent=(NROW(x$methylation)/chrTR) ))))
#summ.df.chrData<-as.data.frame(t(sapply(splc, function(x) summary(1-x$methylation))))
#summ<-cbind(df.chrData, summ.df.chrData)

colnames(summ)<-c("nonCGC", "nonCGC.wC", "percent", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")

write.table(summ, file=paste(args[1],"/","summary.table.txt",sep=''), quote=F, row.names=F, sep="\t", col.names=T)


postscript(file=paste(args[1],"/","BS.ps",sep=''), paper="letter", horizontal=T)
par(mfrow=c(1,1))
     hist(1-chrData[,2],prob=T,xlab="Bisulfite rate", main=paste("Whole Genome Bisulfite Rate",sep=''), breaks=seq(0,1, by=0.05), density=10, angle=45, cex.main=2, cex.lab=1.5, cex.axis=1.5)
#dev.off()


#postscript(file=paste(args[1],"/","BS.chr.ps",sep=''), paper="letter", horizontal=T)
par(mfrow=c(4,6))
for (i in 1:length(splc)) {
     chrName<-splc[[i]][1,1]
     hist(1-splc[[i]][,2],prob=T,xlab="Bisulfite rate", main=paste(splc[[i]][1,1]," Bisulfite Rate",sep=''), breaks=seq(0,1, by=0.05), density=10, angle=45, cex.main=1, cex.lab=1, cex.axis=1)
}
dev.off()


chrCG <- read.table(file=paste(args[1],"/",args[2],".CG.base",sep=''), header=F)
colnames(chrCG)<-c("chr", "pos", "cov", "methylation")
chrCG.cover <- chrCG[chrCG[,3]>0,]


summ2<-c(as.numeric(dim(chrCG)[1]),as.numeric(dim(chrCG.cover)[1]), signif(as.numeric(dim(chrCG.cover)[1])/as.numeric(dim(chrCG)[1]),4), as.numeric(summary(as.numeric(chrCG.cover[,3]))))
summ2<-as.data.frame(t(summ2))

colnames(summ2)<-c("CGC", "CGC.wC", "percent", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")

write.table(summ2, file=paste(args[1],"/","CG.coverage.table.txt",sep=''), quote=F, row.names=F, sep="\t", col.names=T)



postscript(file=paste(args[1],"/","mC.distribution.ps",sep=''), paper="letter", horizontal=T)
par(mfrow=c(1,1))
hist(chrCG.cover[,4], prob=T,xlab="mC level", main="whole genome", breaks=seq(0,1, by=0.05), density=10, angle=45, cex.main=2, cex.lab=1.5, cex.axis=1.5)
splc.CG<-split(chrCG.cover, paste(chrCG.cover$chr))
par(mfrow=c(4,6))
for (i in 1:length(splc.CG)) {
     chrName<-splc.CG[[i]][1,1]
     hist(1-splc.CG[[i]][,2],prob=T,xlab="mC level", main=paste(splc.CG[[i]][1,1]," mC level",sep=''), breaks=seq(0,1, by=0.05), density=10, angle=45, cex.main=1, cex.lab=1, cex.axis=1)
}
dev.off()



