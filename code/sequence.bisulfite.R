# This R-script is used to generate the bisulfite-rate sequencing boxplot
# comparing regions with high and low bisulfite conversion rates

# Note, in the following script, we use 0.05 and 20 as the cut-off values. We do this 
# for the case that one group (e.g., high BS rate regions/group) has dramatically more or less 
# target regions (or lines, e.g., 28,000) than another group (e.g., low BS rate regions/group, e.g., 
# just 60). If this is the case, we will sample the smaller number of (e.g., 60) from the larger one
# (e.g., 28,000) to make the boxplots. 

args<-commandArgs(trailingOnly = TRUE)

if (args[4] < 0.05) {
highBS<-read.table(file=paste(args[1],"/",args[3],".",args[2],".highBS.seq",sep=''), header=T)
lBS<-read.table(file=paste(args[1],"/",args[3],".",args[2],".lowBS.seq",sep=''), header=T)
hBS<-highBS[ sample(NROW(highBS),args[6]), ]
}

if (args[4] > 20) {
hBS<-read.table(file=paste(args[1],"/",args[3],".",args[2],".highBS.seq",sep=''), header=T)
lowBS<-read.table(file=paste(args[1],"/",args[3],".",args[2],".lowBS.seq",sep=''), header=T)
lBS<-lowBS[ sample(NROW(lowBS),args[5]), ]
}

if (args[4] >= 0.05 & args[4] <= 20) {
hBS<-read.table(file=paste(args[1],"/",args[3],".",args[2],".highBS.seq",sep=''), header=T)
lBS<-read.table(file=paste(args[1],"/",args[3],".",args[2],".lowBS.seq",sep=''), header=T)
}


# Below I am creating an index for each percentage column 
# so that I can combine the high/low results on one plot
a.BS.index<-c(rep("high", length(hBS$pa)), rep("low", length(lBS$pa)))
c.BS.index<-c(rep("high", length(hBS$pc)), rep("low", length(lBS$pc)))
g.BS.index<-c(rep("high", length(hBS$pg)), rep("low", length(lBS$pg)))
t.BS.index<-c(rep("high", length(hBS$pt)), rep("low", length(lBS$pt)))
cpg.BS.index<-c(rep("high", length(hBS$pc.g)), rep("low", length(lBS$pc.g)))
cg.BS.index<-c(rep("high", length(hBS$pcgc)), rep("low", length(lBS$pcgc)))
ncg.BS.index<-c(rep("high", length(hBS$pncgc)), rep("low", length(lBS$pncgc)))
lc.BS.index<-c(rep("high", length(hBS$X.lower_count)), rep("low", length(lBS$X.lower_count)))


a.BS.all<-c(hBS$pa, lBS$pa)
c.BS.all<-c(hBS$pc, lBS$pc)
g.BS.all<-c(hBS$pg, lBS$pg)
t.BS.all<-c(hBS$pt, lBS$pt)
cpg.BS.all<-c(hBS$pc.g, lBS$pc.g)
cg.BS.all<-c(hBS$pcgc, lBS$pcgc)
ncg.BS.all<-c(hBS$pncgc, lBS$pncgc)
lc.BS.all<-c(hBS$X.lower_count, lBS$X.lower_count)


# This is the creation of the boxplot
postscript(file=paste(args[1],"/",args[3],".",args[2],".seq.bisulfite.boxplot.ps",sep=''), paper="letter", horizontal=T)
par(mfrow=c(2,4))

# BS
boxplot(a.BS.all~a.BS.index, main=" %A", cex.main=2, cex.lab=1.5, cex.axis=1.5)
boxplot(c.BS.all~c.BS.index, main=" %C", cex.main=2, cex.lab=1.5, cex.axis=1.5)
boxplot(g.BS.all~g.BS.index, main=" %G", cex.main=2, cex.lab=1.5, cex.axis=1.5)
boxplot(t.BS.all~t.BS.index, main=" %T", cex.main=2, cex.lab=1.5, cex.axis=1.5)
boxplot(cpg.BS.all~cpg.BS.index, main=" %C+G", cex.main=2, cex.lab=1.5, cex.axis=1.5)
boxplot(cg.BS.all~cg.BS.index, main=" %CGc", cex.main=2, cex.lab=1.5, cex.axis=1.5)
boxplot(ncg.BS.all~ncg.BS.index, main=" %nonCGc", cex.main=2, cex.lab=1.5, cex.axis=1.5)
boxplot(lc.BS.all~lc.BS.index, main=" %low_count", cex.main=2, cex.lab=1.5, cex.axis=1.5)
dev.off()