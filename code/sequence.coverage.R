# The R-script is used to generate sequencing boxplots comparing regions with 
# high and low values of coverage

# Note, in the following script, we use 0.05 and 20 as the cut-off values. We do this 
# for the case that one group (e.g., high BS rate regions/group) has dramatically more or less 
# target regions (or lines, e.g., 28,000) than another group (e.g., low BS rate regions/group, e.g., 
# just 60). If this is the case, we will sample the smaller number of (e.g., 60) from the larger one
# (e.g., 28,000) to make the boxplots. 

args<-commandArgs(trailingOnly = TRUE)

if (args[4] < 0.05) {
highPer<-read.table(file=paste(args[1],"/",args[3],".",args[2],".highCoverage.seq",sep=''), header=T)
lPer<-read.table(file=paste(args[1],"/",args[3],".",args[2],".lowCoverage.seq",sep=''), header=T)
hPer<-highPer[ sample(NROW(highPer),args[6]), ]
}

if (args[4] > 20) {
hPer<-read.table(file=paste(args[1],"/",args[3],".",args[2],".highCoverage.seq",sep=''), header=T)
lowPer<-read.table(file=paste(args[1],"/",args[3],".",args[2],".lowCoverage.seq",sep=''), header=T)
lPer<-lowPer[ sample(NROW(lowPer),args[5]), ]
}

if (args[4] >= 0.05 & args[4] <= 20) {
hPer<-read.table(file=paste(args[1],"/",args[3],".",args[2],".highCoverage.seq",sep=''), header=T)
lPer<-read.table(file=paste(args[1],"/",args[3],".",args[2],".lowCoverage.seq",sep=''), header=T)
}


# Below I am creating an index for each percentage column 
# so that I can combine the high/low results on one plot
a.Per.index<-c(rep("high", length(hPer$pa)), rep("low", length(lPer$pa)))
c.Per.index<-c(rep("high", length(hPer$pc)), rep("low", length(lPer$pc)))
g.Per.index<-c(rep("high", length(hPer$pg)), rep("low", length(lPer$pg)))
t.Per.index<-c(rep("high", length(hPer$pt)), rep("low", length(lPer$pt)))
cpg.Per.index<-c(rep("high", length(hPer$pc.g)), rep("low", length(lPer$pc.g)))
cg.Per.index<-c(rep("high", length(hPer$pcgc)), rep("low", length(lPer$pcgc)))
ncg.Per.index<-c(rep("high", length(hPer$pncgc)), rep("low", length(lPer$pncgc)))
lc.Per.index<-c(rep("high", length(hPer$X.lower_count)), rep("low", length(lPer$X.lower_count)))


a.Per.all<-c(hPer$pa, lPer$pa)
c.Per.all<-c(hPer$pc, lPer$pc)
g.Per.all<-c(hPer$pg, lPer$pg)
t.Per.all<-c(hPer$pt, lPer$pt)
cpg.Per.all<-c(hPer$pc.g, lPer$pc.g)
cg.Per.all<-c(hPer$pcgc, lPer$pcgc)
ncg.Per.all<-c(hPer$pncgc, lPer$pncgc)
lc.Per.all<-c(hPer$X.lower_count, lPer$X.lower_count)


# This is the creation of the boxplot

postscript(file=paste(args[1],"/",args[3],".",args[2],".seq.coverage.boxplot.ps",sep=''), paper="letter", horizontal=T)
par(mfrow=c(2,4))

# Percent
boxplot(a.Per.all~a.Per.index, main=" %A", cex.main=2, cex.lab=1.5, cex.axis=1.5)
boxplot(c.Per.all~c.Per.index, main=" %C", cex.main=2, cex.lab=1.5, cex.axis=1.5)
boxplot(g.Per.all~g.Per.index, main=" %G", cex.main=2, cex.lab=1.5, cex.axis=1.5)
boxplot(t.Per.all~t.Per.index, main=" %T", cex.main=2, cex.lab=1.5, cex.axis=1.5)
boxplot(cpg.Per.all~cpg.Per.index, main=" %C+G", cex.main=2, cex.lab=1.5, cex.axis=1.5)
boxplot(cg.Per.all~cg.Per.index, main=" %CGc", cex.main=2, cex.lab=1.5, cex.axis=1.5)
boxplot(ncg.Per.all~ncg.Per.index, main=" %nonCGc", cex.main=2, cex.lab=1.5, cex.axis=1.5)
boxplot(lc.Per.all~lc.Per.index, main=" %low_count", cex.main=2, cex.lab=1.5, cex.axis=1.5)
dev.off()




