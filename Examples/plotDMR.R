# function and example code for DMR visualization
# 
# Main function plotRegion()
# 1) generate summarized information for all tumor and normal sampels in targeted region, as a lit. The first element
#    of the list the 
# 2) plot methylation levels for all tumor and normal sampels in targeted region

##############################################################
# example script
##############################################################

install.packages("plotrix")
library("plotrix")

# sample names of two groups
tumor_samples<-c("T1","T2","T3","T4","T5","T6") 
normal_samples<-c("N1","N2","N3","N4","N5","N6")


example.region<-plotRegion(normal, tumor, "chr5", 76249822, 76249807, "DMR.example.ps")

# tumor samples summary
get.summary.table(region_tumor_coverage,region_tumor_mC, tumor_samples)
                 T1    T2   T3    T4    T5     T6
covered_CG    24.00 25.00 23.0 25.00 24.00  26.00
covered_CG(%) 92.00 96.00 88.0 96.00 92.00 100.00
mean_mC        0.76  0.77  0.7  0.66  0.48   0.74
mean_coverage  8.00  7.00  6.8  3.80  3.70   4.20

# normal samples summary
get.summary.table(region_normal_coverage,region_normal_mC, normal_samples)
                 N1    N2     N3    N4     N5     N6
covered_CG    25.00 25.00 25.000 25.00 25.000 25.000
covered_CG(%) 96.00 96.00 96.000 96.00 96.000 96.000
mean_mC        0.13  0.16  0.058  0.07  0.064  0.094
mean_coverage  5.40  5.20  3.600  7.60  6.200  7.200

example.region<-plot.gene(region_tumor_mC, region_tumor_coverage, region_normal_mC, region_normal_coverage,  "example.DMR.ps", tumor_samples, normal_samples , 0.32, main="exmaple DMR")



##############################################################
#  visualization functions
##############################################################

# get.summary.table
# this function summarizes number and percentage of covered CG, mean methylation levekm and mean coverage for provided region

get.summary.table<-function(coverage.mat, mC.mat, sample.name){
  sumary.table<-matrix(NA, ncol=4, nrow=0)
  for (i in 3:dim(coverage.mat)[2]){
     nocov.index<-(1:dim(coverage.mat)[1])[coverage.mat[,i]==0]
     covered_CG<-dim(coverage.mat)[1]-length( nocov.index)
     covered_CG_percent<-((dim(coverage.mat)[1]-length(nocov.index))/dim(coverage.mat)[1])*100
     if ((dim(coverage.mat)[1]-length(nocov.index)) >0){
       if (length(nocov.index)>0) {
           mean_mC<-mean(mC.mat[-nocov.index,i])
           mean_cov<-mean(coverage.mat[-nocov.index,i])    
          }
      else {
           mean_mC<-mean(mC.mat[,i])
           mean_cov<-mean(coverage.mat[,i])    

          }
     }
     else {
        mean_mC<-NA
        mean_cov<-0
     }
     vec<-c(covered_CG, covered_CG_percent,mean_mC, mean_cov)
     sumary.table<-rbind(sumary.table, vec)
  }
  final<-t(sumary.table)
  rownames(final) <- c("covered_CG", "covered_CG(%)", "mean_mC", "mean_coverage")
  colnames(final)<-sample.name
  return(signif(final,2))
}


# plot.gene
# this function plots the methylation level of two groups (tumor and normal samples) for a provied region, showing the difference in methylation patterns of the two groups
plot.gene<-function( region_tumor_mc, region_tumor_coverage, region_normal_mc, region_normal_coverage, ps.name, tumor_samples, normal_samples, pie.radius=0.5, main){
       postscript(ps.name, paper="letter")
       par(mar=c(7,4,5,2))
       plot_CG<-region_tumor_mc[,c(1,2)]
       for (i in 1:length(normal_samples)) {
            CG<-region_normal_mc[,i+2]
            CG_0X<-(1:length(CG))[region_normal_coverage[,i+2]==0]
            CG[CG_0X]<-NaN
            plot_CG<-cbind(plot_CG,CG)
       }
        for (j in 1:length(tumor_samples)) {
            CG<-region_tumor_mc[,j+2]
            CG_0X<-(1:length(CG))[region_tumor_coverage[,j+2]==0]
            CG[CG_0X]<-NaN
            plot_CG<-cbind(plot_CG,CG)
       }


       #plot_CG<-cbind(region_normal_mc[,-2],region_tumor_mc[,-c(1,2)])
       plot(0, ylim=c(1,(length(tumor_samples)+length(normal_samples))), xlim=c(0.5, (dim(region_tumor_mc)[1]+0.5)), type="n", axes=FALSE, ylab="", xlab="position in bp", main=main)
       axis(2, at=c(1:(length(tumor_samples)+length(normal_samples))), labels=rev(c(tumor_samples, normal_samples)), cex.axis=1.3, las=2)
       axis(1, at=c(1:dim(plot_CG)[1]), labels=plot_CG[,2], cex.axis=0.5, las=2)

       # plot normal samples
       for (i in 1:length(normal_samples)) {
        order<-rev(1:length(normal_samples))
          abline(h=order[i], lwd=4 , col="magenta")
          for(j in 1:dim(plot_CG)[1]) {
             if (!is.na(plot_CG[j,i+2]) & as.numeric(plot_CG[j,i+2])==0) {
                  floating.pie(j,order[i], c(1), radius=pie.radius, col=c("white"), startpos=pi/2)
              } else if(!is.na(plot_CG[j,i+2])){
               floating.pie(j,order[i], c(plot_CG[j,i+2], 1-plot_CG[j,i+2]), radius=pie.radius, col=c("magenta", "white"), startpos=pi/2-plot_CG[j,i+2]*2*pi)
             }
           }  
       }

      # plot tumor samples
      for (i in (length(normal_samples)+1):(length(tumor_samples)+length(normal_samples))){
         order<-(length(tumor_samples)+length(normal_samples))-i+length(normal_samples)+1
          abline(h=order, lwd=4 , col="blue")
          for (j in 1:dim(plot_CG)[1]){
            if(!is.na(plot_CG[j,i+2]) & as.numeric(plot_CG[j,i+2])==0){
                floating.pie(j,order, c(1), radius=pie.radius, col=c("white"), startpos=pi/2)
            } else if(!is.na(plot_CG[j,i+2])){
                floating.pie(j,order, c(plot_CG[j,i+2], 1-plot_CG[j,i+2]), radius=pie.radius, col=c("blue", "white"), startpos=pi/2-plot_CG[j,i+2]*2*pi)
            }
         }
     }

     dev.off()
}


