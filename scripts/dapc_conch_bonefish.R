library(adegenet)
library(genepopedit)
library(diveRsity)
library(phangorn)
library(ape)
library(gtools)
library(pegas)

#######################################################################
###################### CUSTOM FUNCTIONS ###############################
#######################################################################

read_structure <- function(x){
  d<-read.table(x, header=F, sep="\t")
  dat <- read.structure(x, 
                        onerowperind=F, 
                        n.ind=(nrow(d)/2), 
                        n.loc=(ncol(d)-2), 
                        col.lab=1,
                        row.marknames=0,
                        col.pop=2,
                        ask=FALSE,
                        quiet=TRUE)
  return(dat)
}

dapc_xval <- function(dat, prefix){
  x<-dat
  mat <- tab(x, NA.method="mean")
  grp<- x$pop
  cx<-xvalDapc(mat, grp, n.pca.max = 300, training.set = 0.9,
               result = "groupMean", center = TRUE, scale = FALSE,
               n.pca = NULL, n.rep = 30, xval.plot = TRUE)
  out<-paste0(prefix, "_xval.pdf")
  pdf(out)
  plot(as.numeric(names(cx[5][[1]])), unlist(cx[5]), main="DAPC Cross-Validation RMSE", xlab="PCA Axes Retained", ylab="RMSE")
  dev.off()
  
  pc_retain<-as.numeric(as.vector(names(cx[5][[1]])))[which.min(unlist(cx[5]))]
  print(paste0("Num PC axes minimizing RMSE: ", pc_retain))
  
  dapc1 <- dapc(x, n.pca=pc_retain, n.da=length(unique(grp)))
  out2<-paste0(prefix, "_dapc.pdf")
  pdf(out2)
  scatter(dapc1, cex=2, cstar=0, scree.da=F)
  scatter(dapc1, xax=2, yax=3, cex=2, cstar=0, scree.da=F)
  scatter(dapc1, cex=2, cstar=1, scree.da=F)
  scatter(dapc1, xax=2, yax=3, cex=2, cstar=1, scree.da=F)
  dev.off()
  
  eig<-dapc1$eig
  pov<-eig/sum(eig)*100
  
  print(paste0("Percentage of variance explained by DA axes:"))
  print(pov)
  #data
  return(dapc1)
}


#######################################################################
######################## CONCH ANALYSIS ###############################
#######################################################################

setwd("~/Dropbox/My Mac (ARSC-7042066.local)/Documents/projects/conch")

#10 axes
#44% PCA variance
dat_sd2 <- read_structure("subsets/conch_sd2.str")
dapc_sd2 <- dapc_xval(dat_sd2, "~/Dropbox/Academic/Manuscripts/Ongoing_Collabs/Conch_Bonefish/conch_sd2")

#5 axes
#50.5% PCA variance
dat_sd3 <- read_structure("subsets/conch_sd3.str")
dapc_sd3 <- dapc_xval(dat_sd3, "~/Dropbox/Academic/Manuscripts/Ongoing_Collabs/Conch_Bonefish/conch_sd3")

#15 axes
#16.7% conserved variance
dat_full <- read_structure("subsets/conch_pop.str")
dapc_full <- dapc_xval(dat_full, "~/Dropbox/Academic/Manuscripts/Ongoing_Collabs/Conch_Bonefish/conch_full")

#Results: In all cases 5 PC axes minimized RMSE

#######################################################################
######################## BONEFISH ANALYSIS ############################
#######################################################################

setwd("~/Dropbox/My Mac (ARSC-7042066.local)/Documents/projects/bonefish")

#1=EFBC; 2=EOHS; 3=EOPC; 4=EOSS; 5=GBFN; 6=GBOS
#  EW       EE      ES1    ES2     GN      GS

#25 axes
#76.9% conserved variance
#35.917301 24.334526 21.185604 12.910580  5.651988
adat_sd2 <- read_structure("alb_sd2.str")
adapc_sd2 <- dapc_xval(adat_sd2, "~/Dropbox/Academic/Manuscripts/Ongoing_Collabs/Conch_Bonefish/alb_sd2")

#15 axes
#66% conserved variance
#39.267541 25.783931 18.944654 12.760828  3.243046
adat_sd3 <- read_structure("alb_sd3.str")
adapc_sd3 <- dapc_xval(adat_sd3, "~/Dropbox/Academic/Manuscripts/Ongoing_Collabs/Conch_Bonefish/alb_sd3")

#45 axes
#73% conserved variance
#41.572933 26.987107 16.728300  7.682078  7.029582
adat_full <- read_structure("alb_pop.str")
adapc_full <- dapc_xval(adat_full, "~/Dropbox/Academic/Manuscripts/Ongoing_Collabs/Conch_Bonefish/alb_full")

#Result
