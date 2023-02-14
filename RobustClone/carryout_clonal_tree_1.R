#lib = "/gpfs/research/fangroup/rk18g/longitudinal/RobustClone/example"
#install.packages("igraph")
#install.packages("shape")
#install.packages("vegan")
#install.packages("RANN")
#install.packages("gmodels")
#install.packages("sva")
#install.packages("reshape")
#install.packages("mgcv")
#install.packages("nlme")
#install.packages("R.matlab")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("sva")

library(igraph)
library(shape)
library(vegan)
library(RANN)
library(gmodels)
library(reshape)
library(sva)
library(mgcv)
library(nlme)
library(MASS)
options(max.print=1000000)
library(R.matlab)

setwd('/gpfs/research/fangroup/rk18g/longitudinal/clustering/RobustClone')
source('/gpfs/research/fangroup/rk18g/longitudinal/clustering/RobustClone/matlab_and_R_scripts/Clustering_EvolutionaryTree_function.R')

args = commandArgs(trailingOnly=TRUE)

start_time <- Sys.time()
## example for SNV data
#AA1_example <- R.matlab:::readMat('/gpfs/research/fangroup/rk18g/longitudinal/RobustClone/example//example_RobustClone.mat')
AA1_example <- R.matlab::readMat(args[1])
AA <- AA1_example[[1]] # Read GTM recovered by RPCA model or extended RPCA model 

robust_clone <- LJClustering(AA) # Louvain-Jaccard clustering
#print(robust_clone)
#write.table(robust_clone,file="Mat_1.csv",row.names = FALSE)
#write.matrix(robust_clone,file="Mat_1.csv")
write.matrix(robust_clone,file=args[2])
clone_gety <- subclone_GTM(AA, robust_clone, 'SNV') # obtain subclonal GTM
#print(clone_gety[[1]])
#save(clone_gety,file = "~/Downloads/RobustClone/clone_gety.rdata")
#write.matrix(clone_gety,file="Mat_2.csv")
write.matrix(clone_gety,file=args[3])

end_time <- Sys.time()
diff_time <- end_time-start_time
print(diff_time)

#MST <- plot_MST(clone_gety, robust_clone, 'SNV', 'exampledata') # calculate and plot clonal MST
#clones_mt <- new_variant_site(clone_gety, MST, 'SNV') # obtain the variant SNV loci each subclone compared with its parent subclone
#print(clones_mt[[1]])
