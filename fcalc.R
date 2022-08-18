rm(list=ls(all=TRUE))
#setwd("fill in here") # the location of your MR files
#setwd("//Users//cms1//Desktop//MRSLIDES2020//") # all my files for MR are here
# read the snp to exposure file
vtas<-read.csv("leadihd21.csv",header=TRUE)
vtas[1,]
# calculate F-statistic for each SNP
vtas$f1<-(vtas$beta*vtas$beta)/(vtas$se*vtas$se)
vtas$f1
#take the mean
mean(vtas$f1)
