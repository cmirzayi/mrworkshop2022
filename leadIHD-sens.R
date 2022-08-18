#remove everything from workspace
rm(list=ls(all=TRUE))
# tell R to use the packages you installed
library(MendelianRandomization)
library(metafor)


#1 read in the leadihd SNPs file you saved from lab 1 
pl<-read.csv("leadihd21.csv",header=TRUE)
nrow(pl)
#2 get estimates for ihd using Mendelian randomization package with fixed effects,random effects, weighted median and MR-Egger
MRInputObject <- mr_input(pl$beta, pl$se, pl$beta.outcome, pl$se.outcome)
mr_ivw(MRInputObject,model="fixed")
mr_ivw(MRInputObject)
mr_median(MRInputObject)
mr_egger(MRInputObject)
#nb exponentiate to get OR and 95% CIs, for OR is exp(Estimate), for CI its exp(Estimate + or - 1.96*SE)


