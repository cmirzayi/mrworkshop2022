rm(list=ls(all=TRUE))

#install MRPRESSO
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
library(MendelianRandomization)
library(metafor)
#read in file of SNPs to exposure and outcome for lead
pl<-read.csv("leadihd21.csv",header=TRUE)
nrow(pl)

#run MR-Presso 
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta", SdOutcome = "se.outcome", SdExposure = "se",
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = pl, NbDistribution = 10000,  SignifThreshold = 0.05)

#identify the SNPs MR-Pressos said were outliers
pl[7,]$SNP


#lets have a look at the original forest plot again
#4 get forest plot with fixed effects, same as Mendelianrandomization IVW with fixed effects 
x<-pl$beta # beta for SNP to exposure
sigmax<-pl$se # its standard errors
y<-pl$beta.outcome # beta for SNP to outcome
sigmay<-pl$se.outcome # its standard errors
pl$Wald<-y/x #Wald estimate
pl$Waldvar<-(sigmay^2/x^2) # using Burgess's method
pl$lab<-paste(pl$SNP, pl$gene, sep=" ")
pl<-pl[pl$SNP!="rs550057",] #uncomment this row to exclude rs550057
dmres<-rma.uni(yi=pl$Wald, vi=pl$Waldvar, slab=pl$lab, method="FE")
dmres
forest(dmres, atransf=exp,xlab=" ", mlab="Ischemic heart disease (OR)", at=log(c(.5, 1,2)),xlim=c(-1.7,1.3),cex=.8)

#2 get estimates for ihd using Mendelian randomization package with fixed effects,random effects, weighted median and MR-Egger
MRInputObject <- mr_input(pl$beta, pl$se, pl$beta.outcome, pl$se.outcome)
mr_allmethods(MRInputObject, method = "all")
mr_conmix(MRInputObject)
mr_lasso(MRInputObject)

