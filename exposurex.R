rm(list=ls(all=TRUE))
#i install MRbase items
install.packages("devtools")
library(devtools)
install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(readr)
library(MendelianRandomization)
# read in SNP to exposure GWAS downloaded from Long et al 2017
#untar("long_et_al_associations.tar.gz", files="Long_et_al_associations/gwas.summary.1E-5.tsv")
#untar("long_et_al_associations.tar", files="Long_et_al_associations/gwas.summary.1E-5.tsv")
#ii read in an extract of SNPs to exposure concerning theobromine only
long<-read.csv("long-theobromine.txt",header=TRUE)

#iii select SNPs with p-value<5x10-6
vk<-long[long$BIOCHEMICAL=="theobromine",] #theobromine
vkf<-vk[vk$pvalue<0.000005,] #select pvalue<5x10-6
nrow(vkf)

#iv only keep uncorrelated SNPs
#put data into format for MRbase
vta<-data.frame(SNP = vkf$rsid,
                beta.exposure = vkf$EffectSize,
                se.exposure = vkf$EffectSizeSE,
                effect_allele.exposure = vkf$alt,
                other_allele.exposure =vkf$ref,
                eaf.exposure = vkf$AltFreq,
                pval.exposure = vkf$pvalue,
                units.exposure = "Units", #not figured out yet
                gene.exposure = "vkf$gene",
                exposure = "theobromine",
                samplesize.exposure = vkf$num_GWAS_genomes)
#iv use MRbase to find indepnednt SNPs with r2<0.05
vtas<- clump_data(vta,clump_r2=0.05)

nrow(vtas)

#v write out a CSV file for future reference
write.csv(vtas,"theobro.csv")

#2a use MR base to get the estimates
outcome_datm <- extract_outcome_data(vtas$SNP, c('ieu-a-7'),) 
datm<- harmonise_data(vtas, outcome_datm) 
zz<-mr(datm)
zz

#2b put them into Mendelian randomization to get extra diagnostics
datm<-datm[datm$mr_keep=="TRUE",]
mr_ivw<-MendelianRandomization::mr_ivw
mr_median<-MendelianRandomization::mr_median
MRInputObject <- mr_input(datm$beta.exposure, datm$se.exposure, datm$beta.outcome, datm$se.outcome)
mr_ivw(MRInputObject) #ivw #remebre
mr_median(MRInputObject)
mr_egger(MRInputObject)
