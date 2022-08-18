install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)

exposure_dat <- extract_instruments("ieu-a-89") #overall height from a non-UK Biobank study
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST006697", maf_threshold = 0.42) #overall lifespan
datm <- harmonise_data(exposure_dat, outcome_dat) # Harmonise the exposure and outcome data, i.e., align the alleles
mr(datm) # perform MR

#also get estimates from MendelianRandomization
library(MendelianRandomization) #make sure each package is using the right function
mr_ivw<-MendelianRandomization::mr_ivw
mr_median<-MendelianRandomization::mr_median
datm<-datm[datm$mr_keep=="TRUE",] #only select the correctly aligned SNPs
#get estimates in years of life gained from Mendelianrandomization
MRInputObject <- mr_input(datm$beta.exposure, datm$se.exposure, datm$beta.outcome*-20, datm$se.outcome*20)
mr_ivw(MRInputObject) #ivw 
mr_median(MRInputObject) # weighted median
mr_egger(MRInputObject) # mr egger
mr_conmix(MRInputObject)  # conmix
