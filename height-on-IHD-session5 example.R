install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)

exposure_dat <- extract_instruments("ieu-a-89") #overall height from non-UK Biobank
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ukb-d-I9_IHD") #IHD
datm <- harmonise_data(exposure_dat, outcome_dat) # Harmonise the exposure and outcome data, i.e. align the alleles
mr(datm) # perform MR

library(MendelianRandomization) #make sure each package is using the right function
mr_ivw<-MendelianRandomization::mr_ivw
mr_median<-MendelianRandomization::mr_median
datm<-datm[datm$mr_keep=="TRUE",] #only select the correctly aligned SNPs
#convert estimates from probability to logodds for outcome
k<-20857/(20857 +340337) # cases to total copied from MR base
datm$beta.outcome<-log(1+datm$beta.outcome/(k*(1-k))) #make outcome beta into logodds
datm$se.outcome<-abs(datm$beta.outcome/qnorm(datm$pval.outcome/2)) #get se from the original p-value 

MRInputObject <- mr_input(datm$beta.exposure, datm$se.exposure, datm$beta.outcome, datm$se.outcome) #convert to years gained
mr_ivw(MRInputObject) #ivw 
mr_median(MRInputObject) # weighted median
mr_egger(MRInputObject) # mr-egger
mr_conmix(MRInputObject) # mr-egger
