#remove everything from workspace
rm(list=ls(all=TRUE))
install.packages(c("MendelianRandomization", "metafor"))
#tell R to use the packages you installed
library(MendelianRandomization)
library(metafor)
#i read snps to exposure into vtas
vtas<-read.table("lead.txt", header=TRUE) #read the file you created and uploaded

#check it is the correct file with the correct column names
vtas[1,]
nrow(vtas)

#ii read snps to the outcome file 
sl<-read.csv("cad-extract.csv")
sl[1,]
#iii rename the columns in the outcome file to be the same for every file
snpto<-data.frame(SNP = sl$markername,
                  beta.outcome = sl$beta,
                  se.outcome = sl$se_dgc,
                  effect_allele.outcome = sl$effect_allele,
                  other_allele.outcome =sl$noneffect_allele,
                  eaf.outcome = sl$effect_allele_freq,
                  pval.outcome = sl$p_dgc,
                  chrpos.outcome = paste(sl$chr,":",sl$bp_hg19,sep="")
)
snpto[1,]

#iv	Merge SNP to exposure and SNP to outcome files
dat<-merge(snpto,vtas,by="SNP",all.x=TRUE) #merge SNP to exposure and SNP to outcome file
nrow(dat)
dat[1,]
#v Align the SNPs on the same effect allele for exposure and outcome
# by changing sign of beta.exposure if effect alleles do not match

for (i in 1:length(dat$beta.outcome))
{
  dat$swapt[i]<-as.numeric(identical(as.character(dat$effect_allele.outcome[i]),as.character(toupper(dat$effect_allele[i])))) #align t to outcome
}
dat$swapt[dat$swapt ==0] <- -1 
dat$beta<-dat$beta*dat$swapt #reverse beta for exposure t if effect alleles do not match

#have a look at the alignment, are there any palindromic SNPs? are they aligned correctly?
data.frame(dat$SNP,dat$effect_allele,dat$other_allele,dat$effect_allele.outcome,dat$other_allele.outcome,dat$swapt,
           dat$eaf,dat$eaf.outcome)


#vi get forest plot with fixed effects, same as Mendelianrandomization IVW with fixed effects 
x<-dat$beta # beta for SNP to exposure
sigmax<-dat$se # its standard errors
y<-dat$beta.outcome # beta for SNP to outcome
sigmay<-dat$se.outcome # its standard errors
dat$Wald<-y/x #Wald estimate
data.frame(dat$SNP,dat$Wald)
dat$Waldvar<-(sigmay^2/x^2) # First term of Feiller's theorem
dat$lab<-paste(dat$SNP, dat$gene, sep=" ")
#pl<-pl[pl$SNP!="rs550057",]
dmres<-rma.uni(yi=dat$Wald, vi=dat$Waldvar, slab=dat$lab, method="FE")
dmres
forest(dmres, atransf=exp,xlab=" ", mlab="Ischemic heart disease (OR)", at=log(c(.5, 1,2)),xlim=c(-1.7,1.3),cex=.8)

#vii get estimates for ihd using Mendelian randomization package with fixed effects
mr_ivw<-MendelianRandomization::mr_ivw
mr_median<-MendelianRandomization::mr_median
MRInputObject <- mr_input(dat$beta, dat$se, dat$beta.outcome, dat$se.outcome)
mr_ivw(MRInputObject,model="fixed")
results<-mr_ivw(MRInputObject,model="fixed")
#exponentiate to get OR and 95% CIs matching the forest plot
round(data.frame(Estimate=exp(mr_ivw(MRInputObject,model="fixed")@Estimate),
           LCI=exp(mr_ivw(MRInputObject,model="fixed")@Estimate-1.96*mr_ivw(MRInputObject,model="fixed")@StdError),
           UCI=exp(mr_ivw(MRInputObject,model="fixed")@Estimate+1.96*mr_ivw(MRInputObject,model="fixed")@StdError)),2)

#viii get MR estimates for ihd using random effects
mr_ivw(MRInputObject)
results<-mr_ivw(MRInputObject)
#exponentiate to get OR and 95% CIs matching the forest plot
round(data.frame(Estimate=exp(mr_ivw(MRInputObject,model="random")@Estimate),
                 LCI=exp(mr_ivw(MRInputObject,model="random")@Estimate-1.96*mr_ivw(MRInputObject,model="random")@StdError),
                 UCI=exp(mr_ivw(MRInputObject,model="random")@Estimate+1.96*mr_ivw(MRInputObject,model="random")@StdError)),2)


#ix save file for future use
write.csv(dat,"leadihd21")

#use MR base code to do the analysis
devtools::install_github('MRCIEU/TwoSampleMR')
library(TwoSampleMR)
ao <- available_outcomes()
exposure_dat <- read_exposure_data(
  filename = 'lead.txt',
  sep = ' ',
  snp_col = 'SNP',
  beta_col = 'beta',
  se_col = 'se',
  effect_allele_col = 'effect_allele',
  phenotype_col = 'Phenotype',
  units_col = 'units',
  other_allele_col = 'other_allele',
  eaf_col = 'eaf',
  samplesize_col = 'samplesize',
  ncase_col = 'ncase',
  ncontrol_col = 'ncontrol',
  gene_col = 'gene',
  pval_col = 'pval'
)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('ieu-a-7'), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results <- mr(dat)
mr_results

