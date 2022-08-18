#remove everything from workspace
rm(list=ls(all=TRUE))
#use MR-Base generated code
#3.ii
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
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
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('ebi-a-GCST006867'), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr(dat)