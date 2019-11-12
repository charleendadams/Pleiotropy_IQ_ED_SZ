#########################################################################
##### Univariable and multivariable MR of intelligence, education, and SZ  
#### and Univariable and multivariable MR of education, bipolar disorder, and SZ
#########################################################################
#Note on r2 and F-stat
#########################################################################
# where:
# b = effect of SNP on exposure
# p=minor allele frequency 
# var = variance of phenotype 
# n= sample size 
# k = number of SNPs

# To obtain standardized betas:
# z = b/se
# b = z/sqrt(2*p*(1-p)*(n+z^2))

# r2 <- (2*b^2)*p*(1-p)/var #var is 1 for exposures in SD units

# When Steiger and hand-calculated differ, use the smaller r2 
# Fform <- r2*(n-1-k)/((1-r2)*k) 
#########################################################################
##### Univariable MR of UK Biobank IQ on SZ
#########################################################################
setwd("C:/Users/charl/Dropbox/COH/IQ_school_sz")
library(TwoSampleMR)
library(MendelianRandomization)
library(RadialMR)
library(MRInstruments)
ao <- available_outcomes()
exposure_dat <- extract_instruments(c('UKB-a:196'))
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('22'), 
  proxies = 1, rsq = 0.8, 
  align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results <- mr(dat)

#Outlier's check
mrinput <- mr_input(bx = dat$beta.exposure, 
  bxse = dat$se.exposure, by = dat$beta.outcome, byse=dat$se.outcome)
mregger <- mr_egger(mrinput, penalized=F) 
raddat <- format_radial(BXG=dat$beta.exposure, BYG=dat$beta.outcome, 
  seBXG=dat$se.exposure, seBYG=dat$se.outcome, RSID=dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=1)
ivwrad2 <- ivw_radial(raddat, alpha=0.05, weights=1)
dim(ivwrad$outliers)[1] 
outliers=ivwrad$outliers[1]
outliers
myvars=outliers$SNP
dat2 <- dat[ ! dat$SNP %in% myvars, ]
dim(dat2) ## 22 SNPs
dim(dat)
mr_results <- mr(dat2)
mr_results=generate_odds_ratios(mr_results)
singlesnp_results=mr_singlesnp(dat2, parameters = default_parameters(),
  single_method = "mr_wald_ratio", 
  all_method = c("mr_ivw",
  "mr_egger_regression", 
  'mr_weighted_mode',
  'mr_weighted_median'))
singlesnp_results=generate_odds_ratios(singlesnp_results)
write.csv(mr_results,'res_OR_IQ_SZ.csv')
write.csv(singlesnp_results,'single_snps_OR_IQ_SZ.csv')
mr_forest_plot(singlesnp_results, exponentiate = TRUE)
mr_scatter_plot(mr_results, dat2)

#Remove the meta-analyses results to obtain the list of SNPs
#used in the MR models (some were dropped from dat2)
dim(dat2)
singlesnp_results.2=singlesnp_results[1:17,]
dim(singlesnp_results.2)
dat2.1=dat2[which(dat2$SNP %in% singlesnp_results.2$SNP),]
dim(dat2.1)
dim(dat2)
dat2=dat2.1
write.csv(dat2,'dat2_17snps_IQ_SZ_SIMEX.csv') #for Supplement

#Generate r2 and Fstat
p=dat2$eaf.exposure 
n=82315 #sample size -- smaller of the two
cases=35476 #cases
cases/n
k=17 #number of SNPs
out <- directionality_test(dat2)
r2=out$snp_r2.exposure #0.005336977
r2

#IQ is normalized -- so the variance is 1 and I don't have to transform the beta
#Check the difference between the 'directionality' test r2 and the manual r2
#use the smaller
r2.1<-(2*(dat2$beta.exposure^2)*p*(1-p))/1 
r2.manual=sum(r2.1)
r2.manual # 0.02402042
Fstat <- r2*(n-1-k)/((1-r2)*k) # 26.47218
Fstat

#Combine all results 
#Already log-transformed above
sin=singlesnp_results
res=mr_results
het<-mr_heterogeneity(dat2)
plt<-mr_pleiotropy_test(dat2)
#sin<-mr_singlesnp(dat2)
all_res<-combine_all_mrresults(res,het,plt,sin,ao_slc=T,
  Exp=T,split.exposure=T,split.outcome=T)
all_res$r2=r2
all_res$lower_intercept=all_res$intercept-1.96*all_res$intercept_se
all_res$upper_intercept=all_res$intercept+1.96*all_res$intercept_se
all_res$Fstat=Fstat
all_res$I2=66.58 #from stata

head(all_res[,c("Method","outcome","exposure","nsnp","b","se","or", "or_lci95",  
 "or_uci95","pval","intercept","intercept_se","intercept_pval",
 "Q","Q_df","Q_pval","consortium","ncase","ncontrol","pmid","population")])
write.csv(all_res,"combined_results_17_snps_IQ_SZ.csv")

#########################################################################
##### Univariable MR of EduAge on SZ
#########################################################################
setwd("C:/Users/charl/Dropbox/COH/IQ_school_sz/ed_SZ")
library(TwoSampleMR)
library(ggplot2)
library(MendelianRandomization)
library(RadialMR)
library(MRInstruments)
ao <- available_outcomes()
exposure_dat <- extract_instruments(c('UKB-a:505'))
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('22'), 
  proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results <- mr(dat)

#Outlier's check
mrinput <- mr_input(bx = dat$beta.exposure, 
  bxse = dat$se.exposure, by = dat$beta.outcome, byse=dat$se.outcome)
mregger <- mr_egger(mrinput, penalized=F) 
raddat <- format_radial(BXG=dat$beta.exposure, BYG=dat$beta.outcome, 
  seBXG=dat$se.exposure, seBYG=dat$se.outcome, RSID=dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=1)
ivwrad2 <- ivw_radial(raddat, alpha=0.05, weights=1)
dim(ivwrad$outliers)[1] 
outliers=ivwrad$outliers[1]
outliers
myvars=outliers$SNP
dat2 <- dat[ ! dat$SNP %in% myvars, ]
dim(dat2)
dim(dat)

#Remove the EduAge SNPs in LD with the IQ SNPs
#Determined from the MV model
dat$SNP_ld=c('rs6736898','rs4731951')
myvars2=dat$SNP_ld
dat_ld <- dat2[! dat2$SNP %in% myvars2, ]
dim(dat_ld)
dat2=dat_ld
mr_results <- mr(dat2)
mr_results=generate_odds_ratios(mr_results)
singlesnp_results=mr_singlesnp(dat2, parameters = default_parameters(),
  single_method = "mr_wald_ratio", 
  all_method = c("mr_ivw",
  "mr_egger_regression", 
  'mr_weighted_mode',
  'mr_weighted_median'))
singlesnp_results=generate_odds_ratios(singlesnp_results)
write.csv(singlesnp_results,'singlesnp_results_edu_SZ.csv')
write.csv(mr_results,'res_ed_SZ.csv')
mr_forest_plot(singlesnp_results, exponentiate = TRUE)
mr_scatter_plot(mr_results, dat2)

#Remove the meta-analyses results to obtain the list of SNPs
#used in the MR models (some were dropped from dat2)
singlesnp_results.2=singlesnp_results[1:9,]
dim(singlesnp_results.2)
dat2.1=dat2[which(dat2$SNP %in% singlesnp_results.2$SNP),]
dim(dat2.1)
dim(dat2)
dat2=dat2.1
write.csv(dat2,'dat2_9snps_EDU_AGE_SZ_SIMEX.csv')

#Generate r2 and Fstat
n= 82315 #sample size of the smaller
cases= 35476 #cases
controls=cases/n
k=9 #number SNPs
p=dat2$eaf.exposure

#Check the difference between the 'directionality' test r2 and the manual r2
#They differed, so I checked a 3rd way with the 'get_r_from_lor' function
#used the LARGER in this case, 
out <- directionality_test(dat2)
r2=out$snp_r2.exposure
r2 #0.001448087
r2.1<-(2*(dat2$beta.exposure^2)*p*(1-p))/1 #
r2.manual=sum(r2.1)
r2.manual #0.000984903
lor=dat2$beta.outcome
af=p
ncase=cases
ncontrol=controls
prevalence = 0.01
r3=get_r_from_lor(lor, af, ncase, ncontrol, prevalence, model="logit", correction=FALSE)
r3.1=sum(r3)
Fstat.r3.1 <- r3.1*(n-1-k)/((1-r3.1)*k) #286.3027; seems too big
Fstat.r.manual <- r2.manual*(n-1-k)/((1-r2.manual)*k) 
Fstat.r.manual #9.015818
Fstat <- r2*(n-1-k)/((1-r2)*k) # 13.26196; seems reasonable 

#Combine results
sin=singlesnp_results #already log-transformed above
res=mr_results
het<-mr_heterogeneity(dat2)
plt<-mr_pleiotropy_test(dat2)
all_res<-combine_all_mrresults(res,het,plt,sin,ao_slc=T,
  Exp=T,split.exposure=T,split.outcome=T)
all_res$r2=r2.manual
all_res$lower_intercept=all_res$intercept-1.96*all_res$intercept_se
all_res$upper_intercept=all_res$intercept+1.96*all_res$intercept_se
all_res$Fstat #22.11383
all_res$I2=79.86 #stata

head(all_res[,c("Method","outcome","exposure","nsnp","b","se","or", "or_lci95",  
  "or_uci95","pval","intercept","intercept_se","intercept_pval",
  "Q","Q_df","Q_pval","consortium","ncase","ncontrol","pmid","population")])
write.csv(all_res,"combined_results_9_snps_LD_clean_EDU_SZ.csv")

#########################################################################
#### Univariable MR of Lee EduYears on SZ
#########################################################################
setwd("C:/Users/charl/Dropbox/COH/IQ_school_sz/Edu_Years_GWAS/edu_years_sz")
library(ggplot2)
library(MendelianRandomization)
library(RadialMR)
library(TwoSampleMR)
ao <- available_outcomes()
exposure_dat <- read_exposure_data(
  filename = 'edu_years_GWAS.csv',
  sep = ',',
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
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('22'), 
  proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results <- mr(dat)

#Outlier's check
mrinput <- mr_input(bx = dat$beta.exposure, 
  bxse = dat$se.exposure, by = dat$beta.outcome, byse=dat$se.outcome)
mregger <- mr_egger(mrinput, penalized=F) 
raddat <- format_radial(BXG=dat$beta.exposure, BYG=dat$beta.outcome, 
  seBXG=dat$se.exposure, seBYG=dat$se.outcome, RSID=dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=1)
ivwrad2 <- ivw_radial(raddat, alpha=0.05, weights=1)
dim(ivwrad$outliers)[1] 
outliers=ivwrad$outliers[1]
outliers
myvars=outliers$SNP
dat2 <- dat[ ! dat$SNP %in% myvars, ]
dim(dat2)
dim(dat)
write.csv(dat2, "EduYears_outliers_gone_277_dat.csv")
mr_results <- mr(dat2)
write.csv(dat2,'dat2_edu_years_SZ.csv')
mr_results=generate_odds_ratios(mr_results)
singlesnp_results=mr_singlesnp(dat2, parameters = default_parameters(),
  single_method = "mr_wald_ratio", 
  all_method = c("mr_ivw",
  "mr_egger_regression", 
  'mr_weighted_mode',
  'mr_weighted_median'))
singlesnp_results=generate_odds_ratios(singlesnp_results)
write.csv(mr_results,'res_OR_Edu_years_SZ.csv')
write.csv(singlesnp_results,'single_snps_OR_Edu_years_SZ.csv')
mr_forest_plot(singlesnp_results, exponentiate = TRUE)
mr_scatter_plot(mr_results, dat2)

#remove the snps that weren't used in the MR
singlesnp_results.2=singlesnp_results[1:238,]
dim(singlesnp_results.2)
dat2.1=dat2[which(dat2$SNP %in% singlesnp_results.2$SNP),]
dim(dat2.1)
dim(dat2)
dat2=dat2.1
write.csv(dat2, '238_EduYears_snps.csv')

#r2 and Fstat
n= 82315#sample size
cases= 35476#cases
cases/n
k=238
p=dat2$eaf.exposure
r2.1<-(2*(dat2$beta.exposure^2)*p*(1-p))/1
r2=sum(r2.1)
r2 #0.01348805
Fstat <- r2*(n-1-k)/((1-r2)*k)
Fstat 
res_raps <- mr(dat2, method_list = c("mr_raps"), 
  parameters = list(over.dispersion = FALSE, loss.function = "l2"))
res_raps

#Combine results
sin=generate_odds_ratios(singlesnp_results)
res=mr_results
het<-mr_heterogeneity(dat2)
plt<-mr_pleiotropy_test(dat2)
#sin<-mr_singlesnp(dat2)
all_res<-combine_all_mrresults(res,het,plt,sin,ao_slc=T,
  Exp=T,split.exposure=T,split.outcome=T)
all_res$r2=r2
all_res$lower_intercept=all_res$intercept-1.96*all_res$intercept_se
all_res$upper_intercept=all_res$intercept+1.96*all_res$intercept_se
all_res$Fstat=Fstat
all_res$I2=77.28 #from stata
head(all_res[,c("Method","outcome","exposure","nsnp","b","se","or", "or_lci95",  
  "or_uci95","pval","intercept","intercept_se","intercept_pval",
  "Q","Q_df","Q_pval","consortium","ncase","ncontrol","pmid","population")])
write.csv(all_res,"combined_results_Lee_Edu_years_SZ.csv")

#########################################################################
#####  Univariable MR of Hill IQ on SZ
#########################################################################
setwd("C:/Users/charl/Dropbox/COH/IQ_school_sz/Edu_Years_GWAS/edu_years_IQ")#where the IQ data are
#Clean the Hill IQ GWAS
IQ=read.csv('IQ-187-dat.csv')
dim(IQ)
IQ$Phenotype="IQ"
write.csv(IQ, "IQ.csv")
#IQ.clean=IQ[!duplicated(IQ$SNP), ]
#dim(IQ.clean)
#head(IQ.clean[,1:5])
library(MendelianRandomization)
library(RadialMR)
library(TwoSampleMR)
ao <- available_outcomes()
exposure_dat <- read_exposure_data(
  filename = 'IQ.csv',
  sep = ',',
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
exposure_dat <- clump_data(exposure_dat)
dim(exposure_dat)
setwd("C:/Users/charl/Dropbox/COH/IQ_school_sz/Edu_Years_GWAS/edu_years_sz")#change back 
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('22'), 
  proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results <- mr(dat)

#Outlier's check
mrinput <- mr_input(bx = dat$beta.exposure, 
  bxse = dat$se.exposure, by = dat$beta.outcome, byse=dat$se.outcome)
mregger <- mr_egger(mrinput, penalized=F) 
raddat <- format_radial(BXG=dat$beta.exposure, BYG=dat$beta.outcome, 
  seBXG=dat$se.exposure, seBYG=dat$se.outcome, RSID=dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=1)
ivwrad2 <- ivw_radial(raddat, alpha=0.05, weights=1)
dim(ivwrad$outliers)[1] 
outliers=ivwrad$outliers[1]
outliers
myvars=outliers$SNP
dat2 <- dat[ ! dat$SNP %in% myvars, ]
dim(dat2) 
dim(dat)
write.csv(dat2, "Hill-SZ-outliers-gone.dat.csv")
mr_results <- mr(dat2)
write.csv(mr_results,'Hill_SZ-outliers-gone_results.csv')
singlesnp_results=mr_singlesnp(dat2, parameters = default_parameters(),
  single_method = "mr_wald_ratio", 
  all_method = c("mr_ivw",
  "mr_egger_regression", 
  'mr_weighted_mode',
  'mr_weighted_median'))
singlesnp_results=generate_odds_ratios(singlesnp_results)
write.csv(singlesnp_results,'Hill_SZ-outliers-gone_singlesnp.csv')

#remove the snps that weren't used in the MR
singlesnp_results.2=singlesnp_results[1:103,]
dim(singlesnp_results.2)
dat2.1=dat2[which(dat2$SNP %in% singlesnp_results.2$SNP),]
dim(dat2.1)
dim(dat2)
dat2=dat2.1
write.csv(dat2, '103_Hill_SZ-snps.csv')

#Generate r2 and Fstat
n= 82315#sample size
cases= 35476#cases
cases/n
n=82315
k=103
p=dat2$eaf.exposure
dat2$samplesize.exposure=293723
r2.1<-(2*(dat2$beta.exposure^2)*p*(1-p))/1
r2=sum(r2.1)
r2 #0.01838603
Fstat <- r2*(n-1-k)/((1-r2)*k)
Fstat #14.94996

#Combine results
sin=singlesnp_results
res=mr_results
het<-mr_heterogeneity(dat2)
plt<-mr_pleiotropy_test(dat2)
all_res<-combine_all_mrresults(res,het,plt,sin,ao_slc=T,
  Exp=T,split.exposure=T,split.outcome=T)
all_res$r2=r2
all_res$lower_intercept=all_res$intercept-1.96*all_res$intercept_se
all_res$upper_intercept=all_res$intercept+1.96*all_res$intercept_se
all_res$Fstat=Fstat
all_res$I2=51.65 #from stata 
head(all_res[,c("Method","outcome","exposure","nsnp","b","se","or", "or_lci95",  
  "or_uci95","pval","intercept","intercept_se","intercept_pval",
  "Q","Q_df","Q_pval","consortium","ncase","ncontrol","pmid","population")])
write.csv(all_res,"combined_results_Hill_SZ.csv")

#########################################################################
##### Multivariable MR of intelligence and EduAge on SZ
#########################################################################
setwd("C:/Users/charl/Dropbox/COH/IQ_school_sz/multi")
library(ggplot2)
library(MendelianRandomization)
library(RadialMR)
library(MRInstruments)
ao <- available_outcomes()
id_exposure=c('UKB-a:505','UKB-a:196')
id_outcome=c('22')
exposure_dat <- mv_extract_exposures(id_exposure)
#remove the snps not in the main models
iq=read.csv("C:/Users/charl/Dropbox/COH/IQ_school_sz/combined_results_17_snps_IQ_SZ.csv")
head(iq)
myvars=iq$SNP
edu=read.csv("C:/Users/charl/Dropbox/COH/IQ_school_sz/combined_results_11_snps_EDU_SZ.csv")
myvars2=edu$SNP
dat3 <- exposure_dat[ exposure_dat$SNP %in% myvars | exposure_dat$SNP %in% myvars2, ]
dim(dat3)
head(dat3)
dim(exposure_dat)
exposure_dat=dat3
dim(exposure_dat)
write.csv(dat3, 'kept_multi_MR_IQ_EduAge_SZ.csv')

#Next, also extract those SNPs from the outcome.
outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome)

#Once the data has been obtained, harmonise so that all are on the same reference allele.
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)

#Finally, perform the multivariable MR analysis
res <- mv_multiple(mvdat)
res
write.csv(res,'multi_res.csv')

#########################################################################
#### Univariable MR of bipolar disorder on SZ
#########################################################################
exposure_dat <- extract_instruments(c('801'))
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('22'), 
  proxies = 1, rsq = 0.8, 
  align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results <- mr(dat)

#Outlier's check
mrinput <- mr_input(bx = dat$beta.exposure, 
 bxse = dat$se.exposure, 
 by = dat$beta.outcome, byse=dat$se.outcome)
mregger <- mr_egger(mrinput, penalized=F) 
raddat <- format_radial(BXG=dat$beta.exposure, BYG=dat$beta.outcome, 
  seBXG=dat$se.exposure, seBYG=dat$se.outcome, RSID=dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=1)
ivwrad2 <- ivw_radial(raddat, alpha=0.05, weights=1)
dim(ivwrad$outliers)[1] 
outliers=ivwrad$outliers[1]
outliers
myvars=outliers$SNP
dat2 <- dat[ ! dat$SNP %in% myvars, ]
dim(dat2) ## 22 SNPs
dim(dat)
mr_results <- mr(dat2)
mr_results=generate_odds_ratios(mr_results)
singlesnp_results=mr_singlesnp(dat, parameters = default_parameters(),
  single_method = "mr_wald_ratio", 
  all_method = c("mr_ivw",
  "mr_egger_regression", 
  'mr_weighted_mode',
  'mr_weighted_median'))
singlesnp_results=generate_odds_ratios(singlesnp_results)
write.csv(mr_results,'res_OR_IQ_SZ.csv')
write.csv(singlesnp_results,'single_snps_OR_bip.csv')
mr_forest_plot(singlesnp_results, exponentiate = TRUE)
mr_scatter_plot(mr_results, dat)

#Remove the meta-analyses results to obtain the list of SNPs
#used in the MR models (some were dropped from dat2)
dim(dat)
singlesnp_results.2=singlesnp_results[1:4,]
dim(singlesnp_results.2)
dat2.1=dat[which(dat$SNP %in% singlesnp_results.2$SNP),]
dim(dat2.1)
dim(dat)
dat2=dat2.1
write.csv(dat2,'dat2_bip_SIMEX.csv') #for Supplement

#Generate r2 and Fstat
p=dat2$eaf.exposure 
n=16731  #sample size -- smaller of the two
cases=7481  #cases
cases/n
k=4 #number of SNPs
out <- directionality_test(dat2)
r2=out$snp_r2.exposure #0.008172205
r2
Fstat <- r2*(n-1-k)/((1-r2)*k) # 34.45364
Fstat

#Combine all results 
#Already log-transformed above
sin=singlesnp_results
res=mr_results
het<-mr_heterogeneity(dat2)
plt<-mr_pleiotropy_test(dat2)
#sin<-mr_singlesnp(dat2)
all_res<-combine_all_mrresults(res,het,plt,sin,ao_slc=T,
  Exp=T,split.exposure=T,split.outcome=T)
all_res$r2=r2
all_res$lower_intercept=all_res$intercept-1.96*all_res$intercept_se
all_res$upper_intercept=all_res$intercept+1.96*all_res$intercept_se
all_res$Fstat=Fstat
all_res$I2=69.15 #from stata
head(all_res[,c("Method","outcome","exposure","nsnp","b","se","or", "or_lci95",  
 "or_uci95","pval","intercept","intercept_se","intercept_pval",
  "Q","Q_df","Q_pval","consortium","ncase","ncontrol","pmid","population")])
write.csv(all_res,"combined_results_bip.csv")

#########################################################################
#### Multivariable MR of Okbay EduYears and IQ on SZ
#########################################################################
setwd("C:/Users/charl/Dropbox/COH/IQ_school_sz/Edu_Years_GWAS/edu_years_sz")
id_exposure=c('1001','UKB-a:196')
id_outcome=c('22')
exposure_dat <- mv_extract_exposures(id_exposure)
dim(exposure_dat)
exposure_dat$original.id.exposure=exposure_dat$id.exposure
exposure_dat$id.exposure='joint_1239_196'
exposure_dat2 <- clump_data(exposure_dat)
dim()
exposure_dat$id.exposure=exposure_dat$original.id.exposure
write.csv(exposure_dat, 'mvmr_okbay_fluid_SZ.csv')

#look at the overlap in genome-wide sig p-values for the two and create a SNP list to exclude
#those SNPs that tag both
myvars_IQ=c(
  'rs7599488',
  'rs7029201',
  'rs1391438',
  'rs11191193',
  'rs11678980',
  'rs12410444',
  'rs28792186',
  'rs9739070')
dat3 <- exposure_dat[ !exposure_dat$SNP %in% myvars_IQ, ]
exposure_dat=dat3

#Remove the snps not in the main models
iq=read.csv("C:/Users/charl/Dropbox/COH/IQ_school_sz/dat2_17snps_IQ_SZ_SIMEX.csv")
head(iq)
myvars_IQ=iq$SNP
length(myvars_IQ)
dim(iq)
edu=read.csv("C:/Users/charl/Dropbox/COH/IQ_school_sz/Edu_Years_GWAS/edu_years_sz/EDU_YEARS_1001_45snps_for_SIMEX.csv")
myvars_EDU_YEARS=edu$SNP
length(myvars_EDU_YEARS)
dat4 <- exposure_dat[ exposure_dat$SNP %in% myvars_IQ| exposure_dat$SNP %in% myvars_EDU_YEARS, ]
dim(dat4)
exposure_dat=dat4

#Extract the outcome data
outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome)

#Harmonize
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)

#Perform MVMR
res <- mv_multiple(mvdat)
res
MV_IQ_EduYears_SZ='res'
MV_IQ_EduYears_SZ=as.data.frame(MV_IQ_EduYears_SZ)
MV_IQ_EduYears_SZ$IQ_SNP=13
MV_IQ_EduYears_SZ$IQ_OR=exp(-0.1809660)
MV_IQ_EduYears_SZ$IQ_lCI=exp(-0.1809660)-1.96*0.06337284
MV_IQ_EduYears_SZ$IQ_uCI=exp(-0.1809660)+1.96*0.06337284
MV_IQ_EduYears_SZ$IQ_pval=4.295870e-03
MV_IQ_EduYears_SZ$IQ_beta=-0.1809660
MV_IQ_EduYears_SZ$IQ_se=0.06337284
MV_IQ_EduYears_SZ$EduYears_SNP=36
MV_IQ_EduYears_SZ$EduYears_OR=exp( 0.6701502)
MV_IQ_EduYears_SZ$EduYears_lCI=exp( 0.6701502)-1.96*0.15912593
MV_IQ_EduYears_SZ$EduYears_uCI=exp( 0.6701502)+1.96*0.15912593
MV_IQ_EduYears_SZ$EduYears_pval=2.537414e-05
MV_IQ_EduYears_SZ$EduYears_beta=0.6701502
MV_IQ_EduYears_SZ$EduYears_se=0.15912593
dim(MV_IQ_EduYears_SZ)
MV_IQ_EduYears_SZ[,2:13]
write.csv(MV_IQ_EduYears_SZ, 'MV_IQ_Okbay_EduYears_SZ_snps_tagging_both_removed.csv')

#########################################################################
#### Multivariable MR of education and bipolar disorder on SZ
#########################################################################
setwd("C:/Users/charl/Dropbox/COH/IQ_school_sz/bip")
library(TwoSampleMR)
library(MendelianRandomization)
library(RadialMR)
library(MRInstruments)
ao <- available_outcomes()
id_exposure=c('1001','801')
id_outcome=c('22')
exposure_dat <- mv_extract_exposures(id_exposure)
exposure_dat2 <- clump_data(exposure_dat)
write.csv(exposure_dat, 'mvmr_SZ.csv')
#Extract the outcome data
outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome)
#Harmonize
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
#Perform MVMR
res <- mv_multiple(mvdat)
res
MV_bp_EduYears_SZ='res'
MV_bp_EduYears_SZ=as.data.frame(MV_bp_EduYears_SZ)
MV_bp_EduYears_SZ$bp_SNP=3
MV_bp_EduYears_SZ$bp_OR=exp(0.1473986)
MV_bp_EduYears_SZ$bp_lCI=exp(0.1473986)-1.96*0.06896995
MV_bp_EduYears_SZ$bp_uCI=exp(0.1473986)+1.96*0.06896995
MV_bp_EduYears_SZ$bp_pval=0.03258638
MV_bp_EduYears_SZ$bp_beta=0.1473986
MV_bp_EduYears_SZ$bp_se=0.06896995
MV_bp_EduYears_SZ$EduYears_SNP=51
MV_bp_EduYears_SZ$EduYears_OR=exp(0.2708158)
MV_bp_EduYears_SZ$EduYears_lCI=exp( 0.2708158)-1.96*0.21448209
MV_bp_EduYears_SZ$EduYears_uCI=exp( 0.2708158)+1.96*0.21448209
MV_bp_EduYears_SZ$EduYears_pval=0.20671504
MV_bp_EduYears_SZ$EduYears_beta=0.2708158
MV_bp_EduYears_SZ$EduYears_se=0.21448209
dim(MV_bp_EduYears_SZ)
MV_bp_EduYears_SZ[,2:15]
write.csv(MV_bp_EduYears_SZ, 'bip_mv.csv')

#########################################################################
###### Univariable MR of Okbay EduYears on SZ
#########################################################################
setwd("C:/Users/charl/Dropbox/COH/IQ_school_sz/Edu_Years_GWAS/edu_years_sz")
library(MendelianRandomization)
library(RadialMR)
library(TwoSampleMR)
ao <- available_outcomes()
exposure_dat <- extract_instruments(c('1001')) #1239 is the larger EDU_YEARs GWAS
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('22'), 
  proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results <- mr(dat)

#Outlier's check
mrinput <- mr_input(bx = dat$beta.exposure, 
  bxse = dat$se.exposure, by = dat$beta.outcome, byse=dat$se.outcome)
mregger <- mr_egger(mrinput, penalized=F) 
raddat <- format_radial(BXG=dat$beta.exposure, BYG=dat$beta.outcome, 
  seBXG=dat$se.exposure, seBYG=dat$se.outcome, RSID=dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=1)
ivwrad2 <- ivw_radial(raddat, alpha=0.05, weights=1)
dim(ivwrad$outliers)[1] 
outliers=ivwrad$outliers[1]
outliers
myvars=outliers$SNP
dat2 <- dat[ ! dat$SNP %in% myvars, ]
dim(dat2) ##  SNPs
dim(dat)
mr_results <- mr(dat2)
mr_results=generate_odds_ratios(mr_results)
dim(dat2)

singlesnp_results=mr_singlesnp(dat2, parameters = default_parameters(),
  single_method = "mr_wald_ratio", all_method = c("mr_ivw",
  "mr_egger_regression", 
  'mr_weighted_mode',
  'mr_weighted_median'))
singlesnp_results=generate_odds_ratios(singlesnp_results)
singlesnp_results.2=singlesnp_results[1:45,]
dim(singlesnp_results.2)
dat2.1=dat2[which(dat2$SNP %in% singlesnp_results.2$SNP),]
dim(dat2.1)
dim(dat2)
dat2=dat2.1
write.csv(dat2, 'EDU_YEARS_1001_45snps_for_SIMEX.csv')

#Generate r2 and Fstat
out <- directionality_test(dat2)
r2=out$snp_r2.exposure
r2 #0.006119558
n=82315 #smaller of the two
ncase=35476
ncontrol=n-ncase
k=45
p=dat2$eaf.exposure
Fstat <- r2*(n-1-k)/((1-r2)*k)
Fstat 

#Combine results
sin=singlesnp_results
res=mr_results
het<-mr_heterogeneity(dat2)
plt<-mr_pleiotropy_test(dat2)
all_res<-combine_all_mrresults(res,het,plt,sin,ao_slc=T,
  Exp=T,split.exposure=T,split.outcome=T)
all_res$r2=r2
all_res$lower_intercept=all_res$intercept-1.96*all_res$intercept_se
all_res$upper_intercept=all_res$intercept+1.96*all_res$intercept_se
all_res$Fstat=Fstat
all_res$I2= 47.25 #from stata 
head(all_res[,c("Method","outcome","exposure","nsnp","b","se","or", "or_lci95",  
  "or_uci95","pval","intercept","intercept_se","intercept_pval",
  "Q","Q_df","Q_pval","consortium","ncase","ncontrol","pmid","population")])
write.csv(all_res,"combined_results_Edu_years_1001_SZ.csv")

#########################################################################
##### Univariable MR of Hill intelligence on Lee EduYears
#### Bidirectional analysis of education and intelligence: part 1
#########################################################################
setwd("C:/Users/charl/Dropbox/COH/IQ_school_sz/Edu_Years_GWAS/edu_years_IQ")
#Clean the Hill IQ GWAS
IQ=read.csv('IQ-187-dat.csv')
dim(IQ)
IQ$Phenotype="IQ"
write.csv(IQ, "IQ.csv")
#IQ.clean=IQ[!duplicated(IQ$SNP), ]
#dim(IQ.clean)
#head(IQ.clean[,1:5])
library(MendelianRandomization)
library(RadialMR)
library(TwoSampleMR)
ao <- available_outcomes()
exposure_dat <- read_exposure_data(
  filename = 'IQ.csv',
  sep = ',',
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
exposure_dat <- clump_data(exposure_dat)
dim(exposure_dat)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('1239'), 
  proxies = 1, rsq = 0.8, align_alleles = 1, 
  palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results <- mr(dat)

#Outlier's check
mrinput <- mr_input(bx = dat$beta.exposure, 
  bxse = dat$se.exposure, by = dat$beta.outcome, byse=dat$se.outcome)
mregger <- mr_egger(mrinput, penalized=F) 
raddat <- format_radial(BXG=dat$beta.exposure, BYG=dat$beta.outcome, 
  seBXG=dat$se.exposure, seBYG=dat$se.outcome, RSID=dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=1)
ivwrad2 <- ivw_radial(raddat, alpha=0.05, weights=1)
dim(ivwrad$outliers)[1] 
outliers=ivwrad$outliers[1]
outliers
myvars=outliers$SNP
dat2 <- dat[ ! dat$SNP %in% myvars, ]
dim(dat2) ##  SNPs
dim(dat)
write.csv(dat2, "Hill_Lee-outliers-gone.dat.csv")
mr_results <- mr(dat2)
write.csv(mr_results,'Hill_Lee-outliers-gone_results.csv')
singlesnp_results=mr_singlesnp(dat2, 
  parameters = default_parameters(),
  single_method = "mr_wald_ratio", 
  all_method = c("mr_ivw",
  "mr_egger_regression", 
  'mr_weighted_mode',
  'mr_weighted_median'))
write.csv(singlesnp_results,'Hill_Lee-outliers-gone_singlesnp.csv')

#remove the snps that weren't used in the MR
singlesnp_results.2=singlesnp_results[1:101,]
dim(singlesnp_results.2)
dat2.1=dat2[which(dat2$SNP %in% singlesnp_results.2$SNP),]
dim(dat2.1)
dim(dat2)
dat2=dat2.1
write.csv(dat2, '101_Hill_IQ_snps.csv')

#determine LD
not_ld=read.csv('unique-hill-lee-346.csv')
dim(not_ld)
newdata <- not_ld[ which(not_ld$study=='hill'), ]
dim(newdata)
dat3=dat2[which(dat2$SNP %in% newdata$SNP),]
dim(dat3)
mr_results <- mr(dat3)
write.csv(mr_results,'Hill-noLD-Lee35snps.csv')
singlesnp_results=mr_singlesnp(dat3, 
  parameters = default_parameters(),
  single_method = "mr_wald_ratio", 
  all_method = c("mr_ivw",
  "mr_egger_regression", 
  'mr_weighted_mode',
  'mr_weighted_median'))
write.csv(singlesnp_results,'Hill-noLD-Lee35snpssinglesnp.csv')

#remove the snps that weren't used in the MR
dim(dat3)
singlesnp_results.3=singlesnp_results[1:35,]
dim(singlesnp_results.3)
dat3.1=dat3[which(dat3$SNP %in% singlesnp_results.3$SNP),]
dim(dat3.1)
dim(dat3)
dat3=dat3.1
dat3$samplesize.exposure=293723
write.csv(dat3,'Hill-noLD-Lee35snps_dat.csv')

#Generate r2 and Fstat
out <- directionality_test(dat3)
r2=out$snp_r2.exposure
r2 #0.005007605
n= 293723#smaller of the two
k= 35
p=dat3$eaf.exposure
Fstat <- r2*(n-1-k)/((1-r2)*k)
Fstat #42.23058

#Combine results
sin=singlesnp_results
res=mr_results
het<-mr_heterogeneity(dat3)
plt<-mr_pleiotropy_test(dat3)

all_res<-combine_all_mrresults(res,het,plt,sin,ao_slc=T,
  Exp=T,split.exposure=T,split.outcome=T)
all_res$r2=r2
all_res$lower_intercept=all_res$intercept-1.96*all_res$intercept_se
all_res$upper_intercept=all_res$intercept+1.96*all_res$intercept_se
all_res$Fstat=Fstat
all_res$I2=31.75 #from stata 

head(all_res[,c("Method","outcome","exposure","nsnp","b","se","or", "or_lci95",  
  "or_uci95","pval","intercept","intercept_se","intercept_pval",
  "Q","Q_df","Q_pval","consortium","ncase","ncontrol","pmid","population")])
write.csv(all_res,"combined_results_Hill_IQ_Lee_EduYears.csv")

#########################################################################
#### Univariable MR Lee EduYears on UK Biobank intelligence
#### Bidirectional analysis of education and intelligence: part 2
#########################################################################
setwd("C:/Users/charl/Dropbox/COH/IQ_school_sz/Edu_Years_GWAS/edu_years_IQ")
library(MendelianRandomization)
library(RadialMR)
library(TwoSampleMR)
ao <- available_outcomes()
exposure_dat <- read_exposure_data(
  filename = 'edu_years_GWAS.csv',
  sep = ',',
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
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('UKB-a:196'), 
  proxies = 1, rsq = 0.8, align_alleles = 1, 
  palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results <- mr(dat)

#Outlier's check
mrinput <- mr_input(bx = dat$beta.exposure, 
  bxse = dat$se.exposure, 
  by = dat$beta.outcome, byse=dat$se.outcome)
mregger <- mr_egger(mrinput, penalized=F) 
raddat <- format_radial(BXG=dat$beta.exposure, 
  BYG=dat$beta.outcome, 
  seBXG=dat$se.exposure, seBYG=dat$se.outcome, RSID=dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=1)
ivwrad2 <- ivw_radial(raddat, alpha=0.05, weights=1)
dim(ivwrad$outliers)[1] 
outliers=ivwrad$outliers[1]
outliers
myvars=outliers$SNP
dat2 <- dat[ ! dat$SNP %in% myvars, ]
dim(dat2)
dim(dat)
write.csv(dat2, "EduYears_outliers_gone_316_on_fluidIQ_dat.csv")
mr_results <- mr(dat2)
write.csv(mr_results,'Lee-outliers-gone_results.csv')
singlesnp_results=mr_singlesnp(dat2, 
  parameters = default_parameters(),
  single_method = "mr_wald_ratio", 
  all_method = c("mr_ivw",
  "mr_egger_regression", 
  'mr_weighted_mode',
  'mr_weighted_median'))
write.csv(singlesnp_results,'Lee-on-fluid-iq-outliers-305-gone_singlesnp.csv')

#remove the snps that weren't used in the MR
singlesnp_results.2=singlesnp_results[1:305,]
dim(singlesnp_results.2)
dat2.1=dat2[which(dat2$SNP %in% singlesnp_results.2$SNP),]
dim(dat2.1)
dim(dat2)
dat2=dat2.1
write.csv(dat2, 'Lee_EduYear_305snps-dat.csv')

#Read in the data for the Hill and Lee instruments that aren't in LD
not_ld=read.csv('unique-hill-lee-346.csv')
dim(not_ld)
newdata <- not_ld[ which(not_ld$study=='lee'), ]
dim(newdata)
dat3=dat2[which(dat2$SNP %in% newdata$SNP),]
dim(dat3)
mr_results <- mr(dat3)
write.csv(mr_results,'Lee-noLD_Hill-299.csv')

singlesnp_results=mr_singlesnp(dat3, 
  parameters = default_parameters(),
  single_method = "mr_wald_ratio", 
  all_method = c("mr_ivw",
  "mr_egger_regression", 
  'mr_weighted_mode',
  'mr_weighted_median'))

write.csv(singlesnp_results,'Lee-noLD_Hill-299singlesnp.csv')

#remove the snps that weren't used in the MR
dim(dat3)
singlesnp_results.3=singlesnp_results[1:299,]
dim(singlesnp_results.3)
dat3.1=dat3[which(dat3$SNP %in% singlesnp_results.3$SNP),]
dim(dat3.1)
dim(dat3)
dat3$samplesize.exposure=1131881
dat3=dat3.1
write.csv(dat3,'Lee-noLD-Hill-299-snpdat.csv')

#Generate r2 and Fstat
out <- directionality_test(dat3)
r2=out$snp_r2.exposure
r2 #0.01559162
n=108818 #smaller of the two
k= 299
p=dat3$eaf.exposure
Fstat <- r2*(n-1-k)/((1-r2)*k)
Fstat # 5.748394

#Combine results
sin=singlesnp_results
res=mr_results
het<-mr_heterogeneity(dat3)
plt<-mr_pleiotropy_test(dat3)
all_res<-combine_all_mrresults(res,het,plt,sin,ao_slc=T,
  Exp=T,split.exposure=T,split.outcome=T)
all_res$r2=r2
all_res$lower_intercept=all_res$intercept-1.96*all_res$intercept_se
all_res$upper_intercept=all_res$intercept+1.96*all_res$intercept_se
all_res$Fstat=Fstat
all_res$I2= 78.82#from stata 

head(all_res[,c("Method","outcome","exposure","nsnp","b","se","or", "or_lci95",  
  "or_uci95","pval","intercept","intercept_se","intercept_pval",
  "Q","Q_df","Q_pval","consortium","ncase","ncontrol","pmid","population")])
write.csv(all_res,"combined_results_Lee_EduYears_IQ.csv")

#########################################################################
#### Main Forest Plot for Univariable and Multivariable results
#########################################################################
setwd("C:/Users/charl/Dropbox/COH/IQ_school_sz/Edu_Years_GWAS/edu_years_IQ/forest")
library(metafor)
library(TwoSampleMR)
library(ggplot2)
library(MendelianRandomization)
library(RadialMR)
ao <- available_outcomes()
library(MRInstruments)
d=read.csv("forest_2.csv", fileEncoding="UTF-8-BOM")

#Assign values for plotting
#labs <- d$alloc
yi   <- d$b
sei  <- d$se
model=d$model

#decrease margins so the full space is used
#par(mfrow = c(1, 1))
#par(mar=c(4,4,1,2))

#meta-analysis of the log risk ratios using a random-effects model
res3 <- rma(yi=yi, sei=sei, method="FE", ai=d$n, bi=d$p, data=d,
            slab=paste(model)
)

forest(res3,
  xlim=c(-21,15),        ### adjust horizontal plot region limits #-2.5,4
  #subset=order(yi),     ### order by size of yi
  #slab=slab=paste(type, model, source, sep=","), annotate=TRUE, ### remove study labels and annotations
  efac=0,                ### remove vertical bars at end of CIs
  pch=15,                ### changing point symbol to filled circle
  col="blue",            ### change color of points/CIs
  psize=.75,             ### increase point size
  cex.lab=1, cex.axis=1, ### increase size of x-axis title/labels
  lty=c("solid","blank"),### remove horizontal line at top of plot
  refline=1, xlab="Odds ratio (95% CI)",
  #mlab="Summary Estimate", 
  transf=exp, 
  ilab=cbind(d$n),        #,d$p
  ilab.xpos=c(5), cex=1.5)#,8

#switch to bold font
par(font=2)
text(c(5), 18, c("SNPs"), cex=1.5) #,8,"P"
text(-21, 18, "Mendelian Randomization Model", pos=4, cex=1.5)
text(15,18, "OR [95% CI]", pos=2, cex=1.5)
#########################################################################
