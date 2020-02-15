#Analysis conducted by Daniel Jones. R version 3.6.2, MR TwoSample Package 0.5.0
#Based on script originally written by Robyn Wootton
#April 2019
######################################################################################################################
#Contents
#1. Load packages

###Lifetime smoking is the exposure
#2. Read in smoking exposure data
#3. Read in IBD outcome data
#4. Harmonise data
#5. Run 2 sample MR
#6. Steiger filtering
#7. Plot results
#8. Regression dilution I2 GX
#9. Simex correction 

######################################################################################################################################
#1. Load packages
######################################################################################################################################
rm(list=ls(all=TRUE))

install.packages("devtools")
library(devtools)
install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
install.packages("ggplot2")
library(ggplot2)
install.packages("knitr")
library(knitr)
install_github('qingyuanzhao/mr.raps')
library(mr.raps)
install.packages("Cairo")
library(Cairo)
install.packages("curl")
library(curl)

######################################################################################################################################
#2. Read in lifetime smoking data or smoking initiation data
######################################################################################################################################

#Lifetime smoking 
lsi_exp_dat <- read_exposure_data(
       filename = "CSI&IBD_Proxies.txt",
       sep = ",",
       snp_col = "SNP",
       beta_col = "beta.exposure",
       se_col = "se.exposure",
       effect_allele_col = "effect_allele.exposure",
       other_allele_col = "other_allele.exposure",
       eaf_col = "eaf.exposure",
       pval_col = "pval.exposure",
   )
lsi_exp_dat$exposure<-"Lifetime_Smoking"

#Standardise the betas and SE
lsi_exp_dat$beta.exposure<-lsi_exp_dat$beta.exposure/0.6940093
lsi_exp_dat$se.exposure<-lsi_exp_dat$se.exposure/0.6940093

#Smoking Initiation 
si_exp_dat <- read_exposure_data(
        filename = "SInit&IBD_Proxies.txt",
        sep = ",",
        snp_col = "rsID",
        beta_col = "Beta",
        se_col = "SE",
        effect_allele_col = "Alternate.Allele",
        other_allele_col = "Reference.Allele",
        eaf_col = "Alternate.Allele.Frequency",
        pval_col = "Pvalue",
        samplesize_col = "N"
  )
si_exp_dat $exposure<-"Smoking_Initiation"

######################################################################################################################################
#3. Read in relevant outcome data
######################################################################################################################################


#UC
uc_outcome_dat <- read_outcome_data(
         filename = "uc_build37_45975_20161107_rsID.txt",
         sep = ",", 
         snp_col = "rsids",
         beta_col = "Effect",
         se_col = "StdErr",
         effect_allele_col = "Allele2",
         other_allele_col = "Allele1",
         pval_col = "P.value")
uc_outcome_dat $outcome<-"UC"

#Crohn's
crohns_outcome_dat <- read_outcome_data(
         filename = "cd_build37_40266_20161107_rsIDs.txt",
         sep = ",",
         snp_col = "SNP",
         beta_col = "Effect",
         se_col = "StdErr",
         effect_allele_col = "Allele2",
         other_allele_col = "Allele1",
         pval_col = "P.value"
              )
crohns_outcome_dat $outcome<-"Crohn's"

#IBD
ibd_outcome_dat <- read_outcome_data(
         filename = "ibd_build37_59957_20161107_rsIDs.txt",
         sep = ",",
         snp_col = "rsids",
         beta_col = "Effect",
         se_col = "StdErr",
         effect_allele_col = "Allele2",
         other_allele_col = "Allele1",
         pval_col = "P.value"
     )
ibd_outcome_dat $outcome<-"IBD"

######################################################################################################################################
#4. Harmonise data
######################################################################################################################################

#LSI-UC
dat1 <- harmonise_data( 
  	exposure_dat = lsi_exp_dat,
  	outcome_dat = uc_outcome_dat,
  	action = 2
  )


#LSI-Crohns
dat2 <- harmonise_data( 
  	exposure_dat = lsi_exp_dat,
  	outcome_dat = crohns_outcome_dat,
  	action = 2
  )


#LSI-IBD
dat3 <- harmonise_data( 
  	exposure_dat = lsi_exp_dat,
  	outcome_dat = ibd_outcome_dat,
  	action = 2
  )


#SI-UC
dat4 <- harmonise_data( 
  	exposure_dat = si_exp_dat,
  	outcome_dat = uc_outcome_dat,
  	action = 2
  )


#SI-Crohns
dat5 <- harmonise_data( 
  	exposure_dat = si_exp_dat,
  	outcome_dat = crohns_outcome_dat,
  	action = 2
  )


#SI-IBD
dat6 <- harmonise_data( 
  	exposure_dat = si_exp_dat,
  	outcome_dat = ibd_outcome_dat,
  	action = 2
  )


######################################################################################################################################
#5. Run MR
######################################################################################################################################
#LSI-UC
mr_het1 <- mr_heterogeneity(dat1)
str(mr_het1)
mr_ruck1 <- mr_rucker(dat1)
res1 <- mr(dat1, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int1 <- mr_pleiotropy_test(dat1)
res1

or1<-generate_odds_ratios(res1)
or1


#LSI-Crohns
mr_het2 <- mr_heterogeneity(dat2)
str(mr_het2)
mr_ruck2 <- mr_rucker(dat2)
res2 <- mr(dat2, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int2 <- mr_pleiotropy_test(dat2)
res2

or2<-generate_odds_ratios(res2)
or2


#LSI-IBD
mr_het3 <- mr_heterogeneity(dat3)
str(mr_het3)
mr_ruck3 <- mr_rucker(dat3)
res3 <- mr(dat3, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int3 <- mr_pleiotropy_test(dat3)
res3

or3<-generate_odds_ratios(res3)
or3



#SI-UC
mr_het4 <- mr_heterogeneity(dat4)
str(mr_het4)
mr_ruck4 <- mr_rucker(dat4)
res4 <- mr(dat4, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int4 <- mr_pleiotropy_test(dat4)
res4

or4<-generate_odds_ratios(res4)
or4



#SI-Crohns
mr_het5 <- mr_heterogeneity(dat5)
str(mr_het5)
mr_ruck5 <- mr_rucker(dat5)
res5 <- mr(dat5, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int5 <- mr_pleiotropy_test(dat5)
res5

or5<-generate_odds_ratios(res5)
or5



#SI-IBD
mr_het6 <- mr_heterogeneity(dat6)
str(mr_het6)
mr_ruck6 <- mr_rucker(dat6)
res6<- mr(dat6, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int6 <- mr_pleiotropy_test(dat6)
res6

or6<-generate_odds_ratios(res6)
or6

#merge results together and save
results <- rbind (or1, or2, or3, or4, or5, or6)
write.csv(results, "Two-sample MR results - smoking to IBD.csv", quote=F, row.names=F)

######################################################################################################################################
#6. Steiger Filtering
######################################################################################################################################

#Ensure relevant data harmonised

#For lifetime smoking and UC
dat1$samplesize.exposure<-463003 
dat1$units.exposure<-"SD"
dat1$units.outcome <- "log odds"
dat1$samplesize.outcome<-59957 
dat1$ncase.outcome<-25042
dat1$ncontrol.outcome<-34915
dat1$prevalence.outcome<-0.01 
dat1$eaf.outcome <- dat1$eaf.exposure
steiger1<- steiger_filtering(dat1)
table(steiger1$steiger_dir)

FALSE  TRUE 
   13   112 

#For smoking initiation and UC
dat4$samplesize.exposure<-1232091 
dat4$units.exposure<-"SD"
dat4$units.outcome <- "log odds"
dat4$samplesize.outcome<-59957 
dat4$ncase.outcome<-25042
dat4$ncontrol.outcome<-34915
dat4$prevalence.outcome<-0.01 
dat4$eaf.outcome <- dat4$eaf.exposure
steiger4 <- steiger_filtering(dat4)
table(steiger4$steiger_dir)

FALSE  TRUE 
   20   354 

#For smoking initiation and Crohn's disease
dat5$samplesize.exposure<-1232091 
dat5$units.exposure<-"SD"
dat5$units.outcome <- "log odds"
dat5$samplesize.outcome<-59957
dat5$ncase.outcome<-25042
dat5$ncontrol.outcome<-34915
dat5$prevalence.outcome<-0.01 
dat5$eaf.outcome <- dat5$eaf.exposure
steiger5 <- steiger_filtering(dat5)
table(steiger5$steiger_dir)

FALSE  TRUE 
   21   352 

#For lifetime smoking and Crohn's
dat2$samplesize.exposure<-463003 
dat2$units.exposure<-"SD"
dat2$units.outcome <- "log odds"
dat2$samplesize.outcome<-59957
dat2$ncase.outcome<-25042
dat2$ncontrol.outcome<-34915
dat2$prevalence.outcome<-0.01 
dat2$eaf.outcome <- dat2$eaf.exposure
steiger2 <- steiger_filtering(dat2)
table(steiger2$steiger_dir)

FALSE  TRUE 
    17   108 

#For lifetime smoking and IBD
dat3$samplesize.exposure<-463003 
dat3$units.exposure<-"SD"
dat3$units.outcome <- "log odds"
dat3$samplesize.outcome<-59957
dat3$ncase.outcome<-25042
dat3$ncontrol.outcome<-34915
dat3$prevalence.outcome<-0.01 
dat3$eaf.outcome <- dat3$eaf.exposure
steiger3 <- steiger_filtering(dat3)
table(steiger3$steiger_dir)

FALSE  TRUE 
    10   115 

#For smoking initiation and IBD
dat6$samplesize.exposure<-1232091 
dat6$units.exposure<-"SD"
dat6$units.outcome <- "log odds"
dat6$samplesize.outcome<-59957
dat6$ncase.outcome<-25042
dat6$ncontrol.outcome<-34915
dat6$prevalence.outcome<-0.01
dat6$eaf.outcome <- dat6$eaf.exposure
steiger6 <- steiger_filtering(dat6)
table(steiger6$steiger_dir)

FALSE  TRUE 
   17   357 

##To run MR analysis without false SNPs
steiger1<-subset(steiger1, steiger_dir=="TRUE", select=SNP:mr_keep)
steiger2<-subset(steiger2, steiger_dir=="TRUE", select=SNP:mr_keep)
steiger3<-subset(steiger3, steiger_dir=="TRUE", select=SNP:mr_keep)
steiger4<-subset(steiger4, steiger_dir=="TRUE", select=SNP:mr_keep)
steiger5<-subset(steiger5, steiger_dir=="TRUE", select=SNP:mr_keep)
steiger6<-subset(steiger6, steiger_dir=="TRUE", select=SNP:mr_keep)

#Repeat section 5
#LSI-UC
mr_het7 <- mr_heterogeneity(steiger1)
str(mr_het7)
mr_ruck7 <- mr_rucker(steiger1)
res7 <- mr(steiger1, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int7 <- mr_pleiotropy_test(steiger1)
res7

or7<-generate_odds_ratios(res7)
or7

#LSI-Crohns
mr_het8 <- mr_heterogeneity(steiger2)
str(mr_het8)
mr_ruck8 <- mr_rucker(steiger2)
res8 <- mr(steiger2, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int8 <- mr_pleiotropy_test(steiger2)
res8

or8<-generate_odds_ratios(res8)
or8

#LSI-IBD
mr_het9 <- mr_heterogeneity(steiger3)
str(mr_het9)
mr_ruck9 <- mr_rucker(steiger3)
res9 <- mr(steiger3, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int9 <- mr_pleiotropy_test(steiger3)
res9

or9<-generate_odds_ratios(res9)
or9

#SI-UC
mr_het10 <- mr_heterogeneity(steiger4)
str(mr_het10)
mr_ruck10 <- mr_rucker(steiger4)
res10 <- mr(steiger4, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int10 <- mr_pleiotropy_test(steiger4)
res10

or10<-generate_odds_ratios(res10)
or10

#SI-Crohns
mr_het11 <- mr_heterogeneity(steiger5)
str(mr_het11)
mr_ruck11 <- mr_rucker(steiger5)
res11 <- mr(steiger5, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int11 <- mr_pleiotropy_test(steiger5)
res11

or11<-generate_odds_ratios(res11)
or11

#SI-IBD
mr_het12 <- mr_heterogeneity(steiger6)
str(mr_het12)
mr_ruck12 <- mr_rucker(steiger6)
res12<- mr(steiger6, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
mr_egger_int12 <- mr_pleiotropy_test(steiger6)
res12

or12<-generate_odds_ratios(res12)
or12

#merge results together and save
results <- rbind (or7, or8, or9, or10, or11, or12)
write.csv(results, "Two-sample MR results - smoking to IBD(Steiger).csv", quote=F, row.names=F)

######################################################################################################################################
#7. Plot results
######################################################################################################################################
#For lifetime smoking and UC
res1_single <- mr_singlesnp(dat1, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res1_loo <- mr_leaveoneout(dat1)
p1 <- mr_scatter_plot(res1, dat1)
ggsave(p1[[1]], file="Lifetime_Smoking&UC_scatter.png", width=7, height=7)
p2 <- mr_forest_plot(res1_single)
ggsave(p2[[1]], file="Lifetime_Smoking&UC_forest.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res1_loo)
ggsave(p3[[1]], file="Lifetime_Smoking&UC_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res1_single)
ggsave(p4[[1]], file="Lifetime_Smoking&UC_funnel.png", width=7, height=7)
mr_report(dat1, output_path="Pathname", output_type = "html",
    author = "Author", study = "Lifetime_Smoking&UC_Report")

#For smoking initiation and UC
res4_single <- mr_singlesnp(dat4, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res4_loo <- mr_leaveoneout(dat4)
p1 <- mr_scatter_plot(res4, dat4)
ggsave(p1[[1]], file="Smoking_Initiation&UC_scatter.png", width=7, height=7)
p2 <- mr_forest_plot(res4_single)
ggsave(p2[[1]], file="Smoking_Initiation&UC_forest.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res4_loo)
ggsave(p3[[1]], file="Smoking_Initiation&UC_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res4_single)
ggsave(p4[[1]], file="Smoking_Initiation&UC_funnel.png", width=7, height=7)
mr_report(dat4, output_path="Pathname", output_type = "html",
    author = "Author", study = "Smoking_Initiation&UC_Report")
    
#For smoking initiation and Crohn's disease
res5_single <- mr_singlesnp(dat5, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res5_loo <- mr_leaveoneout(dat5)
p1 <- mr_scatter_plot(res5, dat5)
ggsave(p1[[1]], file="Smoking_Initiation&Crohns_scatter.png", width=7, height=7)
p2 <- mr_forest_plot(res5_single)
ggsave(p2[[1]], file="Smoking_Initiation&Crohns_forest.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res5_loo)
ggsave(p3[[1]], file="Smoking_Initiation&Crohns_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res5_single)
ggsave(p4[[1]], file="Smoking_Initiation&Crohns_funnel.png", width=7, height=7)
mr_report(dat5, output_path="Pathname", output_type = "html",
    author = "Author", study = "Smoking_Initiation&Crohns_Report")

#For lifetime smoking and Crohn's
res2_single <- mr_singlesnp(dat2, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res2_loo <- mr_leaveoneout(dat2)
p1 <- mr_scatter_plot(res2, dat2)
ggsave(p1[[1]], file="Lifetime_Smoking&Crohns_scatter.png", width=7, height=7)
p2 <- mr_forest_plot(res2_single)
ggsave(p2[[1]], file="Lifetime_Smoking&Crohns_forest.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res2_loo)
ggsave(p3[[1]], file="Lifetime_Smoking&Crohns_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res2_single)
ggsave(p4[[1]], file="Lifetime_Smoking&Crohns_funnel.png", width=7, height=7)
mr_report(dat2, output_path="Pathname", output_type = "html",
    author = "Author", study = "Lifetime_Smoking&Crohns_Report")

#For lifetime smoking and IBD
res3_single <- mr_singlesnp(dat3, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res3_loo <- mr_leaveoneout(dat3)
p1 <- mr_scatter_plot(res3, dat3)
ggsave(p1[[1]], file="Lifetime_Smoking&IBD_scatter.png", width=7, height=7)
p2 <- mr_forest_plot(res3_single)
ggsave(p2[[1]], file="Lifetime_Smoking&IBD_forest.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res3_loo)
ggsave(p3[[1]], file="Lifetime_Smoking&IBD_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res3_single)
ggsave(p4[[1]], file="Lifetime_Smoking&IBD_funnel.png", width=7, height=7)
mr_report(dat3, output_path="Pathname", output_type = "html",
    author = "Author", study = "Lifetime_Smoking&IBD_Report")
    
#For smoking initiation and IBD
res6_single <- mr_singlesnp(dat6, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
res6_loo <- mr_leaveoneout(dat6)
p1 <- mr_scatter_plot(res6, dat6)
ggsave(p1[[1]], file="Smoking_Initiation&IBD_scatter.png", width=7, height=7)
p2 <- mr_forest_plot(res6_single)
ggsave(p2[[1]], file="Smoking_Initiation&IBD_forest.png", width=7, height=7)
p3 <- mr_leaveoneout_plot(res6_loo)
ggsave(p3[[1]], file="Smoking_Initiation&IBD_loo.png", width=7, height=7) 
p4 <- mr_funnel_plot(res6_single)
ggsave(p4[[1]], file="Smoking_Iniation&IBD_funnel.png", width=7, height=7)
mr_report(dat6, output_path="Pathname", output_type = "html",
    author = "Author", study = "Smoking_Iniation&IBD_Report")    
                
######################################################################################################################################
#8. Regression Dilution
######################################################################################################################################
#Isq >0.9 = ok to do MR Egger
#0.3<Isq>0.9 = perform SIMEX correction
#Isq<0.3 do not report MR Egger at all

#Create I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
  }

#calculate Isq wieghted and unweighted in a loop
I2<-c()

###For Smoking Iniation and UC

#Update ouctome data values and then harmonise data
uc_outcome_dat.proxies<-1
uc_outcome_dat.rsq<-0.8
uc_outcome_dat.align_alleles<-1
uc_outcome_dat.maf_threshold<-0.3
uc_outcome_dat.palindromes<-1
dat7 <- harmonise_data(si_exp_dat, uc_outcome_dat, action = 1)

#Rename required columns
for(i in 1:length(dat7)){
   	Vars<-dat7[i]
   }
dat7$BetaXG<-dat7$beta.exposure
dat7$seBetaXG<-dat7$se.exposure
dat7$seBetaYG<-dat7$se.outcome
BetaXG   = dat7$BetaXG
seBetaXG = dat7$seBetaXG 
seBetaYG = dat7$seBetaYG
BXG             = abs(BetaXG)

# Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))

#Save results
output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("Smoking_Initiation", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="Smoking_Initiation&UC_regression_dilution_isq_weighted.csv", row.names = FALSE)

###For Lifetime Smoking and UC
#Harmonise data
dat8 <- harmonise_data(lsi_exp_dat, uc_outcome_dat, action = 1)

#Rename columns
for(i in 1:length(dat8)){
   	Vars<-dat8[i]
   }
dat8$BetaXG<-dat8$beta.exposure
dat8$seBetaXG<-dat8$se.exposure
dat8$seBetaYG<-dat8$se.outcome
BetaXG   = dat8$BetaXG
seBetaXG = dat8$seBetaXG 
seBetaYG = dat8$seBetaYG
BXG             = abs(BetaXG)

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))
I2<-c()

#Save results
output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("Lifetime_Smoking", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="Lifetime_Smoking&UC_regression_dilution_isq_weighted.csv", row.names = FALSE)

###For Smoking Initiation and Crohn's disease
#Update outcome data values and then harmonise data
crohns_outcome_dat.proxies<-1
crohns_outcome_dat.rsq<-0.8
crohns_outcome_dat.align_alleles<-1
crohns_outcome_dat.maf_threshold<-0.3
crohns_outcome_dat.palindromes<-1
dat9 <- harmonise_data(si_exp_dat, crohns_outcome_dat, action = 1)

#Rename required columns
for(i in 1:length(dat9)){
   	Vars<-dat9[i]
   }
dat9$BetaXG<-dat9$beta.exposure
dat9$seBetaXG<-dat9$se.exposure
dat9$seBetaYG<-dat9$se.outcome
BetaXG   = dat9$BetaXG
seBetaXG = dat9$seBetaXG 
seBetaYG = dat9$seBetaYG
BXG             = abs(BetaXG)

# Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))

#Save results
output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-c()
I2<-rbind(I2, output)
colnames(I2) <- c("Smoking_Initiation", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="Smoking_Initiation&Crohns_regression_dilution_isq_weighted.csv", row.names = FALSE)

###For Lifetime Smoking and Crohn's disease
#Harmonise data
dat10 <- harmonise_data(lsi_exp_dat, crohns_outcome_dat, action = 1)

#Rename columns
for(i in 1:length(dat10)){
   	Vars<-dat10[i]
   }
dat10$BetaXG<-dat10$beta.exposure
dat10$seBetaXG<-dat10$se.exposure
dat10$seBetaYG<-dat10$se.outcome
BetaXG   = dat10$BetaXG
seBetaXG = dat10$seBetaXG 
seBetaYG = dat10$seBetaYG
BXG             = abs(BetaXG)

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))
I2<-c()

#Save results
output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("Lifetime_Smoking", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="Lifetime_Smoking&Crohns_regression_dilution_isq_weighted.csv", row.names = FALSE)

###For Smoking Initiation and IBD
#Update outcome data values and then harmonise data
ibd_outcome_dat.proxies<-1
ibd_outcome_dat.rsq<-0.8
ibd_outcome_dat.align_alleles<-1
ibd_outcome_dat.maf_threshold<-0.3
ibd_outcome_dat.palindromes<-1
dat11<- harmonise_data(si_exp_dat, ibd_outcome_dat, action = 1)

#Rename required columns
for(i in 1:length(dat11)){
   	Vars<-dat11[i]
   }
dat11$BetaXG<-dat11$beta.exposure
dat11$seBetaXG<-dat11$se.exposure
dat11$seBetaYG<-dat11$se.outcome
BetaXG   = dat11$BetaXG
seBetaXG = dat11$seBetaXG 
seBetaYG = dat11$seBetaYG
BXG             = abs(BetaXG)

# Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))

#Save results
output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("Smoking_Initiation", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="Smoking_Initiation&IBD_regression_dilution_isq_weighted.csv", row.names = FALSE)

###For Lifetime Smoking and IBD
#Harmonise data
dat12 <- harmonise_data(lsi_exp_dat, ibd_outcome_dat, action = 1)

#Rename columns
for(i in 1:length(dat12)){
   	Vars<-dat12[i]
   }
dat12$BetaXG<-dat12$beta.exposure
dat12$seBetaXG<-dat12$se.exposure
dat12$seBetaYG<-dat12$se.outcome
BetaXG   = dat12$BetaXG
seBetaXG = dat12$seBetaXG 
seBetaYG = dat12$seBetaYG
BXG             = abs(BetaXG)

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG)
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))
I2<-c()

#Save results
output<-cbind(Vars, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("Lifetime_Smoking", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="Lifetime_Smoking&IBD_regression_dilution_isq_weighted.csv", row.names = FALSE)

######################################################################################################################################
#9.Simex corrections
######################################################################################################################################

#Install Packages
install.packages("simex")
library(simex)

###For Smoking Initiation and UC
#create empty dataframe to store output
simexegger<-c()

#Rename required columns
for(i in 1:length(dat7)){
   	Vars<-dat7[i]
   }
dat7$BetaXG<-dat7$beta.exposure
dat7$BetaYG<-dat7$beta.outcome
dat7$seBetaXG<-dat7$se.exposure
dat7$seBetaYG<-dat7$se.outcome
BetaXG <- dat7$BetaXG
BetaYG <- dat7$BetaYG
seBetaXG <- dat7$seBetaXG
seBetaYG <- dat7$seBetaYG
BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG) 

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Run Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

#Extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

#Convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("Smoking_Initiation", "b", "se", "pval")
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

#Extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

#Save results
output<-cbind(Vars, beta1, lci1, uci1, p1, or1, lcior1, ucior1, beta2, lci2, uci2, p2, or2, lcior2, ucior2)
simexegger<-rbind(simexegger, output)
colnames(simexegger) <- c("Smoking_Initiation", "beta_weighted", "lowerCI_weighted", "upperCI_weighted", "p_weighted", "OR_weighted", "lCIOR_weighted", "uCIOR_weighted", "beta_unweighted", "lowerCI_unweighted", "upperCI_unweighted", "p_unweighted", "OR_unweighted", "lCIOR_unweighted", "uCIOR_unweighted")
write.csv(simexegger, file="Smoking_Initiation&UC_mreggersimex_weighted.csv", row.names = FALSE)

###For Lifetime Smoking and UC
#create empty dataframe to store output
simexegger<-c()

#Rename required columns
for(i in 1:length(dat8)){
   	Vars<-dat8[i]
   }
dat8$BetaYG<-dat8$beta.outcome
BetaXG <- dat8$BetaXG
BetaYG <- dat8$BetaYG
seBetaXG <- dat8$seBetaXG
seBetaYG <- dat8$seBetaYG
BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG) 

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Run Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

#Extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

#Convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("Lifetime_Smoking", "b", "se", "pval")
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

#Extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

#Save results
output<-cbind(Vars, beta1, lci1, uci1, p1, or1, lcior1, ucior1, beta2, lci2, uci2, p2, or2, lcior2, ucior2)
simexegger<-rbind(simexegger, output)
colnames(simexegger) <- c("Lifetime_Smoking", "beta_weighted", "lowerCI_weighted", "upperCI_weighted", "p_weighted", "OR_weighted", "lCIOR_weighted", "uCIOR_weighted", "beta_unweighted", "lowerCI_unweighted", "upperCI_unweighted", "p_unweighted", "OR_unweighted", "lCIOR_unweighted", "uCIOR_unweighted")
write.csv(simexegger, file="Lifetime_Smoking&UC_mreggersimex_weighted.csv", row.names = FALSE)

###For Smoking Initiation and Crohn's disease
#Create empty dataframe to store output
simexegger<-c()

#Rename required columns
for(i in 1:length(dat9)){
   	Vars<-dat9[i]
   }
dat9$BetaYG<-dat9$beta.outcome
BetaXG <- dat9$BetaXG
BetaYG <- dat9$BetaYG
seBetaXG <- dat9$seBetaXG
seBetaYG <- dat9$seBetaYG
BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG) 

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Run Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

#Extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

#Convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("Smoking_Initiation", "b", "se", "pval")
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

#Extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

#Save results
output<-cbind(Vars, beta1, lci1, uci1, p1, or1, lcior1, ucior1, beta2, lci2, uci2, p2, or2, lcior2, ucior2)
simexegger<-rbind(simexegger, output)
colnames(simexegger) <- c("Smoking_Initiation", "beta_weighted", "lowerCI_weighted", "upperCI_weighted", "p_weighted", "OR_weighted", "lCIOR_weighted", "uCIOR_weighted", "beta_unweighted", "lowerCI_unweighted", "upperCI_unweighted", "p_unweighted", "OR_unweighted", "lCIOR_unweighted", "uCIOR_unweighted")
write.csv(simexegger, file="Smoking_Iniation&Crohns_mreggersimex_weighted.csv", row.names = FALSE)

###For Lifetime Smoking and Crohn's disease
#Create empty dataframe to store output
simexegger<-c()

#Rename required columns
for(i in 1:length(dat10)){
   	Vars<-dat10[i]
   }
dat10$BetaYG<-dat10$beta.outcome
BetaXG <- dat10$BetaXG
BetaYG <- dat10$BetaYG
seBetaXG <- dat10$seBetaXG
seBetaYG <- dat10$seBetaYG
BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG) 

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Run Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

#Extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

#Convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("Lifetime_Smoking", "b", "se", "pval")
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

#Extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

#Save results
output<-cbind(Vars, beta1, lci1, uci1, p1, or1, lcior1, ucior1, beta2, lci2, uci2, p2, or2, lcior2, ucior2)
simexegger<-rbind(simexegger, output)
colnames(simexegger) <- c("Lifetime_Smoking", "beta_weighted", "lowerCI_weighted", "upperCI_weighted", "p_weighted", "OR_weighted", "lCIOR_weighted", "uCIOR_weighted", "beta_unweighted", "lowerCI_unweighted", "upperCI_unweighted", "p_unweighted", "OR_unweighted", "lCIOR_unweighted", "uCIOR_unweighted")
write.csv(simexegger, file="Lifetime_Smoking&Crohns_mreggersimex_weighted.csv", row.names = FALSE)

###For Smoking Initiation and IBD
#Create empty dataframe to store output
simexegger<-c()

#Rename required columns
for(i in 1:length(dat11)){
   	Vars<-dat11[i]
   }
dat11$BetaYG<-dat11$beta.outcome
BetaXG <- dat11$BetaXG
BetaYG <- dat11$BetaYG
seBetaXG <- dat11$seBetaXG
seBetaYG <- dat11$seBetaYG
BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG) 

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Run Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

#Extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

#Convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("Smoking_Initiation", "b", "se", "pval")
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

#Extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

#Save results
output<-cbind(Vars, beta1, lci1, uci1, p1, or1, lcior1, ucior1, beta2, lci2, uci2, p2, or2, lcior2, ucior2)
simexegger<-rbind(simexegger, output)
colnames(simexegger) <- c("Smoking_Initiation", "beta_weighted", "lowerCI_weighted", "upperCI_weighted", "p_weighted", "OR_weighted", "lCIOR_weighted", "uCIOR_weighted", "beta_unweighted", "lowerCI_unweighted", "upperCI_unweighted", "p_unweighted", "OR_unweighted", "lCIOR_unweighted", "uCIOR_unweighted")
write.csv(simexegger, file="Smoking_Initiation&IBD_mreggersimex_weighted.csv", row.names = FALSE)

###For Lifetime Smoking and IBD
#Create empty dataframe to store output
simexegger<-c()

#Rename required columns
for(i in 1:length(dat12)){
   	Vars<-dat12[i]
   }
dat12$BetaYG<-dat12$beta.outcome
BetaXG <- dat12$BetaXG
BetaYG <- dat12$BetaYG
seBetaXG <- dat12$seBetaXG
seBetaYG <- dat12$seBetaYG
BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG) 

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Run Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

#Extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

#Convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("Lifetime_Smoking", "b", "se", "pval")
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

#Extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

#Save results
output<-cbind(Vars, beta1, lci1, uci1, p1, or1, lcior1, ucior1, beta2, lci2, uci2, p2, or2, lcior2, ucior2)
simexegger<-rbind(simexegger, output)
colnames(simexegger) <- c("Lifetime_Smoking", "beta_weighted", "lowerCI_weighted", "upperCI_weighted", "p_weighted", "OR_weighted", "lCIOR_weighted", "uCIOR_weighted", "beta_unweighted", "lowerCI_unweighted", "upperCI_unweighted", "p_unweighted", "OR_unweighted", "lCIOR_unweighted", "uCIOR_unweighted")
write.csv(simexegger, file="Lifetime_Smoking&IBD_mreggersimex_weighted.csv", row.names = FALSE)