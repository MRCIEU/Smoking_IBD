#Individual level MR analysis in UK Biobank to explore the effects of CHRNA5 on IBD risk
#Robyn Wootton August 2019
#######################################################################################################################
#Contents
#1. Read in the IBD data for UK Biobank
#2. Tidy up data - name variables appropriately
#3. Exclude participants who have withdrawn consent
#4. Read in SNP and other smoking variables
#5. Split into ever and never smokers
#6. Get the demographics per genotype of the SNP
#7. Test Assosiation of genetic score with IBD in smokers
#8. Test Assosiation of genetic score with IBD in non-smokers 
#9. Plot the results

##########################################################################
#1. Read in data
##########################################################################
rm(list=ls(all=TRUE)) #empties your R environment

read <- read.dta13("20170613_sleep_icd10.dta")
str(read)

#subset the IBD diagnosis vars using ICD10 codes
#k50 = Chrons
#k51 = UC

#Identify ICD 10 codes using UK Biobank Variable Number [main diagnoses = 41202]
diagnosis_main<-read[ , grepl( "41202" , names( read ) ) ]
str(diagnosis_main)

#[secondary diagnoses = 41204]
diagnosis_second<-read[ , grepl( "41204" , names( read ) ) ]
str(diagnosis_second)

#stick all of the datasets together
data<-data.frame(read$n_eid, diagnosis_main, diagnosis_second)
str(data)

##########################################################################
#2. Tidy up data - name variables appropriately
##########################################################################
#Set all blank values as missing
data[c(2:816)][data[c(2:816)] == ""] <-NA
str(data)

#restrict all of the ICD 10 codes to the first two numbers only
for (i in 2:816){
	data[, i] <- substr(data[, i], 0, 3)
}

#save dataset
write.csv(data, "ICD10_codes.csv", quote=F, row.names=F)


#make every code NA if it is not UC or chrons
test<-data

for (i in 2:816){
	test[, i] <- ifelse(test[, i]!="K50" & test[, i]!="K51", 0, test[, i])
}
table(test$s_41202_0_0) #only K50 and K51 should remain  
table(test$s_41202_0_10) 

#Now do this for IBD total, just chrons and just UC
IBD <- data
for (i in 2:816){
	IBD[, i] <- ifelse(IBD[, i]!="K50" & IBD[, i]!="K51", 0, 1)
}
table(IBD $s_41202_0_10)  
IBD$ibd_sum <- rowSums(IBD[c(2:816)], na.rm=T)
table(IBD$ibd_sum)
#make all values above 1 into 1
IBD$ibd[IBD$ibd_sum>1] <-1
table(IBD$ibd)

#Now do the same for Chrons only
#k50 = Chrons
Chrons <- data
for (i in 2:816){
	Chrons[, i] <- ifelse(Chrons[, i]!="K50", 0, 1)
}
table(Chrons $s_41202_0_10)
Chrons $chrons_sum <- rowSums(Chrons[c(2:816)], na.rm=T)
table(Chrons $chrons_sum)
#make all values above 1 into 1
Chrons $chrons[Chrons $chrons_sum>1] <-1
table(Chrons $chrons)

#k51 = UC
UC <- data
for (i in 2:816){
	UC[, i] <- ifelse(UC[, i]!="K51", 0, 1)
}
table(UC $s_41202_0_10)
UC $uc_sum <- rowSums(UC[c(2:816)], na.rm=T)
table(UC $uc_sum)
#make all values above 1 into 1
UC $uc[UC $uc_sum>1] <-1
table(UC $uc)

#subset just the relevant columns
IBD2 <- subset(IBD, select=c("read.n_eid", "ibd"))
str(IBD2)
Chrons2 <- subset(Chrons, select=c("read.n_eid", "chrons"))
str(Chrons2)
UC2 <- subset(UC, select=c("read.n_eid", "uc"))
str(UC2)

#merge these three columns together and save
fin <- merge(IBD2, Chrons2, by="read.n_eid", all=T)
str(fin) #N=502621
final <- merge(fin, UC2, by="read.n_eid", all=T)
str(final) #N=502621

#save dataset
write.csv(final, "UKBB_IBDvariables.csv", quote=F, row.names=F)

#########################################################################
#3. Exclude participants who have withdrawn consent
#########################################################################
#Read in the exclusion list
ex<-read.csv("w9142_20181016.csv", header=F)
str(ex)

ex$exclude<-1

#Merge these into the data file
df<-merge(data, ex, by.x="read.n_eid", by.y="V1", all.x=T)
str(df) #N= 502621
dat<-subset(df, is.na(df$exclude))
str(dat) #N=502,543

######################################################################################
#4. Read in SNP and other smoking variables
######################################################################################
#Read in smoking variables (previously subsetted)
data<-read.csv("Feb2018_smoking_mortality.csv", header=T)
str(data)

#merge with the IBD data
final <- merge(dat, data, by.x="read.n_eid", by.y="ID", all.x=T)
head(final)

# Apply the genetic exclusions
table(final $snp)

#Exclude people who are not european ancestry using IEU guide (DOI: 10.5523/bris.3074krb6t2frj29yh2b03x3wxj)
#Read in the genetic data exclusions for BMI - here I have the 500k with exclusions applied and I can subset on this
scores<-read.csv("UKBiobank_BMI_PRS.csv", header=T) 
str(scores) #Check N=337,115

#Read in the ID linker file
library(readstata13)
mergelist <- read.dta13("matching_id.dta")
str(mergelist)

data2 <-merge(final, mergelist, by.x="read.n_eid", by.y="appid9142", all.x=T)
str(data2)

df<-merge(data2, scores, by.x="appid8786", by.y="IID", all.x=T)
str(df) #N=502,543

#Now subset only those who have a score for the 500k
data3<-subset(df, complete.cases(df$SCORE500))
str(data3) #N = 337, 053

######################################################################################
#5. Split into ever and never smokers
######################################################################################
ever<-subset(data3, data3 $nevev=="Ever")
str(ever) #N= 151, 822

never<-subset(data3, data3 $nevev =="Never")
str(never) #N= 184044

current<-subset(data3, data3 $smoking_status =="Current")
str(current) #N= 33357

former <- subset(data3, data3 $smoking_status =="Previous")
str(former) #N= 118465

table(data3$sex)
table(data3$ageatstart)
mean(data3$ageatstart, na.rm=T)
sd(data3$ageatstart, na.rm=T)

######################################################################################
#6. Get the demographics per genotype of the SNP
######################################################################################

#Split the whole dataset into those with 0, 1 and 2 increasing alleles
snp0<-subset(ever, ever$snp ==0)
str(snp0)
snp0.current<-subset(snp0, snp0$smoking_status=="Current")
snp0.former<-subset(snp0, snp0$smoking_status=="Previous")
mean(snp0.current$cpd, na.rm=T)
summary(snp0.current$cpd)
mean(snp0$ageatstart, na.rm=T)
str(snp0.current)
str(snp0.former)

snp1<-subset(ever, ever $snp ==1)
str(snp1)
snp1.current<-subset(snp1, snp1$smoking_status=="Current")
mean(snp1.current$cpd, na.rm=T)
summary(snp1.current$cpd)
mean(snp1$ageatstart, na.rm=T)
str(snp1.current)

snp2<-subset(ever, ever$snp ==2)
str(snp2)
snp2.current<-subset(snp2, snp2$smoking_status=="Current")
str(snp2.current)
mean(snp2.current$cpd, na.rm=T) #check it does indeed increase per allele by about 1 cig per day
summary(snp2.current$cpd)
mean(snp2$ageatstart, na.rm=T)

#p-value for cpd
ever.current<-subset(ever, ever$smoking_status=="Current")
m1<-aov(ever.current $cpd~ ever.current$snp)
summary(m1)

#Use Chisquare to check the frequencies by snp
library(MASS)       # load the MASS package 

#snp by smoking status
tbl<-table(ever$snp, ever$smoking_status)  
prop.table(tbl, 1) #get proportions
tbl<-table(ever $snp, ever $smoking_status, exclude = NULL) #get NA
#p-value for smoking status
m1<-aov(ever$snp ~ ever$smoking_status)
summary(m1)

#snp by sex
tbl<-table(ever$snp, ever$sex.x) 
chisq.test(tbl) 
prop.table(tbl, 1) #get proportions

#snp by education
tbl<-table(ever$snp, ever$edu) 
chisq.test(tbl) #0=primary, 1=secondary, 2=tertiary
prop.table(tbl, 1) #get proportions
table(ever$snp, ever$edu, exclude = NULL) #get NA

#snp by alcohol (where 0=never and 4=daily)
tbl<-table(ever$snp, ever$alcohol) 
chisq.test(tbl) 
prop.table(tbl, 1) #get proportions
table(ever$snp, ever$alcohol, exclude = NULL) #get NA

#snp by IBD status 
tbl<-table(ever $snp, ever $ibd) 
chisq.test(tbl) 
prop.table(tbl, 1) #get proportions
table(ever $snp, ever $ibd, exclude = NULL) #get NA

#snp by UC
tbl<-table(ever $snp, ever $uc) 
chisq.test(tbl) 
prop.table(tbl, 1) #get proportions
table(ever $snp, ever $uc, exclude = NULL) #get NA

#snp by chrons
tbl<-table(ever $snp, ever $chrons) 
chisq.test(tbl) 
prop.table(tbl, 1) #get proportions
table(ever $snp, ever $chrons, exclude = NULL) #get NA


#anova of age and snp to get p value
m1<-aov(ever$ageatstart~ever$snp)
summary(m1)

###Never smokers
#Split the whole dataset into those with 0, 1 and 2 increasing alleles
snp0<-subset(never, never $snp==0)
str(snp0)
mean(snp0$ageatstart)

snp1<-subset(never, never$snp==1)
str(snp1)
mean(snp1$ageatstart)

snp2<-subset(never, never $snp==2)
str(snp2)
mean(snp2$ageatstart)

#snp by sex
tbl<-table(never$snp, never$sex.x) 
chisq.test(tbl) 
prop.table(tbl, 1) #get proportions

#snp by education
tbl<-table(never $snp, never $edu) 
chisq.test(tbl) #0=primary, 1=secondary, 2=tertiary
prop.table(tbl, 1) #get proportions
table(never$snp, never$edu, exclude = NULL) #get NA

#snp by alcohol (where 0=never and 4=daily)
tbl<-table(never$snp, never$alcohol) 
chisq.test(tbl) 
prop.table(tbl, 1) #get proportions
table(never$snp, never$alcohol, exclude = NULL) #get NA


#anova of age and cpd to get p value
m1<-aov(never$ageatstart~ never$snp)
summary(m1)

#snp by IBD status 
tbl<-table(never$snp, never$ibd) 
chisq.test(tbl) 
prop.table(tbl, 1) #get proportions
table(never$snp, never$ibd, exclude = NULL) #get NA

#snp by UC
tbl<-table(never$snp, never$uc) 
chisq.test(tbl) 
prop.table(tbl, 1) #get proportions
table(never$snp, never$uc, exclude = NULL) #get NA

#snp by chrons
tbl<-table(never$snp, never$chrons) 
chisq.test(tbl) 
prop.table(tbl, 1) #get proportions
table(never$snp, never$chrons, exclude = NULL) #get NA

######################################################################################
#7. Test Assosiation of SNP with IBD in smokers
######################################################################################
###ibd
#ever smokers
model1<- glm(ibd ~ snp + ageatstart + sex, data=ever, family=binomial(link="logit"))
mod1<-summary(model1)
str(model1)
or1<-mod1$coefficients[2,1]
se1<-mod1$coefficients[2,2]
p1<-mod1$coefficients[2,4]
row1<-cbind(or1, se1, p1)

#current smokers
current<-subset(data3, data3$smoking_status=="Current")
str(current) #N= 33362
model2<- glm(ibd ~ snp + ageatstart + sex, data=current, family=binomial(link="logit"))
mod2<-summary(model2)
or2<-mod2$coefficients[2,1]
se2<-mod2$coefficients[2,2]
p2<-mod2$coefficients[2,4]
row2<-cbind(or2, se2, p2)

#Former smokers
former<-subset(data3, data3 $smoking_status=="Previous")
str(former) #N= 118478
model3<- glm(ibd ~ snp + ageatstart + sex, data=former, family=binomial(link="logit"))
mod3<-summary(model3)
or3<-mod3$coefficients[2,1]
se3<-mod3$coefficients[2,2]
p3<-mod3$coefficients[2,4]
row3<-cbind(or3, se3, p3)

###UC
#ever smokers
model4<- glm(uc ~ snp + ageatstart + sex, data=ever, family=binomial(link="logit"))
mod4<-summary(model4)
or4<-mod4$coefficients[2,1]
se4<-mod4$coefficients[2,2]
p4<-mod4$coefficients[2,4]
row4<-cbind(or4, se4, p4)

#current smokers
model5<- glm(uc ~ snp + ageatstart + sex, data=current, family=binomial(link="logit"))
mod5<-summary(model5)
or5<-mod5$coefficients[2,1]
se5<-mod5$coefficients[2,2]
p5<-mod5$coefficients[2,4]
row5<-cbind(or5, se5, p5)

#Former smokers
model6<- glm(uc ~ snp + ageatstart + sex, data=former, family=binomial(link="logit"))
mod6<-summary(model6)
or6<-mod6$coefficients[2,1]
se6<-mod6$coefficients[2,2]
p6<-mod6$coefficients[2,4]
row6<-cbind(or6, se6, p6)

###Chrons
#ever smokers
model7<- glm(chrons ~ snp + ageatstart + sex, data=ever, family=binomial(link="logit"))
mod7<-summary(model7)
or7<-mod7$coefficients[2,1]
se7<-mod7$coefficients[2,2]
p7<-mod7$coefficients[2,4]
row7<-cbind(or7, se7, p7)

#current smokers
model8<- glm(chrons ~ snp + ageatstart + sex, data=current, family=binomial(link="logit"))
mod8<-summary(model8)
or8<-mod8$coefficients[2,1]
se8<-mod8$coefficients[2,2]
p8<-mod8$coefficients[2,4]
row8<-cbind(or8, se8, p8)

#Former smokers
model9<- glm(chrons ~ snp + ageatstart + sex, data=former, family=binomial(link="logit"))
mod9<-summary(model9)
or9<-mod9$coefficients[2,1]
se9<-mod9$coefficients[2,2]
p9<-mod9$coefficients[2,4]
row9<-cbind(or9, se9, p9)

######################################################################################
#8. Test Assosiation of SNP with IBD in non-smokers 
######################################################################################
###IBD
#never smokers
model10<- glm(ibd ~ snp + ageatstart + sex, data=never, family=binomial(link="logit"))
mod10<-summary(model10)
or10<-mod10$coefficients[2,1]
se10<-mod10$coefficients[2,2]
p10<-mod10$coefficients[2,4]
row10<-cbind(or10, se10, p10)

###UC
model11<- glm(uc ~ snp + ageatstart + sex, data=never, family=binomial(link="logit"))
mod11<-summary(model11)
or11<-mod11$coefficients[2,1]
se11<-mod11$coefficients[2,2]
p11<-mod11$coefficients[2,4]
row11<-cbind(or11, se11, p11)

###Chrons
model12<- glm(chrons ~ snp + ageatstart + sex, data=never, family=binomial(link="logit"))
mod12<-summary(model12)
or12<-mod12$coefficients[2,1]
se12<-mod12$coefficients[2,2]
p12<-mod12$coefficients[2,4]
row12<-cbind(or12, se12, p12)

#stick all of the log odds together
results<-rbind(row1, row2, row3, row10, row4, row5, row6, row11, row7, row8, row9, row12)
results<-data.frame(results)
names(results)<-c("b", "se", "pval")
str(results)

#convert to odds ratios
library(TwoSampleMR)
or<-generate_odds_ratios(results)
or

######################################################################################
#9. Plot the results
######################################################################################
library(forestplot)

#Betas - create dataframe
cochrane_from_rmeta <- 
  structure(list(
    mean  = c(NA, or$or), 
    lower = c(NA, or$or_lci95),
    upper = c(NA, or$or_uci95)), 
	.Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -13L), 
    class = "data.frame")

tabletext<-cbind(
  c(NA, "IBD", NA, NA, NA, "UC", NA,NA, NA, "Chrons", NA, NA, NA),
  c(NA, "Ever", "Current", "Previous", "Never", "Ever", "Current", "Previous", 
   "Never", "Ever", "Current", "Previous", "Never"),
  c("OR (95% CI)", "1.06 (0.99, 1.13)", "1.02 (0.87, 1.19)", "1.06 (0.99, 1.14)", "0.94 (0.87, 1.01)", "1.10 (1.02, 1.19)", "1.01 (0.82, 1.24)", "1.12 (1.03, 1.21)", "0.93 (0.85, 1.02)", "0.95 (0.85, 1.06)", "1.00 (0.81, 1.24)", "0.93 (0.82, 1.05)", "0.96 (0.86, 1.09)"), 
    c("P value", "0.10", "0.81", "0.09", "0.11", "0.01", "0.93", '0.01', "0.11", "0.34", "0.97", "0.25", '0.54'))
    
pdf.options(reset = TRUE, onefile = FALSE)
pdf("SNP_ibd_forest.pdf", width=10,height=6)
forestplot(tabletext, 
           cochrane_from_rmeta,new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,12)),
           zero=1,
           lineheight = unit(1, "cm"),
           		graphwidth=unit(4, "cm"),
           		boxsize=0.25,
           clip=c(0,2), 
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))
dev.off()


