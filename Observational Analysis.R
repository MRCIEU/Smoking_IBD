#Estimating the observed association between smoking and IBD diagnosis
#Robyn Wootton August 2019
#######################################################################################################################
#Contents
#1. Read in the IBD data
#2. Tidy up data - name variables appropriately
#3. Exclude participants who have withdrawn consent
#4. Read in SNP and other smoking variables
#5. Split into ever and never smokers
#6. Run  the observational associations (controlling for age, sex and SES)

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
#6. Run  the observational associations (controlling for age, sex and SES)
######################################################################################
#Flip the levels of the smoking variables the right way so that it increases if smoking increases
gen<-data3
#Never/Ever
levels(gen$nevev)
gen$nevev<-as.numeric(gen$nevev)
table(gen$nevev)
gen$nevev[gen$nevev==2]<-0
table(gen$nevev)
gen$nevev<-factor(gen$nevev, levels = c(0,1),labels = c("Never", "Ever"))
table(gen$nevev)
levels(ever$smoking_status)

#Current/Previous
levels(gen$smoking_status)
table(gen$smoking_status)
gen$smoking_status <-as.numeric(gen$smoking_status)
table(gen$smoking_status)
gen$smoking_status[gen$smoking_status ==2]<-0
gen$smoking_status[gen$smoking_status ==1]<-2
gen$smoking_status[gen$smoking_status ==3]<-1
table(gen$smoking_status)
gen$smoking_status <-factor(gen$smoking_status, levels = c(0,1,2),labels = c("Never", "Previous", "Current"))
table(gen$smoking_status)

#create ever again with these changes
ever<-subset(gen, gen$nevev=="Ever")
str(ever) #N= 151822

##Smoking status (ever, never)
#IBD
model1<- glm(ibd ~ nevev + ageatstart + sex +SES, data=gen, family=binomial(link="logit"))
mod1<-summary(model1)
or1<-mod1$coefficients[2,1]
se1<-mod1$coefficients[2,2]
p1<-mod1$coefficients[2,4]
n1<-mod1$df[2]
row1<-cbind(or1, se1, p1, n1)
#Chrons
model2<- glm(chrons ~ nevev + ageatstart + sex +SES, data=gen, family=binomial(link="logit"))
mod2<-summary(model2)
or2<-mod2$coefficients[2,1]
se2<-mod2$coefficients[2,2]
p2<-mod2$coefficients[2,4]
n2<-mod2$df[2]
row2<-cbind(or2, se2, p2, n2)
#UC
model3<- glm(uc ~ nevev + ageatstart + sex +SES, data=gen, family=binomial(link="logit"))
mod3<-summary(model3)
or3<-mod3$coefficients[2,1]
se3<-mod3$coefficients[2,2]
p3<-mod3$coefficients[2,4]
n3<-mod3$df[2]
row3<-cbind(or3, se3, p3, n3)

##Smoking status (current/former within ever smokers)
#IBD
model4<- glm(ibd ~ smoking_status + ageatstart + sex +SES, data=ever, family=binomial(link="logit"))
mod4<-summary(model4)
or4<-mod4$coefficients[2,1]
se4<-mod4$coefficients[2,2]
p4<-mod4$coefficients[2,4]
n4<-mod4$df[2]
row4<-cbind(or4, se4, p4, n4)
#crohns
model5<- glm(chrons ~ smoking_status + ageatstart + sex +SES, data=ever, family=binomial(link="logit"))
mod5<-summary(model5)
or5<-mod5$coefficients[2,1]
se5<-mod5$coefficients[2,2]
p5<-mod5$coefficients[2,4]
n5<-mod5$df[2]
row5 <- cbind(or5, se5, p5, n5)
#uc
model6<- glm(uc ~ smoking_status + ageatstart + sex +SES, data=ever, family=binomial(link="logit"))
mod6<-summary(model6)
or6<-mod6$coefficients[2,1]
se6<-mod6$coefficients[2,2]
p6<-mod6$coefficients[2,4]
n6<-mod6$df[2]
row6<-cbind(or6, se6, p6, n6)

##CPD (within ever smokers)
#ibd
model7<- glm(ibd ~ cpd + ageatstart + sex +SES, data=ever, family=binomial(link="logit"))
mod7<-summary(model7)
or7<-mod7$coefficients[2,1]
se7<-mod7$coefficients[2,2]
p7<-mod7$coefficients[2,4]
n7<-mod7$df[2]
row7<-cbind(or7, se7, p7, n7)
#chrons
model8<- glm(chrons ~ cpd + ageatstart + sex +SES, data=ever, family=binomial(link="logit"))
mod8<-summary(model8)
or8<-mod8$coefficients[2,1]
se8<-mod8$coefficients[2,2]
p8<-mod8$coefficients[2,4]
n8<-mod8$df[2]
row8<-cbind(or8, se8, p8, n8)
#uc
model9<- glm(uc ~ cpd + ageatstart + sex +SES, data=ever, family=binomial(link="logit"))
mod9<-summary(model9)
or9<-mod9$coefficients[2,1]
se9<-mod9$coefficients[2,2]
p9<-mod9$coefficients[2,4]
n9<-mod9$df[2]
row9<-cbind(or9, se9, p9, n9)

###CSI
#merge in CSI
csi<-read.table("lifetimesmoking_phen_revision1.txt", header=T)
str(csi)
csi$csi <- as.numeric(csi$csi)
test<-c()
test <- subset(csi, csi$csi>=0.99)
test2 <- subset(test, test$csi<1.01)

#Merge into the main dataset
gen_csi<-merge(gen, csi, by.x="FID", by.y="IID", all.x=T)

#ibd
model10<- glm(ibd ~ csi + ageatstart + sex +SES, data=gen_csi, family=binomial(link="logit"))
mod10<-summary(model10)
or10<-mod10$coefficients[2,1]
se10<-mod10$coefficients[2,2]
p10<-mod10$coefficients[2,4]
n10<-mod10$df[2]
row10<-cbind(or10, se10, p10, n10)
#chrons
model11<- glm(chrons ~ csi + ageatstart + sex +SES, data= gen_csi, family=binomial(link="logit"))
mod11<-summary(model11)
or11<-mod11$coefficients[2,1]
se11<-mod11$coefficients[2,2]
p11<-mod11$coefficients[2,4]
n11<-mod11$df[2]
row11<-cbind(or11, se11, p11, n11)
#uc
model12<- glm(uc ~ csi + ageatstart + sex +SES, data= gen_csi, family=binomial(link="logit"))
mod12<-summary(model12)
or12<-mod12$coefficients[2,1]
se12<-mod12$coefficients[2,2]
p12<-mod12$coefficients[2,4]
n12<-mod12$df[2]
row12<-cbind(or12, se12, p12, n12)

#stick all of the log odds together
results<-rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9, row10, row11, row12)
results<-data.frame(results)
exposure<-c("Ever_Never","Ever_Never","Ever_Never", "Current_Former", "Current_Former","Current_Former","cpd", "cpd","cpd","csi", "csi","csi")
outcome<-c("IBD", "Crohns", "UC", "IBD", "Crohns", "UC","IBD", "Crohns", "UC","IBD", "Crohns", "UC")
results2<-cbind(exposure, outcome, results)
names(results2)<-c("exposure", "outcome", "b", "se", "pval", "N")
str(results2)

#convert to odds ratios
library(TwoSampleMR)
or<-generate_odds_ratios(results2)
or

#save the dataset
write.csv(or, "June2019_IBD_smoking_observed.csv", quote=F, row.names=F)
