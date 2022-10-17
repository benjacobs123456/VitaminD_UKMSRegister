# Vitamin D in MS

This readme contains code used to perform analysis for Kuri et al 2022, where we look at determinants of vitamin D status in people with MS within the UK MS Register.

Code: Ben Jacobs & Ashvin Kuri, 17th Oct 2022


## Download vitamin D GWAS summary stats

````unix
cd /data/Wolfson-UKBB-Dobson/vitd_ukb_20_12_21

# download sunlight from https://drive.google.com/file/d/0BzYDtCo_doHJRUhaWDM4b2dxQVU/view?usp=sharing&resourcekey=0-ZDi12Q_y9uV5NQRj0eXYGw
````

## Process vit D GWAS summary stats
````R
library(dplyr)
library(readr)

sunlight = read_table2("25HydroxyVitaminD_QC.METAANALYSIS1.txt")

# extract SNPs of interest
snps = c("rs2282679",
"rs12785878",
"rs10741657",
"rs17216707",
"rs6538691",
"rs8018720")

sunlight = sunlight %>% filter(MarkerName %in% snps)

sunlight = sunlight %>% select(MarkerName,Allele1,Effect)
sunlight$Allele1 = toupper(sunlight$Allele1)
write_tsv(sunlight,"sunlight.tsv",col_names=F)
````

## Apply vitamin D genetic risk score to UK Biobank

````unix
# load plink
module load plink/2.0-20200328

# extract SNPs
# rs6538691 as proxy for AMDHD1 SNP
for i in 4 11 12 14 20
  do
    plink2 --pfile \
    ../imputed_ukb_genotypes/plink2_files/chr_$i \
    --extract sunlight.tsv \
    --make-pgen --out chr$i
  done

# apply score to imputed genotypes
for i in 4 11 12 14 20
  do
    plink2 --pfile chr$i \
    --score sunlight.tsv variance-standardize cols=+scoresums \
    --out scores_chr$i
  done
````

## Examine validity of genetic risk score in UKB

````R
# libraries
library(dplyr)
library(readr)
setwd("/data/Wolfson-UKBB-Dobson/vitd_ukb_20_12_21")

# read in prs scores
score_df = lapply(c("4","11","12","14","20"),function(i){
  scores = read_tsv(paste0("scores_chr",i,".sscore"))
  scores
})
score_df = do.call("rbind",score_df)

# sum over alleles
score_df = score_df %>% group_by(IID) %>% summarise(GRS = sum(SCORE1_SUM))

# z score
z_normalise = function(x){
 z = (x - mean(x,na.rm=TRUE) ) / sd(x,na.rm=TRUE)
 return(z)
}
score_df$GRS = z_normalise(score_df$GRS)

# read in ukb phenotype data
df = read_tsv("../ukb_pheno_911/ukb_pheno_final_MS_2312")
biochem = read_tsv("../ukb_pheno_2410/ukb_coded_pheno_final.tsv",col_types=cols_only(EID=col_double(),`Vitamin D.0.0`=col_double()))
combo = left_join(df,biochem,by="EID")
rm(df)
rm(biochem)

# filter to high quality vit D recordings
combo = combo %>% filter(`Vitamin D reportability.0.0` == "Reportable at assay and after aliquot correction, if attempted")

# define supplementers
supp = combo %>% filter_at(vars(contains("Vitamin and mineral")),any_vars(.=="Vitamin D")) %>% mutate("VitD_supp"="Yes")
no_supp = combo %>% filter(!EID %in% supp$EID) %>% mutate("VitD_supp"="No")
combo = bind_rows(supp,no_supp)

# look at GRS
combo = combo %>% left_join(score_df %>% rename("EID" = IID) ,by="EID")

# normalise vd
combo$`Vitamin D.0.0` = z_normalise(combo$`Vitamin D.0.0`)

# filter to european-ancestry
combo = combo %>% filter(`Genetic ethnic grouping.0.0` == "Caucasian")

# plots
library(Hmisc)
combo$GRS_decile = cut2(combo$GRS,g=10)
ggplot(combo,aes(GRS_decile,`Vitamin D.0.0`))+geom_violin()+geom_point(alpha=0.1)
ggplot(combo,aes(GRS,`Vitamin D.0.0`))+geom_point()+geom_smooth(method="lm")

# models
model = lm(data=combo,
`Vitamin D.0.0` ~
`Age at recruitment.0.0` + Sex.0.0 + VitD_supp + GRS)
summary(model)

# repeat in non-supplementers
model = lm(data=combo %>% filter(VitD_supp=="No"),
`Vitamin D.0.0` ~
`Age at recruitment.0.0` + Sex.0.0 + GRS)
summary(model)

# repeat in supplementers
model = lm(data=combo %>% filter(VitD_supp=="Yes"),
`Vitamin D.0.0` ~
`Age at recruitment.0.0` + Sex.0.0 + GRS)
summary(model)

model = lm(data=combo,
`Vitamin D.0.0` ~
`Age at recruitment.0.0` + Sex.0.0 + GRS)
summary(model)

# basic demographics
combo %>% count(`Sex.0.0`) %>% mutate(n/sum(n))
combo %>% count(`VitD_supp`) %>% mutate(n/sum(n))
combo %>% count(`Ethnic background.0.0`) %>% mutate(n/sum(n))
combo %>% count(`Country of birth (UK/elsewhere).0.0`) %>% mutate(n/sum(n))
combo %>% summarise_at(vars(`Age at recruitment.0.0`,`Townsend deprivation index at recruitment.0.0`),.funs=c("median","IQR"),na.rm=T)
````

## UK MS Register analysis

````R
# libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(gridExtra)

# set wd
setwd("S:/VitaminDinMS - Vitamin D in MS/Ashvin_Vit_D/")

###################
# read in data
###################

vd_geno = read_csv("vitd_bj_1812.csv")
vd_pheno = read_csv("df_vd_return_for_BEN.csv")

###################
# clean data
###################

# clean genotype file
glimpse(vd_geno)
vd_geno
nrow(vd_geno)

# 684 readings
vd_geno$VitDSerum = as.numeric(vd_geno$VitDSerum)
vd_geno = vd_geno %>% filter(!is.na(VitDSerum))
nrow(vd_geno)

# 579 readings after removal of NAs
summary(vd_geno$VitDSerum)

# check IDs are unique
vd_geno %>% distinct(VDQM_ID) %>% nrow() # they are

# split into alleles

for(colname in colnames(vd_geno)[-c(1:2)] ){
  message("Splitting alleles for ", colname)
  vd_geno <<- vd_geno %>%
    separate(colname,
             pattern="/",
             into=c(paste0(colname,"_A1"),paste0(colname,"_A2")))

}

# merge with pheno file
vd_pheno
nrow(vd_pheno)
# filter to those in geno file
vd_pheno = vd_pheno %>% filter(VDQM_ID %in% vd_geno$VDQM_ID) # 569 remaining
# merge
vd_df = vd_geno %>% left_join(vd_pheno,by="VDQM_ID")

# clean phenotypes
glimpse(vd_df)
# cc status
table(vd_df$Group_1MS_2Con)
vd_df$Group_1MS_2Con = as.factor(vd_df$Group_1MS_2Con)
vd_df = vd_df %>% mutate(case_control_status = ifelse(Group_1MS_2Con==1,"MS","Control"))
table(vd_df$case_control_status) # 259 controls, 320 MS
vd_df = vd_df %>% select(-Group_1MS_2Con)
# vd supplements
table(vd_df$vitd_supplement_yn)
vd_df = vd_df %>% mutate(vitd_supplement_yn = ifelse(vitd_supplement_yn =="NULL",NA, vitd_supplement_yn)) # recode NULL as NA
vd_df = vd_df %>% mutate(vd_supp = ifelse(vitd_supplement_yn=="0" & !is.na(vitd_supplement_yn),"No","Yes"))
table(vd_df$vd_supp)
vd_df = vd_df %>% select(-vitd_supplement_yn)

vd_df$vitdIUperDay = as.numeric(vd_df$vitdIUperDay)

# bmi
vd_df$bmi = as.numeric(vd_df$bmi)
summary(vd_df$bmi) # 17 - 51 which is reasonable

# month
vd_df$Month_blood = as.numeric(vd_df$Month_blood)
vd_df = vd_df %>% mutate(season_of_test = ifelse(Month_blood <=3 | Month_blood >=10,"Winter (Oct - March)","Summer (April - Sept)"))
table(vd_df$season_of_test)

# sex
vd_df$Gender[vd_df$Gender=="NULL"] = NA # recode NULL as NA
table(is.na(vd_df$Gender)) # only 3 missing
vd_df = vd_df %>% filter(!is.na(Gender))
vd_df$Gender = factor(vd_df$Gender)
vd_df = vd_df %>%
  mutate(Gender = ifelse(Gender == "0","F","M"))
table(vd_df$Gender,vd_df$case_control_status)

# age
vd_df$Age_study = as.numeric(vd_df$Age_study)
table(is.na(vd_df$Age_study)) # only 29 missing
vd_df = vd_df %>% filter(!is.na(Age_study))
summary(vd_df$Age_study) # 25 - 85, which is reasonable

# latitude
vd_df$Latitude = as.numeric(vd_df$Latitude)
table(is.na(vd_df$Latitude)) # lots of missingness, 261
summary(vd_df$Latitude) # 50 - 60, which is pretty much cornwall to shetland

# final n
nrow(vd_df)

###################
# descriptive stats
###################

summary(vd_df)

med_iqr = function(x){
  median = median(x,na.rm=TRUE)
  iqr = IQR(x,na.rm=TRUE)
  paste0(round(median,2)," (",round(iqr,2),")")
}
cont_vars = vd_df %>% group_by(case_control_status) %>%
  summarise_at(.vars = c("Age_study",
                         "bmi",
                         "Latitude",
                         "vitdIUperDay",
                         "VitDSerum"),
               .fun="med_iqr") %>% t() %>% data.frame()
colnames(cont_vars) = c("Control","MS")
cont_vars$variable = rownames(cont_vars)
cont_vars = cont_vars[-1,]
rownames(cont_vars) = NULL

cat_vars = c("Gender","vd_supp","season_of_test")
cat_var_summ = lapply(cat_vars,function(x){
tbl = table(vd_df$case_control_status,vd_df[[x]])  
prop_tbl = round(tbl / rowSums(tbl) * 100,2)
joint_tbl = cbind(t(tbl), t(prop_tbl))
paste0(tbl," (",prop_tbl,"%)")
data.frame(
  "Control"=paste0(joint_tbl[,1]," (",joint_tbl[,3],"%)" ),
  "MS"=paste0(joint_tbl[,2]," (",joint_tbl[,4],"%)" ),
  "variable"=paste0(x,"_",colnames(tbl)))
})

cat_var_summ_tbl = do.call("rbind",cat_var_summ)
cat_var_summ_tbl = data.frame(cat_var_summ_tbl)

summary_df = bind_rows(cat_var_summ_tbl,cont_vars)

# save to file
write_csv(summary_df,"./outputs/demographics.csv")

# simple plots

# histograms for continuous variables
cont_vars = c("Age_study",
              "bmi",
              "Latitude",
              "vitdIUperDay",
              "VitDSerum")
labels = c("Age at study entry (years)","BMI","Latitude","Vitamin D supplemented per day (IU)","Raw serum vitamin D (nM)")
hists = lapply(c(1:length(cont_vars)),function(i){
  var = cont_vars[i]
  label = labels[i]
  message("plotting ",var)
  p=ggplot(vd_df,aes(vd_df[[var]]))+
    geom_histogram()+
    labs(x=label)+
    theme_bw()
  print(p)
})


png("./outputs/histograms.png",res=300,units="in",height=8,width=8)
do.call("grid.arrange",c(hists,ncol=2))
dev.off()

vln_ms_v_control = lapply(c(1:length(cont_vars)),function(i){
  var = cont_vars[i]
  label = labels[i]
  message("plotting ",var)
  p=ggplot(vd_df,aes(case_control_status,fill=case_control_status,vd_df[[var]]))+
    geom_violin()+
    geom_jitter(alpha=0.2)+
    labs(y=label,x="MS status",fill="MS status")+
    theme_bw()
  print(p)
})

png("./outputs/violin_plots.png",res=300,units="in",height=8,width=8)
do.call("grid.arrange",c(vln_ms_v_control,ncol=2))
dev.off()


###################
# determinants of vit D
###################

# Normalise vit D
shapiro.test(vd_df$VitDSerum) # not normal
vd_df = vd_df %>% mutate(vitD_Z = ( VitDSerum - mean(VitDSerum) ) / sd (VitDSerum) )

# look at predictors of VitD
predictors = c("bmi","Latitude","Gender","Age_study","case_control_status","vd_supp","season_of_test")
labels = c("BMI","Latitude","Gender","Age","MS status","Vitamin D supplementation?","Season of test")

model_outputs = list()
plots = lapply(c(1:length(predictors)),function(i){
  x = predictors[i]
  label = labels[i]
  filtered_df = vd_df %>% filter(!is.na(vd_df[[x]]))
  message("processing ",x)

  if(is.numeric(filtered_df[[x]])){
    message("Numeric - make scatter plot")
    p <<- ggplot(filtered_df,aes(filtered_df[[x]],vitD_Z))+geom_point()+geom_smooth(method="lm")+labs(x=label,y="Serum Vitamin D Z score")
  } else {
    message("Not numeric - make box plot")
    p <<- ggplot(filtered_df,aes(filtered_df[[x]],vitD_Z))+geom_violin()+geom_jitter(alpha=0.2)+labs(x=label,y="Serum Vitamin D Z score")
  }
  print(p)
})  

png(paste0("./outputs/vitd_predictor_plots.png"),res=300,units="in",height=8,width=12)
do.call("grid.arrange", c(plots,ncol=3))
dev.off()

# models
model_outputs = list()
for(i in c(1:length(predictors))){
  x = predictors[i]
  model=lm(data=filtered_df,
     vitD_Z ~ filtered_df[[x]])
  output = c(x,rownames(summary(model)$coefficients)[2],summary(model)$coefficients[2,])
  model_outputs[[i]] = output
}
library(stringr)

model_df = do.call("rbind",model_outputs) %>% data.frame()
model_df = model_df %>% mutate(V2 = str_remove(V2,"filtered_df\\[\\[x]]")) %>%
  mutate(variable = paste0(V1,"_",V2)) %>%
  select(-V1,-V2,-t.value) %>%
  rename("beta" = Estimate,
         "se" = Std..Error,
         "pval" = Pr...t..) %>%
  select(variable,beta,se,pval) %>%
  mutate(beta = as.numeric(as.character(beta))) %>%
  mutate(se = as.numeric(as.character(se))) %>%
  mutate(pval = as.numeric(as.character(pval))) %>%
  mutate(lower_ci = beta - 1.96*se) %>%
  mutate(upper_ci = beta + 1.96*se) %>%
  mutate(fdr = p.adjust(pval,method="fdr"))

model_df = model_df %>% arrange(fdr)
write_csv(model_df,"./outputs/demographic_predictors_vitd.csv")
# in the cohort predictors of higher vit D were
# supplementation
# having MS
# lower BMI
# being older

# multivariable model including all significant predictors
model=lm(data=vd_df,
         vitD_Z ~ vd_supp + case_control_status + bmi + Age_study
           )
summ_model = data.frame(summ_model$coefficients) %>% mutate(var = rownames(summ_model$coefficients))
write_csv(summ_model,"./outputs/multivariable_model.csv")
message("Adjusted R2 = ",summary(model)$adj.r.squared)


###################
# genetics
###################

sunlight = read_csv("sunlight_snps.csv",col_types=cols(
  Gene = col_character(),
  SNP = col_character(),
  hg19 = col_character(),
  A1 = col_character(),
  A2 = col_character(),
  AF = col_double(),
  Beta = col_double(),
  se = col_double()
))

vd_geno = vd_df %>% select(1,contains("A1"), contains("A2"))
vd_geno = na.omit(vd_geno) # remove any missing genotypes

score_results = list()
allele_stats = list()
ea_dose_list = list()
for(i in c(1:length(sunlight$Gene))){
  gene = sunlight$Gene[i]
  message("Doing ",gene)
  ea = sunlight[sunlight$Gene==gene,]$A1
  beta = sunlight[sunlight$Gene==gene,]$Beta
  genos = vd_geno %>% select(contains(gene))
  ea_doses = cbind(genos[,1]==ea,genos[,2]==ea) %>% rowSums()
  scores = ea_doses * beta
  score_results[[i]] = data.frame(VDQM_ID = vd_geno$VDQM_ID,scores)
  ea_freq = sum(ea_doses)/(length(ea_doses) * 2)
  ea_stats = data.frame(gene,ea,ea_freq)
  allele_stats[[i]] = ea_stats
  ea_dose_list[[i]] = data.frame(VDQM_ID = vd_geno$VDQM_ID,gene,ea_doses)

}

# per SNP tests ~ vD
genotypes = do.call("rbind",ea_dose_list)
genotypes = genotypes %>% left_join(vd_df %>% select(VDQM_ID,vitD_Z,Age_study,Gender),by="VDQM_ID")
coefs = list()
for(i in c(1:length(sunlight$Gene))){
  test_gene = sunlight$Gene[i]
  genotypes_df = genotypes %>% filter(gene == test_gene)
  model = lm(data=genotypes_df,
             vitD_Z ~ Age_study + Gender + ea_doses
             )
  coefs[[i]] = cbind(test_gene,summary(model)$coefficients)
}
coefs = do.call("rbind",coefs)
coefs = data.frame(coefs[rownames(coefs)=="ea_doses",])

# check eafs
allele_df = do.call("bind_rows",allele_stats)
score_df = do.call("bind_rows",score_results)
score_df = score_df %>% group_by(VDQM_ID) %>% summarise(GRS = sum(scores))

# normalise
score_df = score_df %>% mutate(GRS_Z = ( GRS - mean(GRS) ) / sd(GRS) )

# join with main dataset
vd_df = vd_df %>% left_join(score_df,by="VDQM_ID") %>% select(-contains("A1"),-contains("A2"))

###################
# GRS ~ vd
###################

p1=ggplot(vd_df,aes(GRS_Z,vitD_Z))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~vd_supp)+
  theme_classic()+
  labs(x="GRS Z score",y="Vitamin D Z score",facet="Supplementation?")

p2=ggplot(vd_df,aes(GRS_Z,vitD_Z))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~case_control_status)+
  theme_classic()+
  labs(x="GRS Z score",y="Vitamin D Z score",facet="MS status")

png("./outputs/grs_vitd.png",res=300,units="in",height=8,width=12)
grid.arrange(p1,p2,ncol=2)
dev.off()

# models

# just grs
model=lm(data=vd_df,
         vitD_Z ~ GRS_Z
)
summary(model)
message("Adjusted R2 = ",summary(model)$adj.r.squared) # weak association


# controlling for covariates
model=lm(data=vd_df,
         vitD_Z ~ vd_supp + case_control_status + bmi + Age_study + GRS_Z
)
summary(model) # after adjusting for covariates it does predict vit D status

summ_model= summary(model)
anova(model,test="Chisq") # LR test confirms effect P < 0.05
summ_model = data.frame(summ_model$coefficients) %>% mutate(var = rownames(summ_model$coefficients))
write_csv(summ_model,"./outputs/grs_multivariable_model.csv")
message("Adjusted R2 = ",summary(model)$adj.r.squared)

# controlling for covariates w sex
model=lm(data=vd_df,
         vitD_Z ~ vd_supp + case_control_status + bmi + Age_study + GRS_Z + Gender
)
summary(model) # after adjusting for covariates it does predict vit D status

summ_model= summary(model)
anova(model,test="Chisq") # LR test confirms effect P < 0.05
summ_model = data.frame(summ_model$coefficients) %>% mutate(var = rownames(summ_model$coefficients))
write_csv(summ_model,"./outputs/grs_multivariable_model_with_sex.csv")
message("Adjusted R2 = ",summary(model)$adj.r.squared)

# show that ms patients supplement with higher doses
vd_df %>% group_by(case_control_status) %>% summarise(mean(vitdIUperDay,na.rm=TRUE))

# repeat with only those with supplementing dose
model=lm(data=vd_df %>% filter(!is.na(vitdIUperDay)),
         vitD_Z ~ vitdIUperDay + case_control_status + bmi + Age_study + GRS_Z
)
summary(model) # after adjusting for covariates it does predict vit D status

summ_model= summary(model)
anova(model,test="Chisq") # LR test confirms effect P < 0.05
summ_model = data.frame(summ_model$coefficients) %>% mutate(var = rownames(summ_model$coefficients))
write_csv(summ_model,"./outputs/grs_multivariable_model_supplementing_dose.csv")
message("Adjusted R2 = ",summary(model)$adj.r.squared)


# interaction
model=lm(data=vd_df,
         vitD_Z ~ case_control_status + bmi + Age_study + GRS_Z * vd_supp
)
summary(model) # after adjusting for covariates it does predict vit D status
anova(model,test="Chisq")

# interaction with dose - there is a positive interaction
model=lm(data=vd_df,
         vitD_Z ~ case_control_status + bmi + Age_study + GRS_Z * vitdIUperDay
)
summary(model) # after adjusting for covariates it does predict vit D status
anova(model,test="Chisq")

# dose of supplementation
cor.test(vd_df$vitdIUperDay,vd_df$vitD_Z,na.rm=TRUE)

# stratified models
# just grs
model=lm(data=vd_df %>% filter(vd_supp == "No"),
         vitD_Z ~ GRS_Z
)
summary(model)
message("Adjusted R2 = ",summary(model)$adj.r.squared) # weak association

model=lm(data=vd_df %>% filter(vd_supp == "No"),
         vitD_Z ~ Gender + Age_study + GRS_Z
)
summary(model) # after adjusting for covariates it does predict vit D status
anova(model,test="Chisq")
````
