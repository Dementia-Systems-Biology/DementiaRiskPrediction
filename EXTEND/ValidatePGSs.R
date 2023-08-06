# File: ValidatePGSs_EXTEND.R
# Author: Jarno Koetsier
# Date: August 6, 2023

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)
library(caret)

# Load data
load("EXTEND/Data/df_list.RData")         # PGSs
load("EXTEND/Data/metaData_ageFil.RData") # Meta data

# Prepare data
df_Result_PGS <- df_list$bayesr
rownames(dat) <- dat$ID
dat_fil <- dat[rownames(df_Result_PGS),]

# Systolic blood pressure
R2(pred = df_Result_PGS$SBPauto[!is.na(as.numeric(dat_fil$MeanSysBP))], 
   obs = as.numeric(dat_fil$MeanSysBP)[!is.na(as.numeric(dat_fil$MeanSysBP))])

# Total cholesterol
R2(pred = df_Result_PGS$TC[!is.na(as.numeric(dat_fil$Chol))], 
   obs = as.numeric(dat_fil$Chol)[!is.na(as.numeric(dat_fil$Chol))])

# BMI
R2(pred = df_Result_PGS$BMI[!is.na(as.numeric(dat_fil$BMI))], 
   obs = as.numeric(dat_fil$BMI)[!is.na(as.numeric(dat_fil$BMI))])

# HDL cholesterol
R2(pred = df_Result_PGS$HDL[!is.na(as.numeric(dat_fil$HDL))], 
   obs = as.numeric(dat_fil$HDL)[!is.na(as.numeric(dat_fil$HDL))])

# Education
Education <- ifelse(dat_fil$None.of.the.above == 1,"Yes","No")
roc(Education[!is.na(Education)],
    df_Result_PGS$EA22[!is.na(Education)])

# Physical (in)acitvity
Physical <- ifelse(dat_fil$Exercise.increased.pulse.more.than.2halfhrsawk == 1,"No","Yes")
roc(Physical[!is.na(Physical)],
    df_Result_PGS$MVPA[!is.na(Physical)])

# Dietary intake (healthy diet)
Diet <- ifelse(as.numeric(dat_fil$Fruit) + as.numeric(dat_fil$Vegtables) > 3, "No", "Yes")
roc(Diet[!is.na(Diet)],
    df_Result_PGS$DC2[!is.na(Diet)])

# Depression
Depression <- ifelse(dat_fil$Depression==1,"Yes","No")
roc(Depression[!is.na(Depression)],
    df_Result_PGS$MDD[!is.na(Depression)])

# Type II diabetes
Diabetes <-ifelse(dat_fil$T2.Diabetes==1,"Yes","No")
roc(Diabetes[!is.na(Diabetes)],
    df_Result_PGS$T2D[!is.na(Diabetes)])

# Heart disease
HeartDisease <- ifelse(dat_fil$Heart.Disease==1,"Yes","No")
roc(HeartDisease[!is.na(HeartDisease)],
    df_Result_PGS$CAD[!is.na(HeartDisease)])

# Low-moderate alcohol intake
dat$LtoMAlcohol <- NULL
dat[dat$How.often.do.you.drink.alcohol==0 ,"LtoMAlcohol"] <- 1
dat[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==3 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==4 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==5 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1

dat[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-1

dat[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-1

# High alcohol intake
dat[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==3  &  dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==3 &   dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==2 &   dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==3 &   dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day=="." &    dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0
dat_fil <- dat[rownames(df_Result_PGS),]

Alcohol <- ifelse(dat_fil$LtoMAlcohol==1,"Yes","No")
roc(Alcohol[!is.na(Alcohol)],
    df_Result_PGS$Alcohol[!is.na(Alcohol)])