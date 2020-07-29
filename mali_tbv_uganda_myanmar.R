---
title: "NIH - Mali TBV Trial Comparison Studies"
---
  
#Emily Reichert (PATH), July 2019
#This script pulls in data from two studies looking at the Alere PF uRDT conducted in 
  #Myanmar (n= 1847) and Uganda (n = 607)
  #Permission was obtained from study authors and dataset owners to use this data
  #in comparison to Mali TBV Trial uRDT performance.
  
#Note: Must run 'mali_tbv_exploratory.R' script before this to clean and prepare 'NIH' dataset.

library(tidyr)
library(readr)
library(cutpointr)
library(readxl)
library(dplyr)
library(ggplot2)

#Retrieve two datasets for comparison
setwd("~/Documents/2019 PATH research")
SMRU_Myanmar <- read_excel("smru_myanmar.xlsx")
BIOME1_uganda_1_ <- read_csv("BIOME1_uganda.csv")

#I. Format SMRU Myanmar Study data
#Make HRP2 and LDH numeric, treating concentrations beyond LOQ as LLOQ/2 and ULOQ*2
SMRU_Myanmar$HRP2_pg_ml[SMRU_Myanmar$HRP2_pg_ml == '< 1.07'] <- 1.07/2
SMRU_Myanmar$HRP2_pg_ml[SMRU_Myanmar$HRP2_pg_ml == "> 16500.00"] <- 16500*2
SMRU_Myanmar$HRP2_pg_ml <- as.numeric(SMRU_Myanmar$HRP2_pg_ml)

SMRU_Myanmar$LDH_Pan_pg_ml[SMRU_Myanmar$LDH_Pan_pg_ml == "< 14.41"] <- 14.41/2
SMRU_Myanmar$LDH_Pan_pg_ml[SMRU_Myanmar$LDH_Pan_pg_ml == "> 525700.00"] <- 525700*2
SMRU_Myanmar$LDH_Pan_pg_ml <- as.numeric(SMRU_Myanmar$LDH_Pan_pg_ml)

#make binary 0/1 variable for HSRDT result (HSRDT = highly-sensitive RDT = uRDT)
SMRU_Myanmar$HSRDT_num <- NA
SMRU_Myanmar$HSRDT_num[SMRU_Myanmar$HSRDT == "Neg"] <- 0
SMRU_Myanmar$HSRDT_num[SMRU_Myanmar$HSRDT == "PF"] <- 1

#II. Format BIOME Uganda Study Data
#make binary 0/1 variables for HRP2 and HSRDT results in Uganda
BIOME1_uganda_1_$HRP2_result <- as.numeric(BIOME1_uganda_1_$`Quansys HRP2 designation` == "Pos")

BIOME1_uganda_1_$HS_RDT_Result <- as.numeric(BIOME1_uganda_1_$`HS RDT pos/neg` == "Pos")

#Make HRP2 numeric,treating concentrations beyond LOQ as LLOQ/2 and ULOQ*2
BIOME1_uganda_1_$`Quansys HRP2 avg conc (pg/mL)`[BIOME1_uganda_1_$`Quansys HRP2 avg conc (pg/mL)`
                                                 == '<0.1'] <- 0.1/2
BIOME1_uganda_1_$`Quansys HRP2 avg conc (pg/mL)`[BIOME1_uganda_1_$`Quansys HRP2 avg conc (pg/mL)`
                                                 == '>14600'] <- 14600*2
BIOME1_uganda_1_$`Quansys HRP2 avg conc (pg/mL)`<- 
  as.numeric(BIOME1_uganda_1_$`Quansys HRP2 avg conc (pg/mL)`)

#III. Compare HRP2 distributions, uRDT performance by study
hist(log10(SMRU_Myanmar$HRP2_pg_ml))
hist(log10(BIOME1_uganda_1_$`Quansys HRP2 avg conc (pg/mL)`))
hist(log10(NIH$HRP2_pg_ml))

#Calculate HRP2 concentration where 50% uRDT+ achieved, Mali
mod_hrp2 = glm(uRDT_result ~ log10(HRP2_pg_ml), family = binomial, data = NIH)
lod50 <- dose.p(mod_hrp2, p=0.50)  #Calculates LD50
lod50.ci <- lod50 + attr(lod50, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1) #calculate 95% CI
lod50.ci <- 10^lod50.ci

#Calculate HRP2 concentration where 50% uRDT+ achieved, Myanmar
mod_hrp2 = glm(HSRDT_num ~ log10(HRP2_pg_ml), family = binomial, data = sub_myanmar)
lod50 <- dose.p(mod_hrp2, p=0.50)  #Calculates LD50
lod50.ci <- lod50 + attr(lod50, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1) #calculate 95% CI
lod50.ci <- 10^lod50.ci

#Calculate HRP2 concentration where 50% uRDT+ achieved, Uganda
colnames(BIOME1_uganda_1_)[colnames(BIOME1_uganda_1_)=="Quansys HRP2 avg conc (pg/mL)"] <- 
  "HRP2_pg_ml"
mod_hrp2 = glm(HS_RDT_Result ~ log10(HRP2_pg_ml), 
               family = binomial, data = BIOME1_uganda_1_)
lod50 <- dose.p(mod_hrp2, p=0.50)  #Calculates LD50
lod50.ci <- lod50 + attr(lod50, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1) #calculate 95% CI
lod50.ci <- 10^lod50.ci

