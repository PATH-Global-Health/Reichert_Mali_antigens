---
title: "NIH - Mali TBV Trial Exploratory Data Cleaning & Analysis"
---

#This script cleans and outputs initial statistics from the NIH Mali
#TBV Trial dataset.
#by Emily Reichert (PATH) June 2019

setwd("~/Documents/2019 PATH research")

#load necessary libraries
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)
library(EnvStats)
library(binom)
library(reportROC)

#import dataset of WB specimen data (n = 622)
NIH <- read_csv("NIH.csv")
NIH <- tbl_df(NIH)

#I. Data Cleaning

#Make binary uRDT result variable
#defines discordant RDT results as NA
NIH <- NIH %>%
  mutate(uRDT_result = NULL)

NIH$uRDT_result[NIH$uRDT_overall=='Discordant'] <- NA
NIH$uRDT_result[NIH$uRDT_overall=='+'] <- 1
NIH$uRDT_result[NIH$uRDT_overall=='-'] <- 0
NIH$uRDT_result <- as.integer(NIH$uRDT_result)

#Make binary variable for HRP2, Pan LDH, and Pv LDH positivity
NIH <- NIH %>%
  mutate(HRP2_pg_ml_pos = 
           as.numeric(HRP2_pg_ml_pos == "positive")) %>%
  mutate(LDH_Pan_pg_ml_pos = 
           as.numeric(LDH_Pan_pg_ml_pos == "positive")) %>%
  mutate(LDH_Pv_pg_ml_pos = 
           as.numeric(LDH_Pv_pg_ml_pos == "positive"))

#Replace missing values with 0 per 1000 WBC's for microscopy results
#Make binary +/- variable for P. falciparum micro result
#Then create variable for total no. of Pf parasites in p/uL
#using CDC conversion of 8000 WBC's per uL
#note - treat those Po/Pm infected as microscopy neg

NIH <- NIH %>%
  mutate_at(vars(`P. falciparum`, 
                 `Avg. Pf Gam`, 
                 `P. malariae`, 
                 `P. ovale`), ~replace_na(., 0)) %>%
  mutate(pf_parasites = `P. falciparum` + `Avg. Pf Gam`) %>%
  mutate(micro_result = ifelse(pf_parasites != 0, '+', '-')) %>%
  mutate(micro_result_num = as.numeric(micro_result == '+')) %>%
  mutate(parasites_uL = pf_parasites*8)

#Clean up HRP2 (pg/mL) and pLDH (pg/mL) Quansys variables
#Convert values beyond LOQ to LLOQ/2 and ULOQ*2
NIH$HRP2_pg_ml[NIH$HRP2_pg_ml == '< 1.07'] <- 1.07/2
NIH$HRP2_pg_ml[NIH$HRP2_pg_ml == '> 16500.00'] <- 16500.00*2
NIH$HRP2_pg_ml<- as.numeric(NIH$HRP2_pg_ml)

NIH$LDH_Pan_pg_ml[NIH$LDH_Pan_pg_ml == '< 14.41'] <- 14.41/2
NIH$LDH_Pan_pg_ml[NIH$LDH_Pan_pg_ml == '> 525700.00'] <- 525700.00*2
NIH$LDH_Pan_pg_ml<- as.numeric(NIH$LDH_Pan_pg_ml)

hist(log10(NIH$HRP2_pg_ml))
hist(log10(NIH$LDH_Pan_pg_ml))

#convert bleed day (i.e. study day) to numeric variable
NIH$`Bleed Day`[NIH$`Bleed Day` == 'D547'] = 547
NIH$`Bleed Day`[NIH$`Bleed Day` == 'D554'] = 554
NIH$`Bleed Day`[NIH$`Bleed Day` == 'D568'] = 568
NIH$`Bleed Day`[NIH$`Bleed Day` == 'D582'] = 582
NIH$`Bleed Day` <- as.numeric(NIH$`Bleed Day`)

#II. Exploratory Data Analysis

#descriptive analyses

mean(NIH$mean_age, na.rm = T)
sd(NIH$mean_age, na.rm = T)
sum(NIH$micro_result_num == 1, na.rm = T)
sum(NIH$uRDT_result == 1, na.rm = T)

#calculate no. of study participants + by microscopy and uRDT
#at least once during sampling window
NIH %>% group_by(`Volunteer ID`) %>%
  summarise(uRDT_pos = sum(uRDT_result, na.rm = T),
            micro_pos = sum(micro_result_num, na.rm = T)) %>%
  summarise(sum(uRDT_pos > 0), sum(micro_pos > 0))

#Filter out non-Pf infections from further analyses
NIH <- NIH %>%
  mutate(pm_only = (`P. malariae` != 0 & pf_parasites == 0)) %>%
  mutate(po_only = (`P. ovale` != 0 & pf_parasites == 0)) %>%
  filter(pm_only == FALSE & po_only == FALSE)

#uRDT Performance vs. HRP2 Quansys Reference, pLDH Quansys,
#and microsopy result

reportROC(gold=NIH$HRP2_pg_ml_pos,
          predictor=as.numeric(NIH$uRDT_result),
          important="se",plot=TRUE)

reportROC(gold=NIH$LDH_Pan_pg_ml_pos,
          predictor=as.numeric(NIH$uRDT_result),
          important="se",plot=TRUE)

reportROC(gold=NIH$micro_result_num,
          predictor=as.numeric(NIH$uRDT_result),
          important="se",plot=TRUE)

#uRDT performance statistics

# % of infections both microscopy and HRP2+ that uRDT detects
NIH %>% filter(HRP2_pg_ml_pos == 1 & micro_result_num == 1) %>%
  summarise(n = n(), uRDT_pos = sum(uRDT_result == 1, na.rm = T),
            percent_uRDT = uRDT_pos/n*100)

# % of infections HRP2+ only that uRDT detects
NIH %>% filter(HRP2_pg_ml_pos == 1) %>%
  summarise(n = n(), uRDT_pos = sum(uRDT_result == 1, na.rm = T),
            percent_uRDT = uRDT_pos/n*100)

# % of infections pLDH+ that uRDT detects
NIH %>% filter(LDH_Pan_pg_ml_pos == 1) %>%
  summarise(n = n(), uRDT_pos = sum(uRDT_result == 1, na.rm = T),
            percent_uRDT = uRDT_pos/n*100)

#HRP2 and pLDH summary statistics
geoMean(NIH$HRP2_pg_ml, na.rm = T)
geoSD(NIH$HRP2_pg_ml, na.rm = T)

geoMean(NIH$LDH_Pan_pg_ml, na.rm = T)
geoSD(NIH$LDH_Pan_pg_ml, na.rm = T)

#linear model between # Pf parasites (including gametocytes)
#and HRP2, pLDH concentrations (pg/mL)
mod_hrp2 <- lm(log10(HRP2_pg_ml) ~ log10(parasites_uL+1),
               data = NIH)
summary(mod_hrp2)

mod_ldh <- lm(log10(LDH_Pan_pg_ml) ~ log10(parasites_uL+1),
              data = NIH)
summary(mod_ldh)

#HRP2 and pLDH by uRDT, microscopy results
NIH %>% group_by(uRDT_result) %>%
  summarise(geoMean(HRP2_pg_ml, na.rm = T),
            geoMean(LDH_Pan_pg_ml, na.rm = T))

NIH %>% group_by(micro_result_num) %>%
  summarise(geoMean(HRP2_pg_ml, na.rm = T),
            geoMean(LDH_Pan_pg_ml, na.rm = T))

#III. Venn Diagram Visualization (Figure 1)
#visualize overlap of detection by various malaria diagnostic tests
# (Microscopy, uRDT, HRP2 by Quansys, and pLDH by Quansys)
#output as PDF file

library(eulerr)
vd <- NIH %>%
  select(HRP2_pg_ml_pos, LDH_Pan_pg_ml_pos, uRDT_result, micro_result_num) %>%
  drop_na()

n1 = sum(vd$HRP2_pg_ml_pos == 1)
n2 = sum(vd$LDH_Pan_pg_ml_pos == 1)
n3 = sum(vd$uRDT_result == 1)
n4 = sum(vd$micro_result_num == 1)
n12 = sum(vd$HRP2_pg_ml_pos == 1 & vd$LDH_Pan_pg_ml_pos == 1)
n13 = sum(vd$HRP2_pg_ml_pos == 1 & vd$uRDT_result == 1)
n14 = sum(vd$HRP2_pg_ml_pos == 1 & vd$micro_result_num == 1)
n23 = sum(vd$LDH_Pan_pg_ml_pos == 1 & vd$uRDT_result == 1)
n24 = sum(vd$LDH_Pan_pg_ml_pos == 1 & vd$micro_result_num == 1)
n34 = sum(vd$uRDT_result == 1 & vd$micro_result_num == 1)
n123 = sum(vd$HRP2_pg_ml_pos == 1 & vd$LDH_Pan_pg_ml_pos == 1 &
             vd$uRDT_result == 1)
n124 = sum(vd$HRP2_pg_ml_pos == 1 & vd$LDH_Pan_pg_ml_pos == 1 &
             vd$micro_result_num == 1)
n134 = sum(vd$HRP2_pg_ml_pos == 1 & vd$uRDT_result == 1 &
             vd$micro_result_num == 1)
n234 = sum(vd$LDH_Pan_pg_ml_pos == 1 & vd$uRDT_result == 1 &
             vd$micro_result_num == 1)
n1234 = sum(vd$HRP2_pg_ml_pos == 1 & vd$LDH_Pan_pg_ml_pos == 1
                    & vd$uRDT_result == 1 & vd$micro_result_num == 1)

Venn <- euler(c(A = n1, B = n2, C = n3, D= n4 ,
                "A&B" = n12, "A&C" = n13, "A&D"=n14, 
                "B&C" = n23, "B&D"=n24, "C&D"=n34,
                "A&B&C" = n123, "A&B&D"=n124, "A&C&D"=n134, "B&C&D"=n234,
                "A&B&C&D"=n1234),input = c("union"),
              shape = c("ellipse"))

pdf(file = "Mali_Venn.pdf", width = 8.6, height = 6.3)
plot(Venn, legend = list(labels = c("HRP2", "pLDH", "uRDT","Microscopy"), cex = 1.3),
     fill = c("lightblue", "tomato", "violetred3", "slateblue1"),
     labels = c("HRP2", "pLDH", "uRDT","Microscopy"), 
     Count=TRUE, edges = FALSE, fontsize = 4, quantities = list(cex = 1.3))
dev.off()
