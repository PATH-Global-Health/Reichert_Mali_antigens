---
title: "NIH - Mali TBV Trial Survival Analysis"
---
  
#Survival Curve Analysis of uRDT+ persistence post antimalarial treatment
#Emily Reichert, Aug. 2019

#Note: Run 'mali_tbv_clinical.Rmd' script first to initialize and clean NIH_clinical dataset
#The following script generates a survival curve and median time to uRDT negativity estimate 
#for the NIH Mali TBV data.
  
#load necessary packages
library(dplyr)
library(survival)
library(survminer)

#Only look at individuals treated within past 70 days with NO persistent infection
#i.e. no Pf parasites by microscopy at any point after 1 day posttreatment
survival <- NIH_clinical %>%
  select(time_from_drug, uRDT_overall, parasites_uL) %>%
  filter(time_from_drug > 0 & time_from_drug < 70) %>%
  mutate(recurrence = ifelse((time_from_drug > 1 & parasites_uL != 0), 1,0)) %>%
  filter(recurrence != 1) %>%
  select(time_from_drug, uRDT_overall)

#Visualize survival curve
pdf(file = "Mali_surv.pdf", width = 7, height = 6)
survminer::ggsurvplot(
  fit = survival::survfit(survival::Surv(time_from_drug, uRDT_overall) ~ 1, data = survival), 
  color = 'springgreen3',
  xlab = "Days Posttreatment", 
  ylab = "Probability of Positive uRDT",
  xlim = c(0,70), ggtheme = theme_survminer(text = element_text(size = 12)))
dev.off()

#calculate median time to uRDT negativity
survival::survfit(survival::Surv(time_from_drug, uRDT_overall) ~ 1, data = survival)

