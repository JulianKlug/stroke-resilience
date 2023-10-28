library(foreign)
library("readxl")
library(lme4)
library(lmerTest) #for p-values in lmer mods
library(emmeans)

# only eglob
data_path <- '/Users/jk1/Downloads/data_attack_Eglob_C1_vs_P1_vs_P2_vs_P3.sav'
df <- read.spss(data_path, to.data.frame = TRUE)

group_model_fit <- lmer(order_parameter_auc ~ subject_timepoint + (1|subject), df)
summary(group_model_fit)

emmeans(group_model_fit, pairwise ~ TP, adjust="none")



gmf.rg <- ref_grid(group_model_fit)
tp.emm <- emmeans(gmf.rg, "TP")
pairs(tp.emm, adjust = "bonferroni")
