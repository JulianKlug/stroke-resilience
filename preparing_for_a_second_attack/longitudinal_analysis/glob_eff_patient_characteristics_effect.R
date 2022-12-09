library(lme4)
library(lmerTest)
library(sjPlot) #for plotting lmer and glmer mods
library(readr)

df_path = '/Users/jk1/stroke_research/resilience_stroke/longitudinal_analysis/glob_eff_auc_df_with_patient_characteristics.csv'

df_full <- read_csv(df_path)
df_full$timepoint = as.factor(df_full$timepoint) # converting to categorical

df_st <- subset(df_full, group == "st")

# build a linear mixed model to check for an effect of lesion volume
patient_covar_model_fit <- lmer(glob_eff_auc ~ lesion_side + lesion_volume + (1|subject) + (1|timepoint), df_st)

# display results of linear regression
summary(patient_covar_model_fit)
plot(patient_covar_model_fit, main='Analysis of patient covariates')
res<-resid(patient_covar_model_fit)
plot(res, main='Analysis of patient covariates')

sjPlot:: tab_model(patient_covar_model_fit, title = 'Analysis of patient covariates')

qqnorm(resid(patient_covar_model_fit), main = 'QQ plot, Analysis of patient covariates')
anova(patient_covar_model_fit)
