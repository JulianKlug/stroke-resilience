library(lme4)
library(lmerTest)
library(sjPlot) #for plotting lmer and glmer mods
library(readr)

df_path = '/Users/jk1/stroke_research/resilience_stroke/longitudinal_analysis/mean_degree_auc_df_with_patient_characteristics.csv'
df_full <- read_csv(df_path)
# df_full$timepoint = as.factor(df_full$timepoint) # converting to categorical
# df_full$lesion_volume = as.factor(df_full$lesion_volume) # converting to categorical

df_st <- subset(df_full, group == "st")

# build a linear mixed model to check for an effect of timepoint on mean degree auc
# mean_degree_auc_model_fit <- lmer(mean_degree_auc ~ timepoint + (1|subject) , df_st)
mean_degree_auc_model_fit_corrected <- lmer(mean_degree_auc ~ timepoint + lesion_volume + lesion_side + (1|subject), df_st)

mean_degree_auc_model_fit <- mean_degree_auc_model_fit_corrected

# display results of linear regression
summary(mean_degree_auc_model_fit)
plot(mean_degree_auc_model_fit, main='Analysis of patient covariates')
res<-resid(mean_degree_auc_model_fit)
plot(res, main='Analysis of patient covariates')

sjPlot:: tab_model(mean_degree_auc_model_fit, title = 'Analysis of patient covariates')

qqnorm(resid(mean_degree_auc_model_fit), main = 'QQ plot, Analysis of patient covariates')
anova(mean_degree_auc_model_fit)