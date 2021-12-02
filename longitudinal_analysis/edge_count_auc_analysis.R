library(lme4)
library(lmerTest)
library(sjPlot) #for plotting lmer and glmer mods
library(readr)

df_path = '/Users/jk1/temp/stroke_resilience/output/analysis/edge_count_auc.csv'
df_full <- read_csv(df_path)
df_full$timepoint = as.factor(df_full$timepoint) # converting to categorical

df_st <- subset(df_full, group == "ST")

# build a linear mixed model to check for an effect of timepoint on edge count
edge_count_model_fit <- lmer(edge_count_auc ~ timepoint + (1|subject) , df_st)

# display results of linear regression
summary(edge_count_model_fit)
plot(edge_count_model_fit, main='Analysis of patient covariates')
res<-resid(edge_count_model_fit)
plot(res, main='Analysis of patient covariates')

sjPlot:: tab_model(edge_count_model_fit, title = 'Analysis of patient covariates')

qqnorm(resid(edge_count_model_fit), main = 'QQ plot, Analysis of patient covariates')
anova(edge_count_model_fit)