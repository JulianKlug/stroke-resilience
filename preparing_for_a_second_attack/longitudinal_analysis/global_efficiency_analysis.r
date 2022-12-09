library(lme4)
library(lmerTest)
library(sjPlot) #for plotting lmer and glmer mods


df_full <- read_csv("longitudinal_analysis/glob_eff_df.csv")
df_full$timepoint = as.factor(df_full$timepoint) # converting to categorical
df_full$time_group = with(df_full, interaction(timepoint,  group))

df_st <- subset(df_full, group == "st")

# get parameter estimates from a linear regression with random effects
# probably better to separate into two models

# Model 1: testing time point hypothesis
timepoint_model_fit <- lmer(glob_eff ~ timepoint + (1|subject) + (1|density_bin), df_st)

# display results of linear regression
summary(timepoint_model_fit)
plot(timepoint_model_fit, main='Analysis of time points')
res<-resid(timepoint_model_fit)
plot(res, main='Analysis of time points')

sjPlot:: tab_model(timepoint_model_fit, title = 'Analysis of time points')

qqnorm(resid(timepoint_model_fit), main = 'QQ plot, Analysis of time points')

# Model 2: group difference at individual time points hypothesis
group_model_fit <- lmer(glob_eff ~ time_group + (1|subject)  + (1|density_bin), df_full)

# display results of linear regression
summary(group_model_fit)

plot(group_model_fit, main='Analysis of group difference at time points')
res<-resid(group_model_fit)
plot(res, main='Analysis of group difference at time points')

sjPlot:: tab_model(group_model_fit, title = 'Analysis of group difference at time points')

qqnorm(resid(group_model_fit), main = 'QQ plot, Analysis of group difference at time points')

