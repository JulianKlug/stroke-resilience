library(lme4)
library(lmerTest)
library(sjPlot) #for plotting lmer and glmer mods
library(tidyverse)

# Analysis of mean global efficiency after random and targeted attacks
# 1. Mean was taken over sequential removal of random/targeted nodes
# 2. Here we perform
# - A. a mixed effects model with timepoint as fixed effects and subject and density bin as random effects
# - B. a mixed effects model with time-group interaction as fixed effects and subject and density bin as random effects

df_full <- read_csv("/Users/jk1/stroke_research/resilience_stroke/attack_analysis/mean_glob_eff_after_attack_df.csv")
df_full$timepoint = as.factor(df_full$timepoint) # converting to categorical
df_full$time_group = with(df_full, interaction(timepoint,  group))

# creating subsets for attack types
df_random <- subset(df_full, attack_type == "random")
df_targeted <- subset(df_full, attack_type == "mean_degree_targeted")

# creating subsets for groups
df_random_st <- subset(df_random, group == "st")
df_targeted_st <- subset(df_targeted, group == "st")
df_st <- subset(df_full, group == "st")


# Model 1: testing time point hypothesis for random attacks
random_attack_timepoint_model_fit <- lmer(mean_glob_eff ~ timepoint + (1|subject) + (1|density_bin), df_random_st)
# display results
summary(random_attack_timepoint_model_fit)
plot(random_attack_timepoint_model_fit, main='Analysis of time points after random attack')
res<-resid(random_attack_timepoint_model_fit)
plot(res, main='Analysis of time points after random attack')
sjPlot:: tab_model(random_attack_timepoint_model_fit, title = 'Analysis of time points after random attack')
qqnorm(resid(random_attack_timepoint_model_fit), main = 'QQ plot, Analysis of time points after random attack')

# Model 2: group difference at individual time points hypothesis for random attacks
random_attack_group_model_fit <- lmer(mean_glob_eff ~ time_group + (1|subject) + (1|density_bin), df_random)
# display resultsn
summary(random_attack_group_model_fit)
plot(random_attack_group_model_fit, main='Analysis of group difference at time points after random attack')
res<-resid(random_attack_group_model_fit)
plot(res, main='Analysis of group difference at time groups after random attack')
sjPlot:: tab_model(random_attack_group_model_fit, title = 'Analysis of group difference at time points after random attack')
qqnorm(resid(random_attack_group_model_fit), main = 'QQ plot, Analysis of group difference at time points after random attack')


# Model 1: testing time point hypothesis for targeted attacks
targeted_attack_timepoint_model_fit <- lmer(mean_glob_eff ~ timepoint + (1|subject) + (1|density_bin), df_targeted_st)
# display results
summary(targeted_attack_timepoint_model_fit)
plot(targeted_attack_timepoint_model_fit, main='Analysis of time points after targeted attack')
res<-resid(targeted_attack_timepoint_model_fit)
plot(res, main='Analysis of time points after targeted attack')
sjPlot:: tab_model(targeted_attack_timepoint_model_fit, title = 'Analysis of time points after targeted attack')
qqnorm(resid(targeted_attack_timepoint_model_fit), main = 'QQ plot, Analysis of time points after targeted attack')

# Model 2: group difference at individual time points hypothesis for targeted attacks
targeted_attack_group_model_fit <- lmer(mean_glob_eff ~ time_group + (1|subject) + (1|density_bin), df_targeted)
# display resultsn
summary(targeted_attack_group_model_fit)
plot(targeted_attack_group_model_fit, main='Analysis of group difference at time points after targeted attack')
res<-resid(targeted_attack_group_model_fit)
plot(res, main='Analysis of group difference at time groups after targeted attack')
sjPlot:: tab_model(targeted_attack_group_model_fit, title = 'Analysis of group difference at time points after targeted attack')
qqnorm(resid(targeted_attack_group_model_fit), main = 'QQ plot, Analysis of group difference at time points after targeted attack')



