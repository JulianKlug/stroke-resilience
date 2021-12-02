library(lme4)
library(sjPlot) #for plotting lmer and glmer mods
library(tidyverse)

# Analysis of mean global efficiency after random and targeted attacks
# 1. Mean was taken over sequential removal of random/targeted nodes
# 2. AUC was taken over selected density thresholds (0.3-1)
# 3. Here we perform
# - A. a mixed effects model with timepoint as fixed effects and subject as random effects
# - B. a mixed effects model with time-group interaction as fixed effects and subject as random effects

df_full <- read_csv("/Users/jk1/stroke_research/resilience_stroke/attack_analysis/mean_glob_eff_auc_after_attack_df_with_patien_characteristics.csv")
# df_full$timepoint = as.factor(df_full$timepoint) # converting to categorical
df_full$lesion_volume = as.factor(df_full$lesion_volume) # converting to categorical
df_full$time_group = with(df_full, interaction(timepoint,  group))

# creating subsets for attack types
df_random <- subset(df_full, attack_type == "random")
df_targeted <- subset(df_full, attack_type == "mean_degree_targeted")

# creating subsets for groups
df_random_st <- subset(df_random, group == "st")
df_targeted_st <- subset(df_targeted, group == "st")
df_st <- subset(df_full, group == "st")

# Model 1: testing time point hypothesis for both attack types
overall_timepoint_model_fit <- lmer(mean_glob_eff_auc ~ timepoint + (1|subject) + (1|attack_type), df_st)
# display results
summary(overall_timepoint_model_fit)
plot(overall_timepoint_model_fit, main='Analysis of time points in all attack types (random/targeted)')
res<-resid(overall_timepoint_model_fit)
plot(res, main='Analysis of time points in all attack types (random/targeted)')
sjPlot:: tab_model(overall_timepoint_model_fit, title = 'Analysis of time points in all attack types (random/targeted)')
qqnorm(resid(overall_timepoint_model_fit), main = 'QQ plot, Analysis of time points in all attack types (random/targeted)')

# Model 2: group difference at individual time points hypothesis for both attack types
overall_group_model_fit <- lmer(mean_glob_eff_auc ~ time_group + (1|subject)  + (1|attack_type), df_full)
# display resultsn
summary(overall_group_model_fit)
plot(overall_group_model_fit, main='Analysis of group difference at time points in all attack types (random/targeted)')
res<-resid(overall_group_model_fit)
plot(res, main='Analysis of group difference at time points')
sjPlot:: tab_model(overall_group_model_fit, title = 'Analysis of group difference at time points in all attack types (random/targeted)')
qqnorm(resid(overall_group_model_fit), main = 'QQ plot, Analysis of group difference at time points in all attack types (random/targeted)')


# Model 1: testing time point hypothesis for random attacks
# random_attack_timepoint_model_fit <- lmer(mean_glob_eff_auc ~ timepoint + (1|subject), df_random_st)
# correct for patient initial lesion volume and side
random_attack_timepoint_model_fit_corrected <- lmer(mean_glob_eff_auc ~ timepoint + (1|lesion_volume) + (1|subject) + (1|lesion_side), df_random_st)
random_attack_timepoint_model_fit <- random_attack_timepoint_model_fit_corrected
# display results
summary(random_attack_timepoint_model_fit)
plot(random_attack_timepoint_model_fit, main='Analysis of time points after random attack')
res<-resid(random_attack_timepoint_model_fit)
plot(res, main='Analysis of time points after random attack')
sjPlot:: tab_model(random_attack_timepoint_model_fit, title = 'Analysis of time points after random attack')
qqnorm(resid(random_attack_timepoint_model_fit), main = 'QQ plot, Analysis of time points after random attack')
anova(random_attack_timepoint_model_fit)

# Model 2: group difference at individual time points hypothesis for random attacks
random_attack_group_model_fit <- lmer(mean_glob_eff_auc ~ time_group + (1|subject), df_random)
# display resultsn
summary(random_attack_group_model_fit)
plot(random_attack_group_model_fit, main='Analysis of group difference at time points after random attack')
res<-resid(random_attack_group_model_fit)
plot(res, main='Analysis of group difference at time groups after random attack')
sjPlot:: tab_model(random_attack_group_model_fit, title = 'Analysis of group difference at time points after random attack')
qqnorm(resid(random_attack_group_model_fit), main = 'QQ plot, Analysis of group difference at time points after random attack')


# Model 1: testing time point hypothesis for targeted attacks
# targeted_attack_timepoint_model_fit <- lmer(mean_glob_eff_auc ~ timepoint + (1|subject), df_targeted_st)
targeted_attack_timepoint_model_fit_corrected <- lmer(mean_glob_eff_auc ~ timepoint + (1|lesion_volume) + (1|subject) + (1|lesion_side), df_targeted_st)
targeted_attack_timepoint_model_fit <- targeted_attack_timepoint_model_fit_corrected
# display results
summary(targeted_attack_timepoint_model_fit)
plot(targeted_attack_timepoint_model_fit, main='Analysis of time points after targeted attack')
res<-resid(targeted_attack_timepoint_model_fit)
plot(res, main='Analysis of time points after targeted attack')
sjPlot:: tab_model(targeted_attack_timepoint_model_fit, title = 'Analysis of time points after targeted attack')
qqnorm(resid(targeted_attack_timepoint_model_fit), main = 'QQ plot, Analysis of time points after targeted attack')
anova(targeted_attack_timepoint_model_fit)

# Model 2: group difference at individual time points hypothesis for targeted attacks
targeted_attack_group_model_fit <- lmer(mean_glob_eff_auc ~ time_group + (1|subject), df_targeted)
# display resultsn
summary(targeted_attack_group_model_fit)
plot(targeted_attack_group_model_fit, main='Analysis of group difference at time points after targeted attack')
res<-resid(targeted_attack_group_model_fit)
plot(res, main='Analysis of group difference at time groups after targeted attack')
sjPlot:: tab_model(targeted_attack_group_model_fit, title = 'Analysis of group difference at time points after targeted attack')
qqnorm(resid(targeted_attack_group_model_fit), main = 'QQ plot, Analysis of group difference at time points after targeted attack')



