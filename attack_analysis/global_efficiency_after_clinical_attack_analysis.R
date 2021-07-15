library(lme4)
library(sjPlot) #for plotting lmer and glmer mods

# Analysis of global efficiency after clinical attacks

df_full <- read_csv("/Users/jk1/stroke_research/resilience_stroke/attack_analysis/glob_eff_after_attack_df.csv")
df_full$timepoint = as.factor(df_full$timepoint) # converting to categorical
df_full$time_group = with(df_full, interaction(timepoint,  group))

df_st <- subset(df_full, group == "st")

# get parameter estimates from a linear regression with random effects
# probably better to separate into two models

# Model 1: testing time point hypothesis
timepoint_model_fit <- lmer(glob_eff ~ timepoint * attack_location + (1|subject) + (1|density_bin), df_st)

# display results of linear regression
summary(timepoint_model_fit)
plot(timepoint_model_fit, main='Analysis of time points & attack location')
res<-resid(timepoint_model_fit)
plot(res, main='Analysis of time points & attack location')

sjPlot:: tab_model(timepoint_model_fit, title = 'Analysis of time points & attack location')

qqnorm(resid(timepoint_model_fit), main = 'QQ plot, Analysis of time points & attack location')


# Model 2: group difference at individual time points hypothesis
group_model_fit <- lmer(glob_eff ~ time_group * attack_location + (1|subject)  + (1|density_bin), df_full)

# display results of linear regression
summary(group_model_fit)

plot(group_model_fit, main='Analysis of group difference at time points & attack locations')
res<-resid(group_model_fit)
plot(res, main='Analysis of group difference at time points')

sjPlot:: tab_model(group_model_fit, title = 'Analysis of group difference at time points & attack locations')

qqnorm(resid(group_model_fit), main = 'QQ plot, Analysis of group difference at time points & attack locations')

