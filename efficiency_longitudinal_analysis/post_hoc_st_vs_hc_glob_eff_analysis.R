library(lme4)
library(sjPlot) #for plotting lmer and glmer mods


df_full <- read_csv("/Users/jk1/stroke_research/resilience_stroke/efficiency_longitudinal_analysis/glob_eff_auc_df.csv")
df_full$timepoint = as.factor(df_full$timepoint) # converting to categorical


## 1. Compare hc vs st at time point 0
tp0_df = df_full[df_full[, "timepoint"] == 0,]

# get parameter estimates from a linear regression with random effects
tp0_model_fit <- lm(glob_eff_auc ~ group, tp0_df)

# display results of linear regression
summary(tp0_model_fit)

plot(tp0_model_fit, main='Time point 0')
res<-resid(tp0_model_fit)
plot(res, main='Time point 0')

sjPlot:: tab_model(tp0_model_fit, title = 'Time point 0')

qqnorm(resid(tp0_model_fit), main = 'QQ plot, Time point 0')


## 2. Compare hc vs st at time point 1
tp1_df = df_full[df_full[, "timepoint"] == 1 | df_full[, "group"] == "hc",]

# get parameter estimates from a linear regression with random effects
tp1_model_fit <- lm(glob_eff_auc ~ group, tp1_df)

# display results of linear regression
summary(tp1_model_fit)

plot(tp1_model_fit, main='Time point 1')
res<-resid(tp1_model_fit)
plot(res, main='Time point 1')

sjPlot:: tab_model(tp1_model_fit, title = 'Time point 1')

qqnorm(resid(tp1_model_fit), main = 'QQ plot, Time point 1')


## 1. Compare hc vs st at time point 3
tp2_df = df_full[df_full[, "timepoint"] == 2 | df_full[, "group"] == "hc",]

# get parameter estimates from a linear regression with random effects
tp2_model_fit <- lm(glob_eff_auc ~ group, tp2_df)

# display results of linear regression
summary(tp2_model_fit)

plot(tp2_model_fit, main='Time point 2')
res<-resid(tp2_model_fit)
plot(res, main='Time point 2')

sjPlot:: tab_model(tp2_model_fit, title = 'Time point 2')

qqnorm(resid(tp2_model_fit), main = 'QQ plot, Time point 2')


