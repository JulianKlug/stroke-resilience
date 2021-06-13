library(lme4)
library(lmerTest)
library(readr)
library(tidyverse) #for all data wrangling
library(cowplot) #for manuscript ready figures
library(sjPlot) #for plotting lmer and glmer mods
library(lattice)
library(sjstats) #use for r2 functions
library(rstatix)

df_full <- read_csv("/Users/jk1/stroke_research/resilience_stroke/efficiency_longitudinal_analysis/loc_eff_auc_df.csv")
df_full$timepoint = as.factor(df_full$timepoint) # converting to categorical
df_full$roi = as.factor(df_full$roi) # converting to categorical


# get parameter estimates from a linear regression with random effects
model_fit <- lmer(loc_eff_auc ~ group * timepoint + (1|subject) + (1|roi), df_full)
# display results of linear regression
summary(model_fit)

plot(model_fit)
res<-resid(model_fit)
plot(res)

sjPlot::plot_model(model_fit)
sjPlot:: tab_model(model_fit)
qqmath(model_fit, id = 0.05)
anova(model_fit)

qqnorm(resid(model_fit))

