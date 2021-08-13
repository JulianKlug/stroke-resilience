library(lme4)
library(lmerTest)
# library(readr)
# library(tidyverse) #for all data wrangling
# library(cowplot) #for manuscript ready figures
library(sjPlot) #for plotting lmer and glmer mods
# library(lattice)
# library(sjstats) #use for r2 functions
# library(rstatix)

# read in data from 4-way ANOVA with between-subject and within-subject factors
df_full <- read_csv("longitudinal_analysis/glob_eff_auc_df.csv")
df_full$timepoint = as.factor(df_full$timepoint) # converting to categorical
df_full$time_group = with(df_full, interaction(timepoint,  group))

df_st <- subset(df_full, group == "st")

# get parameter estimates from a linear regression with random effects

# this model might be wrong as measures from hc count into timepoint 
overall_model_fit <- lmer(glob_eff_auc ~ group * timepoint + (1|subject), df_full)


# probably better to separate into two models
timepoint_model_fit <- lmer(glob_eff_auc ~ timepoint + (1|subject), df_st)
group_model_fit <- lmer(glob_eff_auc ~ time_group + (1|subject), df_full)

model_fit <- timepoint_model_fit

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
