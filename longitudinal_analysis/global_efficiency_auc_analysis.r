library(lme4)
library(lmerTest)
library(readr)
# library(tidyverse) #for all data wrangling
# library(cowplot) #for manuscript ready figures
library(sjPlot) #for plotting lmer and glmer mods
library(lattice)
# library(sjstats) #use for r2 functions
# library(rstatix)

df_path = '/Users/jk1/stroke_research/resilience_stroke/longitudinal_analysis/glob_eff_auc_df_with_patient_characteristics.csv'

# read in data from 4-way ANOVA with between-subject and within-subject factors
df_full <- read_csv(df_path)
# df_full$timepoint = as.factor(df_full$timepoint) # converting to categorical
df_full$lesion_volume = as.factor(df_full$lesion_volume) # converting to categorical
df_full$time_group = with(df_full, interaction(timepoint,  group))

df_st <- subset(df_full, group == "st")

# get parameter estimates from a linear regression with random effects

# this model might be wrong as measures from hc count into timepoint 
overall_model_fit <- lmer(glob_eff_auc ~ group * timepoint + (1|subject), df_full)


# probably better to separate into two models
timepoint_model_fit <- lmer(glob_eff_auc ~ timepoint + (1|subject), df_st)
# timepoint_model_fit_corrected <- lmer(glob_eff_auc ~ timepoint + (lesion_volume|subject) + (lesion_volume||lesion_side), df_st)
timepoint_model_fit_corrected <- lmer(glob_eff_auc ~ timepoint + (1|lesion_volume) + (1|subject) + (1|lesion_side), df_st)
group_model_fit <- lmer(glob_eff_auc ~ time_group + (1|subject), df_full)

model_fit <- timepoint_model_fit_corrected

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

