library(lme4)
library(sjPlot) #for plotting lmer and glmer mods
library(lmerTest) #for p-values in lmer mods
library(readxl)
library(foreign)
library(emmeans)

# Analysis of global efficiency after clinical attacks
# get parameter estimates from a linear regression with random effects

# df_full <- read.spss('/Users/jk1/temp/resilience_part2/atttacks/data_globalEff_C1_P2.sav', to.data.frame = TRUE)


df_full <- read_xlsx("/Users/jk1/temp/resilience_part2/atttacks/resilience2_attacks.xlsx")

df_full$TP[df_full$TP == 'TP1'] <- 1
df_full$TP[df_full$TP == 'TP2'] <- 2
df_full$TP[df_full$TP == 'TP3'] <- 3
df_full$TP[df_full$TP == 'TPC1'] <- 1
df_full$TP[df_full$TP == 'TPC2'] <- 2

# # replace value TPC2 in column TP with TPC1
# df_full$TP[df_full$TP == 'TPC2'] <- 'TPC1'

# drop TP2 for controls in column TP
df_full <- df_full[df_full$TP != 'TPC2',]
# only check a single TP for subjects
df_full <- df_full[df_full$TP != 'TP1',]
df_full <- df_full[df_full$TP != 'TP3',]


# t test between both groups
# t.test(df_full$deltaglobalEff[df_full$group == 'patient '], df_full$deltaglobalEff[df_full$group == 'ctrl    '])


# Model 1: group difference at individual time points hypothesis
(group_model_fit <- lmer(attack_delta_eglob ~ group + (1|subject), df_full))
# group_model_fit <- lmer(deltaEglob ~ Group + (1|Subject), df_full)

# display results of linear regression
summary(group_model_fit)


# gmf.rg <- ref_grid(group_model_fit)
# group.emm <- emmeans(gmf.rg, "group")
# pairs(group.emm, reverse = TRUE)
#
# tp.emm <- emmeans(gmf.rg, "TP")
# pairs(tp.emm, reverse = TRUE)



plot(group_model_fit, main='Analysis of group difference at time points & attack locations')
res<-resid(group_model_fit)
plot(res, main='Analysis of group difference at time points')
sjPlot:: tab_model(group_model_fit, title = 'Analysis of group difference at time points & attack locations')
qqnorm(resid(group_model_fit), main = 'QQ plot, Analysis of group difference at time points & attack locations')
