library(foreign)
library(lme4)
library(lmerTest) #for p-values in lmer mods
library(emmeans)

data_path <- '/Users/jk1/temp/resilience_part2/atttacks/data_globalEff_mediators_P2.sav'

df <- read.spss(data_path, to.data.frame = TRUE)


model_fit <- lmer(globEff ~ lesion_side
                            + wmd + wmgm
                            + (1|patient), df)
summary(model_fit)
