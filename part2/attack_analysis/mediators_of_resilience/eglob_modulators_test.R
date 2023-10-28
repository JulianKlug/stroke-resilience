library(foreign)
library(lme4)
library(lmerTest) #for p-values in lmer mods
library(emmeans)

## Analyzing global efficiency (delta post attack)

all_modulators_data_path <- '/Users/jk1/Downloads/data_attack_Eglob_P1_vs_P2_vs_P3_modulators_size_site_TP1.sav'
all_mod_df <- read.spss(all_modulators_data_path, to.data.frame = TRUE)

all_mod_df$lesion_site <- as.factor(all_mod_df$lesion_site)

test_model_fit <- lmer(Eglob ~ lesion_site
                            + (1|Subject), all_mod_df)
summary(test_model_fit)

anova(test_model_fit)


library(nlme)
subject <- factor(all_mod_df$Subject)
m3 <- lme(Eglob ~ lesion_site, random = 1|subject,
          data = all_mod_df)
summary(m3)