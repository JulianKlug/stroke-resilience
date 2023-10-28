library(foreign)
library(lme4)
library(lmerTest) #for p-values in lmer mods
library(emmeans)

## Analyzing global efficiency (delta post attack)

all_modulators_data_path <- '/Users/jk1/temp/resilience_part2/atttacks/mediators/data_attack_Eglob_P1_vs_P2_vs_P3_modulators.sav'
all_mod_df <- read.spss(all_modulators_data_path, to.data.frame = TRUE)

selected_modulators_data_path <- '/Users/jk1/temp/resilience_part2/atttacks/mediators/data_attack_Eglob_P1_vs_P2_vs_P3_modulators_FINAL.sav'
selected_df <- read.spss(selected_modulators_data_path, to.data.frame = TRUE)

# all modulators & all timepoints
print("Fitting model with all modulators and all timepoints")
all_mod_model_fit <- lmer(Eglob ~ TP + age + gender + nihss_hospital + handedness
                            + lesion_size + lesion_side
                            + cgm + sgm + wmd + wmgm
                            + (1|Subject), all_mod_df)
summary(all_mod_model_fit)

# all modulators (only TP1)
print("Fitting model with all modulators and only TP1")
all_mod_model_fit_TP1 <- lmer(Eglob ~ age + gender + nihss_hospital + handedness
                            + lesion_size + lesion_side
                            + cgm + sgm + wmd + wmgm
                            + (1|Subject), all_mod_df, subset = TP == "TP1")
summary(all_mod_model_fit_TP1)

# selected modulators at all timepoints
print("Fitting model with selected modulators and all timepoints")
selected_mod_model_fit <- lmer(Eglob ~ TP + age + nihss_hospital
                            + lesion_size
                            + cgm + sgm + wmgm
                            + (1|Subject), selected_df)

summary(selected_mod_model_fit)

# only TP1
print("Fitting model with selected modulators and only TP1")
selected_mod_model_fit_TP1 <- lmer(Eglob ~ age + nihss_hospital
                            + lesion_size
                            + cgm + sgm + wmgm
                            + (1|Subject), selected_df, subset = TP == "TP1")
summary(selected_mod_model_fit_TP1)


# Using only WMGM as modulator, at all timepoints
wmgm_model_fit <- lmer(Eglob ~ TP + wmgm
                            + (1|Subject), all_mod_df)
summary(wmgm_model_fit)


# Using subset of modulators (size, cgm, wmgm) at TP1
subset_modulators_model_fit <- lmer(Eglob ~ lesion_size
                            + cgm + wmgm
                            + (1|Subject), all_mod_df, subset = TP == "TP1")

summary(subset_modulators_model_fit)

# clinical modulators, TP1 and TP2
clinical_modulators_model_fit <- lmer(Eglob ~ TP + age + nihss_hospital
                            + (1|Subject), all_mod_df, subset = TP == "TP1" | TP == "TP2")
summary(clinical_modulators_model_fit)