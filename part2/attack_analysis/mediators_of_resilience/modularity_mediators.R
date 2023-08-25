library(foreign)
library(lme4)
library(lmerTest) #for p-values in lmer mods
library(emmeans)

## Analyzing modularity (delta post attack)

all_modulars_data_path <- '/Users/jk1/temp/resilience_part2/atttacks/mediators/data_attack_modularity_P1vsP2vsP3_modulators.spv.sav'
all_modulators_df <- read.spss(all_modulars_data_path, to.data.frame = TRUE)

selected_modulators_data_path <- '/Users/jk1/temp/resilience_part2/atttacks/mediators/data_attack_modularity_P1vsP2vsP3_modulators_FINAL.spv.sav'
selected_df <- read.spss(selected_modulators_data_path, to.data.frame = TRUE)

# all modulators & all timepoints
print("Fitting model with all modulators and all timepoints")

all_modulators_model_fit <- lmer(modul ~ TP + age + gender + nihss_hospital + handedness
                            + lesion_size + lesion_side
                            + cgm + sgm + wmd + wmgm
                            + (1|subj), all_modulators_df)
summary(all_modulators_model_fit)

# all modulators (only TP1)
print("Fitting model with all modulators and only TP1")

all_modulators_model_fit_TP1 <- lmer(modul ~ age + gender + nihss_hospital + handedness
                            + lesion_size + lesion_side
                            + cgm + sgm + wmd + wmgm
                            + (1|subj), all_modulators_df, subset = TP == "TP1")
summary(all_modulators_model_fit_TP1)


## Using only selected modulators
# rename columns
colnames(selected_df)[colnames(selected_df) == "VAR00004"] <- "modul"
colnames(selected_df)[colnames(selected_df) == "VAR00002"] <- "TP"
colnames(selected_df)[colnames(selected_df) == "VAR00003"] <- "subj"

# selected modulators at all timepoints
print("Fitting model with selected modulators and all timepoints")

selected_modulators_model_fit <- lmer(modul ~ TP + age + nihss_hospital
                            + lesion_size
                            + cgm + sgm + wmgm
                            + (1|subj), selected_df)
summary(selected_modulators_model_fit)

# only TP1
print("Fitting model with selected modulators and only TP1")

selected_modulators_model_fit_TP1 <- lmer(modul ~ age + nihss_hospital
                            + lesion_size
                            + cgm + sgm + wmgm
                            + (1|subj), selected_df, subset = TP == "TP1")
summary(selected_modulators_model_fit_TP1)


# Using only WMGM as modulator, at all timepoints
wmgm_model_fit <- lmer(modul ~ TP + wmgm
                            + (1|subj), all_modulators_df)
summary(wmgm_model_fit)