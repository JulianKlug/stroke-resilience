# demo1  <- read.csv("https://stats.idre.ucla.edu/stat/data/demo1.csv")
# ## Convert variables to factor
# demo1_fact <- within(demo1, {
#   group <- factor(group)
#   time <- factor(time)
#   id <- factor(id)
# })
#
# par(cex = .6)
#
# with(demo1_fact, interaction.plot(time, group, pulse,
#   ylim = c(5, 20), lty= c(1, 12), lwd = 3,
#   ylab = "mean of pulse", xlab = "time", trace.label = "group"))
#
# demo1.aov <- aov(pulse ~ group * time + Error(id), data = demo1)
# summary(demo1.aov)

library(foreign)
data_path <- '/Users/jk1/Downloads/data_attack_Eglob_C1_vs_P1_vs_P2_vs_P3.sav'

df <- read.spss(data_path, to.data.frame = TRUE)

df_fact <- within(df, {
  Group <- factor(Group)
  TP <- factor(TP)
  Subject <- factor(Subject)
})

with(df_fact, interaction.plot(TP, Group, Eglob,
  ylab = "mean of Eglob", xlab = "time", trace.label = "group"))

df.aov <- aov(Eglob ~ Group + Error(Subject), data = df_fact)
summary(df.aov)
