data <- read.csv("Seminole.Nutri.BO.phenotypes.V2.csv",header=T,na.string = "-999")

library(Anova)
lm.hcw <- lm(hcw ~ as.factor(CG)  ,data=na.omit(data))
summary(lm.hcw)
anova(lm.hcw)

library(car)
Anova(lm.hcw, type=3)