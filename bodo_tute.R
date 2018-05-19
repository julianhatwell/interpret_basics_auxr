pitch = c(233,204,242,130,112,142)
sex = c(rep("female",3),rep("male",3))
my.df = data.frame(sex,pitch)
my.df
xmdl = lm(pitch ~ sex, my.df)
summary(xmdl)
mean(my.df[my.df$sex=="female",]$pitch)
mean(my.df[my.df$sex=="male",]$pitch)
coef(xmdl)[1] + coef(xmdl)[2]

age = c(14,23,35,48,52,67)
pitch = c(252,244,240,233,212,204)
my.df = data.frame(age,pitch)
xmdl = lm(pitch ~ age, my.df)
summary(xmdl)

my.df$age.c = my.df$age - mean(my.df$age)
xmdl = lm(pitch ~ age.c, my.df)
summary(xmdl)
plot(fitted(xmdl),residuals(xmdl))
plot(rnorm(100),rnorm(100))
hist(residuals(xmdl))
qqnorm(residuals(xmdl))
dfbeta(xmdl)
library(lme4)
politeness=
  read.csv("http://www.bodowinter.com/tutorial/politeness_data.csv")
which(is.na(politeness)==TRUE)
boxplot(frequency ~ attitude*gender,
        col=c("white","lightgray"),politeness)

politeness.model = lmer(frequency ~ attitude +
                          (1|subject) + (1|scenario), data=politeness)
politeness.model
summary(politeness.model)

politeness.model = lmer(frequency ~ attitude +
                          gender + (1|subject) +
                          (1|scenario), data=politeness)
politeness.model
summary(politeness.model)

# to compare models
politeness.null = lmer(frequency ~ gender +
                         (1|subject) + (1|scenario), data=politeness,
                       REML=FALSE) # note REML
politeness.model = lmer(frequency ~ attitude +
                          gender +
                          (1|subject) + (1|scenario), data=politeness,
                        REML=FALSE)
anova(politeness.null,politeness.model)

coef(politeness.model) # note static slopes

politeness.model = lmer(frequency ~ attitude +
                          gender + (1+attitude|subject) +
                          (1+attitude|scenario),
                        data=politeness,
                        REML=FALSE)

coef(politeness.model) # different slopes for attitude

# new null model
politeness.null = lmer(frequency ~ gender +
                         (1+attitude|subject) + (1+attitude|scenario),
                       data=politeness, REML=FALSE)
anova(politeness.null,politeness.model)
