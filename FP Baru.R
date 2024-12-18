setwd("D:/Statistika ITS/SMST 7/Analisis Survival")
getwd()

## Import Library
library(flexsurv)
library(survival)
library(ggsurvfit)
library(tidycmprsk)
library(KMsurv)

data("kidrecurr")
attach(kidrecurr)
patient = c(patient,patient)
time = c(time1,time2)
infect = c(infect1,infect2)
age = c(age,age)
gender = c(gender,gender)
gn = c(gn,gn)
an = c(an,an)
pkd = c(pkd,pkd)

# Create the data frame
df <- data.frame(
  patient = patient,
  time = time,
  infect = infect,
  age = age,
  gender = gender,
  gn = gn,
  an = an,
  pkd = pkd
)
head(df)

# Menambahkan kolom disease
df$disease <- with(df, ifelse(gn == 1, 1,         # gn
                              ifelse(an == 1, 2,       # an
                                     ifelse(pkd == 1, 3,   # pkd
                                            0))))          # no disease
print(df)

#Statistika deskriptif
summary(df)

# Frekuensi untuk variabel kategorik
table(df$gender)
table(df$gn)
table(df$an)
table(df$pkd)

# Distribusi frekuensi dengan proporsi
prop.table(table(df$gender))
prop.table(table(df$gn))
prop.table(table(df$an))
prop.table(table(df$pkd))

# Define survival time
y <- Surv(df$time, df$infect)

#### Kaplan Meier Curve ####
#age
df$age_group <- cut(df$age,
                    breaks = c(-Inf, 40, 50, Inf), 
                    labels = c("<40", "40-50", ">50"))
df$age_group

survdiff(y~df$age_group)

kmfit_age <- survfit(y~df$age_group)
kmfit_age

kmfit_age %>%
  ggsurvfit() +
  labs(
    x = "Time",
    y = "Survival probabilities"
  ) + 
  add_confidence_interval() +
  ggtitle("Kaplan-Meier Survival Curve with Confidence Interval 95%")

#gender
kmfit_gender <- survfit(y~df$gender)
kmfit_gender

survdiff(y~df$gender)

kmfit_gender %>%
  ggsurvfit() +
  labs(
    x = "Time",
    y = "Survival probabilities"
  ) + 
  add_confidence_interval() +
  ggtitle("Kaplan-Meier Survival Curve with Confidence Interval 95%")


#gn
kmfit_gn <- survfit(y~df$gn)
kmfit_gn

survdiff(y~df$gn)

kmfit_gn %>%
  ggsurvfit() +
  labs(
    x = "Time",
    y = "Survival probabilities"
  ) + 
  add_confidence_interval() +
  ggtitle("Kaplan-Meier Survival Curve with Confidence Interval 95%")

#an
kmfit_an <- survfit(y~df$an)
kmfit_an

survdiff(y~df$an)

kmfit_an %>%
  ggsurvfit() +
  labs(
    x = "Time",
    y = "Survival probabilities"
  ) + 
  add_confidence_interval() +
  ggtitle("Kaplan-Meier Survival Curve with Confidence Interval 95%")

#pkd
kmfit_pkd <- survfit(y~df$pkd)
kmfit_pkd

survdiff(y~df$pkd)

kmfit_pkd %>%
  ggsurvfit() +
  labs(
    x = "Time",
    y = "Survival probabilities"
  ) + 
  add_confidence_interval() +
  ggtitle("Kaplan-Meier Survival Curve with Confidence Interval 95%")

survdiff(y~df$disease)

#### Cox Proportional Hazard Model ####
#### Loglikelihood Ratio Test ####
model1 <- coxph(y~df$age)
summary(model1)

model2 <- coxph(y~df$age + df$gender)
summary(model2)

LRT1 <- -2*model1$loglik[2] - (-2*model2$loglik[2]); LRT1

model3 <- coxph(y~df$age + df$gender + df$gn)
summary(model3)

LRT2 <- -2*model2$loglik[2] - (-2*model3$loglik[2]); LRT2 # model 2 terpilih

model4 <- coxph(y~df$age + df$gender + df$an)
summary(model4)

LRT3 <- -2*model2$loglik[2] - (-2*model4$loglik[2]); LRT3 # model 2 terpilih

model5 <- coxph(y~df$age + df$gender + df$pkd)
summary(model5)

LRT4 <- -2*model2$loglik[2] - (-2*model5$loglik[2]); LRT4 # model 5 terpilih

#### Log-Log Plot #####
#age
plot(kmfit_age, fun = function(x) -log(-log(x)),
     xlab ="Time", 
     ylab = "-log-log S", 
     lty = c("solid", "dashed"),
     col = c("black", "red", "blue"),
     main = "Log-Log Plot of Age",
     cex.main = 0.8)
legend("topright", c("<40", "40-50", ">50"),
       lty=c("solid"), 
       col = c("black", "red", "blue"))

#gender
plot(kmfit_gender, fun = function(x) -log(-log(x)),
     xlab ="Time", 
     ylab = "-log-log S", 
     lty = c("solid", "dashed"),
     col = c("blue", "red"),
     main = "Log-Log Plot of Gender",
     cex.main = 0.8)
legend("topright", c("female", "male"),
       lty=c("solid"), 
       col = c("red", "blue"))

#gn
plot(kmfit_gn, fun = function(x) -log(-log(x)),
     xlab ="Time", 
     ylab = "-log-log S", 
     lty = c("solid", "dashed"),
     col = c("blue", "red"),
     main = "Log-Log Plot of Disease Type GN",
     cex.main = 0.8)
legend("topright", c("Yes", "No"),
       lty=c("solid"), 
       col = c("red", "blue"))

#an
plot(kmfit_an, fun = function(x) -log(-log(x)),
     xlab ="Time", 
     ylab = "-log-log S", 
     lty = c("solid", "dashed"),
     col = c("blue", "red"),
     main = "Log-Log Plot of Disease Type AN",
     cex.main = 0.8)
legend("topright", c("Yes", "No"),
       lty=c("solid"), 
       col = c("red", "blue"))

#pkd
plot(kmfit_pkd, fun = function(x) -log(-log(x)),
     xlab ="Time", 
     ylab = "-log-log S", 
     lty = c("solid", "dashed"),
     col = c("blue", "red"),
     main = "Log-Log Plot of Disease Type PKD",
     cex.main = 0.8)
legend("topright", c("Yes", "No"),
       lty=c("solid"), 
       col = c("red", "blue"))

#### Goodness of Fit ####
check_ph1 <- cox.zph(model1, transform = rank);check_ph1$table
check_ph2 <- cox.zph(model2, transform = rank);check_ph2$table
check_ph3 <- cox.zph(model3, transform = rank);check_ph3$table
check_ph4 <- cox.zph(model4, transform = rank);check_ph4$table
check_ph5 <- cox.zph(model5, transform = rank);check_ph5$table # gender tidak signifikan
check_ph5$y

#### Stratified Cox Model ##############
# strata non interaction
model6 <- coxph(y~df$age_group + strata(df$gender) + df$pkd)
summary(model6)

# strata with interaction 
model7 <- coxph(y~df$age_group + strata(df$gender) + df$pkd + strata(df$gender)*df$age_group + strata(df$gender)*df$pkd)
summary(model7)

# uji Stratified vs Non Stratified Model
model5.female <- coxph(y~df$age + df$gender + df$pkd, subset=(gender == 1))
model5.male <- coxph(y~df$age + df$gender + df$pkd, subset=(gender == 0))
LL.nostrat <- model5.female$loglik[2] + model5.male$loglik[2]; LL.nostrat
LL.strat <- model5$loglik[2]; LL.strat
chisq_stat <- -2*(LL.strat - LL.nostrat); chisq_stat
ngroups <- length(unique(df$gender)); ngroups
ncoef <- length(model5$coefficients); ncoef
pvalue <- 1 - pchisq(chisq_stat, df = (ngroups-1)*ncoef); pvalue # model SC lebih baik (model6)

# Uji SC Model No Interaction vs With Interaction
LRT5 <- -2*model6$loglik[2] - (-2*model7$loglik[2]); LRT5
-2*model6$loglik[2]
-2*model7$loglik[2]
ncoef.dif <- length(model7$coefficients) - length(model5$coefficients); ncoef.dif 
pvalue <- 1 - pchisq(LRT5, df = ncoef.dif); pvalue # strata non interaction lebih baik (model6)

# Estimasi Hazard Ratio (exp(coef))
exp(coef(model6))

# Confidence Interval 95% untuk HR
confint(model6)

# P-value untuk setiap variabel
summary(model6)$coefficients[,4]


