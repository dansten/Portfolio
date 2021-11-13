plasma <- read.delim("Data/plasma.txt")
head(plasma)

library(ggplot2)

plasma$sex <- factor(plasma$sex,
                     levels = c(1, 2),
                     labels = c("Male", "Female"))
plasma$smokstat <- factor(plasma$smokstat,
                          levels = c(1, 2, 3),
                          labels = c("Never", "Former", "Current Smoker"))
plasma$bmicat <- factor(plasma$bmicat,
                        levels = c(1, 2, 3, 4),
                        labels = c("Underweight", "Normal", "Overweight", "Obese"))
plasma$vituse <- factor(plasma$vituse,
                        levels = c(1, 2, 3),
                        labels = c("Yes, fairly often", "Yes, not often", "No"))

plasma <- plasma[plasma$betaplasma > 0,]

plasma$bmicat <- relevel(plasma$bmicat, "Normal")
plasma$sex <- relevel(plasma$sex, "Female")
plasma$smokstat <- relevel(plasma$smokstat, "Never")

#Part 0 ####
#0a ####
#0.42 mikromol/l -> ng/ml

# g/mol 
a <- 225.498

#0b ####
plasma$lowplasma <- as.numeric(plasma$betaplasma < a)
plasma$plasmacat <- factor(plasma$lowplasma,
                           levels = c(0, 1),
                           labels = c("high", "low"))
ggplot(plasma, aes(age, lowplasma, color = plasmacat)) +
  geom_point()+
  xlab("Age") +
  ylab("Low plasma beta-carotene") +
  labs(title = "Low PBC (=1) or high PBC (=0) vs age") +
  theme(text = element_text(size = 14))
         
(lowcount <- length(which(plasma$lowplasma == 1)))
#loiw = 235
(highcount <- length(which(plasma$lowplasma == 0))) 
#high = 80

#percentage low
(problow <- lowcount/length(plasma$plasmacat))
#0.7460317

#Part 1 ####
#1a ####
(probtable <- table(plasma$smokstat, plasma$lowplasma))
#probability for low, and each smokstat
(pneverlow <- probtable[1,2]/sum(probtable[1,]))
(pformerlow <- probtable[2,2]/sum(probtable[2,]))
(pcurrlow <- probtable[3,2]/sum(probtable[3,]))
(p <- c(pneverlow, pformerlow, pcurrlow))
#odds ratio = p / 1 + p
(oddsneverlow <- p[1]/(1-p[1]))
(oddsformerlow <- p[2]/(1-p[2]))
(oddscurrlow <- p[3]/(1-p[3]))
(odds <- c(oddsneverlow, oddsformerlow, oddscurrlow)) 


#odds ratios
(odds1 <- oddsneverlow/oddsneverlow)
(odds2 <- oddsformerlow/oddsneverlow)
(odds3 <- oddscurrlow/oddsneverlow)
(oddsratio <- c(odds1, odds2, odds3)) 
(probtable <- cbind(probtable, p, odds, oddsratio))

#logistic model
model.log <- glm(lowplasma ~ smokstat, family="binomial", data = plasma)
summary(model.log)$deviance

plasma$smokstatCurrent <- as.numeric(plasma$smokstat == "Current Smoker")
plasma$smokstatFormer <- as.numeric(plasma$smokstat == "Former")
plasma$smokstatNever <- as.numeric(plasma$smokstat == "Never")
#Deviancereduced - deviancefull
model.log.1 <- glm(lowplasma ~ smokstatFormer , family="binomial", data = plasma)
summary(model.log.1)
model.log.2 <- glm(lowplasma ~ smokstatCurrent, family="binomial", data = plasma)

model.log.3 <- glm(lowplasma ~ smokstatNever , family="binomial", data = plasma)
summary(model.log.3)

(Diff.Former <- summary(model.log.1)$deviance - summary(model.log)$deviance)
(Diff.Current <- summary(model.log.2)$deviance - summary(model.log)$deviance)
(Diff.Never <- summary(model.log.3)$deviance - summary(model.log)$deviance)

#GLobal likklihood ratio test
model.null <- glm(lowplasma ~ 1, family = "binomial", data = plasma)
(D_diff <- summary(model.null)$deviance - summary(model.log)$deviance)
(df_diff <- summary(model.null)$df.null - summary(model.log)$df.residual)

qchisq(1 - 0.05, df_diff)
# or P-value:
pchisq(D_diff, df_diff, lower.tail = FALSE)

# beta: log-odds(ratio) with c.i.:
model.log$coefficients
(ci.beta <- confint(model.log))

# Odds (exp(beta0)) and OR, odds ratio, exp(beta1)
exp(model.log$coefficients)
(ci.or <- exp(ci.beta))
probtable

plasma.pred <- cbind(
  plasma,
  phat = predict(model.log, type ="response"),
  logit = predict(model.log, se.fit = TRUE)
)
head(plasma.pred)

# An unnecessary variable:
plasma.pred$logit.residual.scale <- NULL

# Calculate confidence intervals for the log odds
# standard normal quantile:
(lambda <- qnorm(1 - 0.05/2))
plasma.pred$logit.lwr <- plasma.pred$logit.fit - lambda*plasma.pred$logit.se.fit
plasma.pred$logit.upr <- plasma.pred$logit.fit + lambda*plasma.pred$logit.se.fit
head(plasma.pred)

# transform the log-odds intervals into C.I. for odds
plasma.pred$odds.lwr <- exp(plasma.pred$logit.lwr)
plasma.pred$odds.upr <- exp(plasma.pred$logit.upr)
head(plasma.pred)

# transform the odds intervals into C.I. for p
plasma.pred$p.lwr <- plasma.pred$odds.lwr/(1 + plasma.pred$odds.lwr)
plasma.pred$p.upr <- plasma.pred$odds.upr/(1 + plasma.pred$odds.upr)

#C.I for each smokstat category
(plasma.pred[,c(16,23,24)])

summary(model.log)

#1b ####
model.age <- glm(lowplasma ~ age, family = "binomial", data = plasma)
summary(model.age)

model.age$coefficients
(ci.beta <- confint(model.age))

# Odds (exp(beta0)) and OR, odds ratio, exp(beta1)
exp(model.age$coefficients)
(ci.or <- exp(ci.beta))

plasma.pred.age <- cbind(
  plasma,
  phat = predict(model.age, type ="response"),
  logit = predict(model.age, se.fit = TRUE)
)
head(plasma.pred.age)

# An unnecessary variable:
plasma.pred.age$logit.residual.scale <- NULL

# Calculate confidence intervals for the log odds
# standard normal quantile:
plasma.pred.age$logit.lwr <- plasma.pred.age$logit.fit - lambda*plasma.pred.age$logit.se.fit
plasma.pred.age$logit.upr <- plasma.pred.age$logit.fit + lambda*plasma.pred.age$logit.se.fit
head(plasma.pred.age)

# transform the log-odds intervals into C.I. for odds
plasma.pred.age$odds.lwr <- exp(plasma.pred.age$logit.lwr)
plasma.pred.age$odds.upr <- exp(plasma.pred.age$logit.upr)
head(plasma.pred.age)

# transform the odds intervals into C.I. for p
plasma.pred.age$p.lwr <- plasma.pred.age$odds.lwr/(1 + plasma.pred.age$odds.lwr)
plasma.pred.age$p.upr <- plasma.pred.age$odds.upr/(1 + plasma.pred.age$odds.upr)

(plot.p <- ggplot(plasma.pred.age, aes(age, lowplasma)) +
  geom_point() +
  geom_line(aes(y = phat), color = "red", size = 1) +
  geom_ribbon(aes(ymin = p.lwr, ymax = p.upr), alpha = 0.2) +
  xlab("Age") +
  ylab("Low plasma beta-carotene") +
  labs(title = "Low PBC (=1) or high PBC (=0) vs age",
       caption = "red = fitted line, with 95% confidence interval") +
  theme(text = element_text(size = 14)))

plasma.pred.age

coeff <- model.age$coefficients
coeffCi <- confint(model.age)
tsm <- cbind(
  coeff,
  coeffCi
)
tsm <- tsm[2,]
tsm
ages <- c(30, 31, 70, 71)
tsm*ages
tsm30 <- 30*tsm
tsm31 <- 31*tsm
tsm70 <- 70*tsm
tsm71 <- 71*tsm

tsmAges <- rbind(tsm30, tsm31, tsm70, tsm71)
tsmAgesExp <- exp(tsmAges)
tsmAgesExp

(diff3031 <- tsmAgesExp[1,] - tsmAgesExp[2,])
(diff7071 <- tsmAgesExp[3,] - tsmAgesExp[4,])
#1c ####
plasma.pred.age2 <- cbind(plasma,
                   fit = predict(model.age),
                   v = influence(model.age)$hat)
#saving plot to use later
(plot.v <- ggplot(plasma.pred.age2, aes(age, v)) + 
    geom_point() +
    geom_hline(yintercept = 2*length(model.age$coefficients)/nrow(plasma), 
               color = "red", size = 1) +
    theme(text = element_text(size = 14)))

plasma <- plasma[plasma$betaplasma > 0,]
model.lin <- lm(log(betaplasma) ~ age, data = plasma)
plasma.pred.lin <- cbind(plasma,
                          fit = predict(model.lin),
                          v = influence(model.lin)$hat)
plot.v +
  geom_point(data = plasma.pred.lin, 
             color = "blue") +
  xlab("Age") +
  ylab("v, Leverage") +
  labs(title = "Leverage vs. Age",
       caption = "blue = linear model, black = logistic model")

(v.max <- max(plasma.pred.age2$v))
v.largest <- which(plasma.pred.age2$v == v.max)

#add largest leverage to plot.p
plot.p +
  geom_point(data = plasma.pred.age2[v.largest,], color = "magenta", size = 2) + 
  labs(title = "Low PBC (=1) or high PBC (=0) vs age",
       caption = "red = fitted line, with 95% confidence interval, magenta = largest leverage")
#1d ####
plasma.pred.age2$devres <- influence(model.age)$dev.res
plasma.pred.age2$devstd <- plasma.pred.age2$devres/sqrt(1 - plasma.pred.age2$v)
head(plasma.pred.age2)


(devstd.largest <- which(abs(plasma.pred.age2$devstd) > 2))

ggplot(plasma.pred.age2, aes(age, devstd, color = plasmacat)) +
  geom_point() +
  geom_point(data = plasma.pred.age2[v.largest,], color = "magenta", size = 2) + 
  geom_point(data = plasma.pred.age2[devstd.largest,], color = "blue", size = 2) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", size = 1) +
  geom_hline(yintercept = c(-3, 3), linetype = "dotted", size = 1) +
  labs(title = "Standardized deviance residuals vs age",
       color = "Y",
       caption = "magenta = largest leverage, blue = largest standardized deviance residual") +
  xlab("Age") +
  ylab("Standardized deviance residuals") +
  theme(text = element_text(size = 14))

#1e ####
plasma.pred.age2$Dcook <- cooks.distance(model.age)
head(plasma.pred.age2)

Dcook.largest <- which(plasma.pred.age2$Dcook > 0.02)
ggplot(plasma.pred.age2, aes(age, Dcook, color = plasmacat)) +
  geom_point() +
  geom_point(data = plasma.pred.age2[v.largest,], color = "magenta", size = 2) + 
  geom_point(data = plasma.pred.age2[Dcook.largest,], color = "green", size = 2) +
  geom_point(data = plasma.pred.age2[devstd.largest,], color = "blue", size = 2) +
  geom_hline(yintercept = 4/nrow(plasma), linetype = "dotted",
             size = 1) +
  labs(title = "Cook's distance vs age",
       color = "Y", 
       caption = "magenta = largest leverage, blue = largest standardized deviance residual, green = second largest Cook's D") +
  xlab("Age") +
  ylab("Cook's D") +
  theme(text = element_text(size = 14))

#Part 2 ####
#2a ####
#Null-model
(model.0 <- glm(lowplasma ~ 1, family = "binomial", data = plasma))
summary(model.0)

#tried with all background variables, quetelet led to largest decrease in AIC score, we add it
(model.1 <- glm(lowplasma ~ quetelet, family = "binomial", data = plasma))
summary(model.1)

(model.2 <- glm(lowplasma ~ quetelet + smokstat, family = "binomial", data = plasma))
summary(model.2)

(model.3 <- glm(lowplasma ~ quetelet + smokstat + age, family = "binomial", data = plasma))
summary(model.3)

(model.4 <- glm(lowplasma ~ quetelet + smokstat + age + sex, family = "binomial", data = plasma))
summary(model.4)

(collect.AIC <- data.frame(
  AIC(model.0, model.1, model.2, model.3, model.4)
))

model.background <- model.3

background.estim <- model.background$coefficients
background.estim
background.coef <- cbind(
  background.estim,
  confint(model.background)
)
background.coef
exp(background.coef)

head(model.background)
#divide into three categories
plasma$agecat <- cut(plasma$age, breaks = c(0, 40, 55, 100))

plasma.pred.background <- cbind(
  plasma,
  phat = predict(model.background, type ="response"),
  logit = predict(model.background, se.fit = TRUE),
  fit = predict(model.background),
  v = influence(model.background)$hat,
  DCook = cooks.distance(model.background)
)
head(plasma.pred.background)

plasma.pred.background$devres <- influence(model.background)$dev.res
plasma.pred.background$devstd <- plasma.pred.background$devres/sqrt(1 - plasma.pred.background$v)

# An unnecessary variable:
plasma.pred.background$logit.residual.scale <- NULL

# Calculate confidence intervals for the log odds
# standard normal quantile:
plasma.pred.background$logit.lwr <- plasma.pred.background$logit.fit - lambda*plasma.pred.background$logit.se.fit
plasma.pred.background$logit.upr <- plasma.pred.background$logit.fit + lambda*plasma.pred.background$logit.se.fit
head(plasma.pred.background)

# transform the log-odds intervals into C.I. for odds
plasma.pred.background$odds.lwr <- exp(plasma.pred.background$logit.lwr)
plasma.pred.background$odds.upr <- exp(plasma.pred.background$logit.upr)
head(plasma.pred.background)

# transform the odds intervals into C.I. for p
plasma.pred.background$p.lwr <- plasma.pred.background$odds.lwr/(1 + plasma.pred.background$odds.lwr)
plasma.pred.background$p.upr <- plasma.pred.background$odds.upr/(1 + plasma.pred.background$odds.upr)
head(plasma.pred.background)

ggplot(plasma.pred.background, aes(age,lowplasma)) +
  geom_point() +
  geom_point(aes(y = phat, color = bmicat), size = 1) +
  facet_wrap(~smokstat) +
  labs(title = "Low PBC (=1) or high PBC (=0) vs age, for different smoking categories", 
       caption = "Predicted probabilites in colors representing different BMI categories") +
  xlab("Age") +
  ylab("Low plasma beta-carotene") +
  theme(text = element_text(size = 14))

ggplot(plasma.pred.background, aes(quetelet,lowplasma)) +
  geom_point() +
  geom_point(aes(y = phat, color = agecat), size = 1) +
  facet_wrap(~smokstat) +
  labs(title = "Low PBC (=1) or high PBC (=0) vs quetelet, for different smoking categories", 
       caption = "Predicted probabilites in colors representing different Age categories") +
  xlab("Quetelet") +
  ylab("Low plasma beta-carotene") +
  theme(text = element_text(size = 14))

#2b ####
plasma.pred.backround.lev <- cbind(plasma,
                                   fit = predict(model.background),
                                   v = influence(model.background)$hat)

ggplot(plasma.pred.background, aes(age, v, color = bmicat)) +
  geom_point() +
  facet_wrap(~smokstat) +
  geom_hline(yintercept = 2*length(model.background$coefficients)/nrow(plasma), 
             color = "red", size = 1) +
  theme(text = element_text(size = 14))

ggplot(plasma.pred.background, aes(quetelet, v, color = agecat)) +
  geom_point() +
  facet_wrap(~smokstat) +
  geom_hline(yintercept = 2*length(model.background$coefficients)/nrow(plasma), 
             color = "red", size = 1) +
  theme(text = element_text(size = 14))

ggplot(plasma.pred.background, aes(age, devstd, color = bmicat)) +
  geom_point() +
  facet_wrap(~smokstat) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", size = 1) +
  geom_hline(yintercept = c(-3, 3), linetype = "dotted", size = 1) +
  theme(text = element_text(size = 14))

ggplot(plasma.pred.background, aes(quetelet, devstd, color = agecat)) +
  geom_point() +
  facet_wrap(~smokstat) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", size = 1) +
  geom_hline(yintercept = c(-3, 3), linetype = "dotted", size = 1) +
  theme(text = element_text(size = 14))

ggplot(plasma.pred.background, aes(quetelet, DCook, color = agecat)) +
  geom_point() +
  facet_wrap(~smokstat) +
  geom_hline(yintercept = 4/nrow(plasma), linetype = "dotted",
             size = 1) +
  theme(text = element_text(size = 14))

ggplot(plasma.pred.background, aes(age, DCook, color = bmicat)) +
  geom_point() +
  facet_wrap(~smokstat) +
  geom_hline(yintercept = 4/nrow(plasma), linetype = "dotted",
             size = 1) +
  theme(text = element_text(size = 14))

#2c ####
(model.d.0 <- glm(lowplasma ~ vituse + calories + fat + fiber + alcohol + cholesterol + betadiet, family = "binomial", data = plasma))
summary(model.d.0)

(model.d.1 <- glm(lowplasma ~ vituse + calories + fat + fiber + alcohol + betadiet, family = "binomial", data = plasma))
summary(model.d.1)

(model.d.2 <- glm(lowplasma ~ vituse + calories + fiber + betadiet, family = "binomial", data = plasma))
summary(model.d.2)
plasma$vituseNo <- as.numeric(plasma$vituse == "No")

(model.d.3 <- glm(lowplasma ~ vituseNo + calories + fiber + betadiet, family = "binomial", data = plasma))
summary(model.d.3)

(collect.AIC <- data.frame(
  AIC(model.d.0, model.d.1, model.d.2, model.d.3)
))
model.dietary <- model.d.3

(sum.dietary <- summary(model.dietary))
(beta.dietary <- cbind(model.dietary$coefficients, ci = confint(model.dietary)))
(exp(beta.dietary))

#mcfadden pseduo, unadjusted, adjusted r2
(collect.AIC <- data.frame(
  AIC(model.0, model.age, model.background, model.dietary)
))

collect.AIC$loglik <- 
  c(logLik(model.0)[1],
    logLik(model.age)[1],
    logLik(model.background)[1],
    logLik(model.dietary)[1])

# calculate R2_McF;
collect.AIC$R2McF <- 1 - collect.AIC$loglik/lnL0
# calculate R2_McF,adj. Note that p+1 = df (and df.1):
collect.AIC$R2McF.adj <- 1 - (collect.AIC$loglik - (collect.AIC$df - 1)/2)/lnL0

# Show them as % with one decimal value:
round(100*collect.AIC[, c("R2McF", "R2McF.adj")], digits = 1)

#2d ####
#Start with combined variables from diet and background, remove least significant
(model.f.1 <- glm(lowplasma ~ vituseNo + calories + fiber + betadiet + quetelet + smokstat + age, family = "binomial", data = plasma))
summary(model.f.1)

#remove smokstatFormer
plasma$smokstatCurrent <- as.numeric(plasma$smokstat == "Current Smoker")
(model.f.2 <- glm(lowplasma ~ vituseNo + calories + fiber + betadiet + quetelet + smokstatCurrent + age, family = "binomial", data = plasma))
summary(model.f.2)

#removing of fiber, made it worse, add interaction
(model.f.3 <- glm(lowplasma ~ vituseNo + calories + betadiet + quetelet + smokstatCurrent + age + fiber*(1 + betadiet + quetelet), family = "binomial", data = plasma))
summary(model.f.3)
plasma$fiberBetadiet <- plasma$fiber*plasma$betadiet
plasma$fiberQuetelet <- plasma$fiber*plasma$quetelet

(model.f.4 <- glm(lowplasma ~ vituseNo + calories + smokstatCurrent + age + fiber + fiberBetadiet + fiberQuetelet, family = "binomial", data = plasma))
summary(model.f.4)

(model.f.5 <- glm(lowplasma ~ vituseNo + age + fiber + fiberBetadiet + fiberQuetelet + smokstatCurrent*fiber, family = "binomial", data = plasma))
summary(model.f.5)
plasma$fibersmokstatCurrent <- plasma$fiber*plasma$smokstatCurrent

(model.f.6 <- glm(lowplasma ~ vituseNo + age + fiber + fiberBetadiet + fiberQuetelet + fibersmokstatCurrent, family = "binomial", data = plasma))
summary(model.f.6)

(collect.AIC.2 <- data.frame(
  AIC(model.0, model.dietary, model.background, model.f.1, model.f.2, model.f.3, model.f.4, model.f.5, model.f.6)
))

collect.AIC.2$loglik <- 
  c(logLik(model.0)[1],
    logLik(model.background)[1],
    logLik(model.dietary)[1],
    logLik(model.f.1)[1],
    logLik(model.f.2)[1],
    logLik(model.f.3)[1],
    logLik(model.f.4)[1],
    logLik(model.f.5)[1],
    logLik(model.f.6)[1])

# calculate R2_McF;
collect.AIC.2$R2McF <- 1 - collect.AIC.2$loglik/lnL0
# calculate R2_McF,adj. Note that p+1 = df (and df.1):
collect.AIC.2$R2McF.adj <- 1 - (collect.AIC.2$loglik - (collect.AIC.2$df - 1)/2)/lnL0

# Show them as % with one decimal value:
round(100*collect.AIC.2[, c("R2McF", "R2McF.adj")], digits = 1)

model.final <- model.f.4
(sum.final <- summary(model.final))
(beta.final <- cbind(model.final$coefficients, ci = confint(model.final)))
(exp(beta.final))

#Part 3 ####
#3a ####
#cut-off value
pi <- 0.5

pred.phat <- cbind(
  plasma,
  p.age = predict(model.age, type = "response"),
  p.background = predict(model.background, type = "response"),
  p.diet = predict(model.dietary, type = "response"),
  p.final = predict(model.final, type = "response")
)
head(pred.phat)

pred.phat$yhat.age <- as.numeric(pred.phat$p.age > pi)
pred.phat$yhat.background <- as.numeric(pred.phat$p.background > pi)
pred.phat$yhat.diet <- as.numeric(pred.phat$p.diet > pi)
pred.phat$yhat.final <- as.numeric(pred.phat$p.final > pi)
head(pred.phat)
#age
(row.01 <- table(plasma$lowplasma))

(col.01.age <- table(pred.phat$yhat.age))
(confusion.age <- table(pred.phat$lowplasma, pred.phat$yhat.age))
(spec.age <- confusion.age[1, 1] / row.01[1])
(sens.age <- confusion.age[2, 2] / row.01[2])
(accu.age <- sum(diag(confusion.age)) / sum(confusion.age))
(prec.age <- confusion.age[2, 2] / col.01.age[2])

#background
(col.01.background <- table(pred.phat$yhat.background))
(confusion.background <- table(pred.phat$lowplasma, pred.phat$yhat.background))
(spec.background <- confusion.background[1, 1] / row.01[1])
(sens.background <- confusion.background[2, 2] / row.01[2])
(accu.background <- sum(diag(confusion.background)) / sum(confusion.background))
(prec.background <- confusion.background[2, 2] / col.01.background[2])

#diet
(col.01.diet <- table(pred.phat$yhat.diet))
(confusion.diet <- table(pred.phat$lowplasma, pred.phat$yhat.diet))
(spec.diet <- confusion.diet[1, 1] / row.01[1])
(sens.diet <- confusion.diet[2, 2] / row.01[2])
(accu.diet <- sum(diag(confusion.diet)) / sum(confusion.diet))
(prec.diet <- confusion.diet[2, 2] / col.01.diet[2])

#final
(col.01.final <- table(pred.phat$yhat.final))
(confusion.final <- table(pred.phat$lowplasma, pred.phat$yhat.final))
(spec.final <- confusion.final[1, 1] / row.01[1])
(sens.final <- confusion.final[2, 2] / row.01[2])
(accu.final <- sum(diag(confusion.final)) / sum(confusion.final))
(prec.final <- confusion.final[2, 2] / col.01.final[2])


#3b ####
library(pROC)
library(ResourceSelection)

roc.age <- roc(lowplasma ~ p.age, data = pred.phat)
roc.df.age <- coords(roc.age, transpose = FALSE)
roc.df.age$model <- "age"
roc.background <- roc(lowplasma ~ p.background, data = pred.phat)
roc.df.background <- coords(roc.background, transpose = FALSE)
roc.df.background$model <- "background"
roc.diet <- roc(lowplasma ~ p.diet, data = pred.phat)
roc.df.diet <- coords(roc.diet, transpose = FALSE)
roc.df.diet$model <- "diet"
roc.final <- roc(lowplasma ~ p.final, data = pred.phat)
roc.df.final <- coords(roc.final, transpose = FALSE)
roc.df.final$model <- "final"

roc.df <- rbind(roc.df.age, roc.df.background, roc.df.diet, roc.df.final)

# Plot all the curves, in different colors:
ggplot(roc.df, aes(specificity, sensitivity,
                   color = model)) +
  geom_path(size = 1) +
  coord_fixed() +       # square plotting area
  scale_x_reverse() +   # Reverse scale on the x-axis!
  labs(title = "ROC-curves for all the models") +
  theme(text = element_text(size = 14))

(aucs <- 
    data.frame(
      model = c("Age", "Background", "Diet", "Final"),
      auc = c(auc(roc.age), auc(roc.background), auc(roc.diet),
              auc(roc.final)),
      lwr = c(ci(roc.age)[1], ci(roc.background)[1],
              ci(roc.diet)[1],
              ci(roc.final)[1]),
      upr = c(ci(auc(roc.age))[3], ci(auc(roc.background))[3],
              ci(auc(roc.diet))[3],
              ci(auc(roc.final))[3])))
#compare auc-values
roc.test(roc.age, roc.final)
roc.test(roc.background, roc.final)
roc.test(roc.diet, roc.final)
roc.test(roc.age, roc.background)
roc.test(roc.diet, roc.age)
roc.test(roc.diet, roc.background)

#3c ####
#find optimal cutoff for each model
cutoff.age <- 0.568
roc.df.age[roc.df.age$sensitivity > cutoff.age & 
           roc.df.age$specificity > cutoff.age, ]
(threshold.age <- roc.df.age[roc.df.age$sensitivity > cutoff.age & 
                              roc.df.age$specificity > cutoff.age, ]$threshold[1])

cutoff.background <- 0.6538
roc.df.background[roc.df.background$sensitivity > cutoff.background & 
                    roc.df.background$specificity > cutoff.background, ]
(threshold.background <- roc.df.background[roc.df.background$sensitivity > cutoff.background & 
                                             roc.df.background$specificity > cutoff.background, ]$threshold[1])

cutoff.diet <- 0.6499
roc.df.diet[roc.df.diet$sensitivity > cutoff.diet & 
              roc.df.diet$specificity > cutoff.diet, ]
(threshold.diet <- roc.df.diet[roc.df.diet$sensitivity > cutoff.diet & 
              roc.df.diet$specificity > cutoff.diet, ]$threshold[1])

cutoff.final <- 0.71249
roc.df.final[roc.df.final$sensitivity > cutoff.final & 
               roc.df.final$specificity > cutoff.final, ]
(threshold.final <- roc.df.final[roc.df.final$sensitivity > cutoff.final & 
               roc.df.final$specificity > cutoff.final, ]$threshold[1])

pred.phat$yhat.age <- as.numeric(pred.phat$p.age > threshold.age)
pred.phat$yhat.background <- as.numeric(pred.phat$p.background > threshold.background)
pred.phat$yhat.diet <- as.numeric(pred.phat$p.diet > threshold.diet)
pred.phat$yhat.final <- as.numeric(pred.phat$p.final > threshold.final)
head(pred.phat)
#age
(row.01 <- table(plasma$lowplasma))

(col.01.age <- table(pred.phat$yhat.age))
(confusion.age <- table(pred.phat$lowplasma, pred.phat$yhat.age))
(spec.age <- confusion.age[1, 1] / row.01[1])
(sens.age <- confusion.age[2, 2] / row.01[2])
(accu.age <- sum(diag(confusion.age)) / sum(confusion.age))
(prec.age <- confusion.age[2, 2] / col.01.age[2])

#background
(col.01.background <- table(pred.phat$yhat.background))
(confusion.background <- table(pred.phat$lowplasma, pred.phat$yhat.background))
(spec.background <- confusion.background[1, 1] / row.01[1])
(sens.background <- confusion.background[2, 2] / row.01[2])
(accu.background <- sum(diag(confusion.background)) / sum(confusion.background))
(prec.background <- confusion.background[2, 2] / col.01.background[2])

#diet
(col.01.diet <- table(pred.phat$yhat.diet))
(confusion.diet <- table(pred.phat$lowplasma, pred.phat$yhat.diet))
(spec.diet <- confusion.diet[1, 1] / row.01[1])
(sens.diet <- confusion.diet[2, 2] / row.01[2])
(accu.diet <- sum(diag(confusion.diet)) / sum(confusion.diet))
(prec.diet <- confusion.diet[2, 2] / col.01.diet[2])

#final
(col.01.final <- table(pred.phat$yhat.final))
(confusion.final <- table(pred.phat$lowplasma, pred.phat$yhat.final))
(spec.final <- confusion.final[1, 1] / row.01[1])
(sens.final <- confusion.final[2, 2] / row.01[2])
(accu.final <- sum(diag(confusion.final)) / sum(confusion.final))
(prec.final <- confusion.final[2, 2] / col.01.final[2])

#3d ####
pred.sort.age <- pred.phat[order(pred.phat$p.age), ]
pred.sort.age$rank <- seq(1, nrow(pred.sort.age))
head(pred.sort.age)

pred.sort.background <- pred.phat[order(pred.phat$p.background), ]
pred.sort.background$rank <- seq(1, nrow(pred.sort.background))
head(pred.sort.background)

pred.sort.diet <- pred.phat[order(pred.phat$p.diet), ]
pred.sort.diet$rank <- seq(1, nrow(pred.sort.diet))
head(pred.sort.diet)

pred.sort.final <- pred.phat[order(pred.phat$p.final), ]
pred.sort.final$rank <- seq(1, nrow(pred.sort.final))
head(pred.sort.final)

#age
length(model.age$coefficients)
#g > 2
(HL.age <- hoslem.test(pred.sort.age$lowplasma, pred.sort.age$p.age, g = 4))
#g: 3 - 8
HL.age$expected

# Collect the data in a useful form for plotting:
(HL.df.age <- data.frame(group = seq(1, 4),
                       Obs0 = HL.age$observed[, 1],
                       Obs1 = HL.age$observed[, 2],
                       Exp0 = HL.age$expected[, 1],
                       Exp1 = HL.age$expected[, 2]))

#smallest expected 9.144366

#background
length(model.background$coefficients)
#g > 5
(HL.background <- hoslem.test(pred.sort.background$lowplasma, pred.sort.background$p.background, g = 7))
#g: 7, p + 2 
HL.background$expected

# Collect the data in a useful form for plotting:
(HL.df.background <- data.frame(group = seq(1, 7),
                         Obs0 = HL.background$observed[, 1],
                         Obs1 = HL.background$observed[, 2],
                         Exp0 = HL.background$expected[, 1],
                         Exp1 = HL.background$expected[, 2]))
#smallest expected 2.117457 

#diet
length(model.dietary$coefficients)
#g > 5
(HL.diet <- hoslem.test(pred.sort.diet$lowplasma, pred.sort.diet$p.diet, g = 7))
#g: 7, p + 2 
HL.diet$expected

# Collect the data in a useful form for plotting:
(HL.df.diet <- data.frame(group = seq(1, 7),
                                Obs0 = HL.diet$observed[, 1],
                                Obs1 = HL.diet$observed[, 2],
                                Exp0 = HL.diet$expected[, 1],
                                Exp1 = HL.diet$expected[, 2]))
#smallest expected 2.787815 
#final
length(model.final$coefficients)
#g > 8
(HL.final <- hoslem.test(pred.sort.final$lowplasma, pred.sort.final$p.diet, g = 10))
#g: 7, p + 2 
HL.final$expected

# Collect the data in a useful form for plotting:
(HL.df.final <- data.frame(group = seq(1, 13),
                                Obs0 = HL.final$observed[, 1],
                                Obs1 = HL.final$observed[, 2],
                                Exp0 = HL.final$expected[, 1],
                                Exp1 = HL.final$expected[, 2]))
#smallest expected 1.717411
#age
ggplot(HL.df.age, aes(x = group)) +
  geom_line(aes(y = Obs0, linetype = "observed", color = "Y = 0"), size = 1) +
  geom_line(aes(y = Obs1, linetype = "observed", color = "Y = 1"), size = 1) +
  geom_line(aes(y = Exp0, linetype = "expected", color = "Y = 0"), size = 1) +
  geom_line(aes(y = Exp1, linetype = "expected", color = "Y = 1"), size = 1) +
  labs(title = "Model Age: Observed and expected in each group",
       y = "number of observations") +
  scale_x_continuous(breaks = seq(1, 11)) +
  theme(text = element_text(size = 14))

#background
ggplot(HL.df.background, aes(x = group)) +
  geom_line(aes(y = Obs0, linetype = "observed", color = "Y = 0"), size = 1) +
  geom_line(aes(y = Obs1, linetype = "observed", color = "Y = 1"), size = 1) +
  geom_line(aes(y = Exp0, linetype = "expected", color = "Y = 0"), size = 1) +
  geom_line(aes(y = Exp1, linetype = "expected", color = "Y = 1"), size = 1) +
  labs(title = "Model Background: Observed and expected in each group",
       y = "number of observations") +
  scale_x_continuous(breaks = seq(1, 11)) +
  theme(text = element_text(size = 14))

#diet
ggplot(HL.df.diet, aes(x = group)) +
  geom_line(aes(y = Obs0, linetype = "observed", color = "Y = 0"), size = 1) +
  geom_line(aes(y = Obs1, linetype = "observed", color = "Y = 1"), size = 1) +
  geom_line(aes(y = Exp0, linetype = "expected", color = "Y = 0"), size = 1) +
  geom_line(aes(y = Exp1, linetype = "expected", color = "Y = 1"), size = 1) +
  labs(title = "Model Diet: Observed and expected in each group",
       y = "number of observations") +
  scale_x_continuous(breaks = seq(1, 11)) +
  theme(text = element_text(size = 14))

#final
ggplot(HL.df.final, aes(x = group)) +
  geom_line(aes(y = Obs0, linetype = "observed", color = "Y = 0"), size = 1) +
  geom_line(aes(y = Obs1, linetype = "observed", color = "Y = 1"), size = 1) +
  geom_line(aes(y = Exp0, linetype = "expected", color = "Y = 0"), size = 1) +
  geom_line(aes(y = Exp1, linetype = "expected", color = "Y = 1"), size = 1) +
  labs(title = "Model Final: Observed and expected in each group",
       y = "number of observations") +
  scale_x_continuous(breaks = seq(1, 11)) +
  theme(text = element_text(size = 14))
