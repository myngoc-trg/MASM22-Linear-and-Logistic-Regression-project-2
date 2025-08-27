## ----------------------------------------------------------------------------------------------------------------
library(tidyverse); theme_set(theme_bw() + theme(text = element_text(size = 18)))
getwd()  # Check current working directory
list.files("Data")  # Check if 'carotene.xlsx' is inside the "Data" folder
install.packages("readxl")  # Only run once if not installed
library(readxl)
library(dplyr)  # Load dplyr for select()
# Lecture 3 - Elasticity - multiple linear regression
library(caret)
library(pROC)
library(car)
library(rstatix)
library(GGally)
df <- read_excel("MASM22/Data/carotene.xlsx")
head(df)


## ----------------------------------------------------------------------------------------------------------------
plasmab_molar_mass <- 40*12.010788+56*1.0079477
alpha <- 0.42*plasmab_molar_mass
# 225.4882
summary(df)
# 3rd Qu.: 230.0  


## ----------------------------------------------------------------------------------------------------------------
#?ifelse
#plasmab$lowplasma_01 <- ifelse(plasmab$betaplasma < alpha, 1, 0)
glimpse(df)

df |> mutate(
lowplasma_01 = as.numeric(betaplasma < alpha),
lowplasma_hl = factor(lowplasma_01,
levels = c(0, 1),
labels = c("high", "low"))) -> df


## ----------------------------------------------------------------------------------------------------------------
glimpse(df)
ggplot(df, aes(age, lowplasma_01, color = lowplasma_hl)) +
  geom_point() +
  #geom_smooth() +
  xlab("") +
  ylab("Low concentration of plasma β-carotene") +
  labs(title = "Low plasma β-carotene(=1) or Not low β-carotene (=0) vs x",
       color = "Beta-carotene") 


## ----------------------------------------------------------------------------------------------------------------
lowplasma <- as.data.frame(
  count(df, lowplasma_01))
lowplasma$percentage <- round(lowplasma$n*100/sum(lowplasma$n),1)
lowplasma


## ----------------------------------------------------------------------------------------------------------------
plasma_vituse <- df |> count(vituse, lowplasma_hl) |>
pivot_wider(id_cols = vituse,
names_from = lowplasma_hl,
values_from = n)
plasma_vituse


## ----------------------------------------------------------------------------------------------------------------
plasma_vituse
# Prob. for having low plasma B
# Correspnding odds
# odds = Pr(success) / Pr(failure)
plasma_vituse <- plasma_vituse |>
  mutate(
    prob_low = low*100/(low + high),
    odds = prob_low / (100-prob_low)
  )
plasma_vituse


## ----------------------------------------------------------------------------------------------------------------
mutate(df,
       vituse = factor(vituse,
                         levels = c(1, 2, 3),
                         labels = c("often", "notoften", "no"))) -> df
# odd ratios = OR = oddsj/ odd ref
glimpse(df)

plasma_vituse <- plasma_vituse |>
  mutate(
   OR = odds / plasma_vituse$odds[1]
  )
plasma_vituse





## ----------------------------------------------------------------------------------------------------------------
vituse_glm <- glm(lowplasma_01 ~ vituse, family = "binomial", data = df)
vituse_glm
vituse_sum <- summary(vituse_glm)
vituse_sum
vituse_sum$coefficients
bhat <- vituse_glm$coefficients
ci.beta <- confint(vituse_glm)
or = exp(bhat)
ci.or <- exp(ci.beta)
cbind('beta'=bhat,ci.beta, `exp(beta)` = or, ci.or) |> round(digits = 2)


## ----------------------------------------------------------------------------------------------------------------
vituse_x0 <- data.frame(vituse = c("often", "notoften", "no"))

vituse_pred <- predict(vituse_glm, vituse_x0, se.fit = TRUE)
vituse_pred
lambda <- qnorm(1 - 0.05/2)

vituse_x0 |> mutate(
  logit = vituse_pred$fit,
  logit.lwr = vituse_pred$fit - lambda*vituse_pred$se.fit,
  logit.upr = vituse_pred$fit + lambda*vituse_pred$se.fit,
  
  "odds" = exp(vituse_pred$fit),
  odds.lwr = exp(logit.lwr),
  odds.upr = exp(logit.upr),
  
  "prob_low" = exp(vituse_pred$fit)/ (1+exp(vituse_pred$fit)),
  p.lwr = odds.lwr/(1 + odds.lwr),
  p.upr = odds.upr/(1 + odds.upr)) -> vituse_x0

vituse_x0


## ----------------------------------------------------------------------------------------------------------------
D_diff <- vituse_sum$null.deviance - vituse_sum$deviance
df_diff <- vituse_sum$df.null - vituse_sum$df.residual
chi2_alpha <- qchisq(p = 1 - 0.05, df = df_diff)
Pvalue <- pchisq(q = D_diff, df = df_diff, lower.tail = FALSE)

cbind(D_diff, df_diff, chi2_alpha, Pvalue)


## ----------------------------------------------------------------------------------------------------------------
highlightcolors <- c("|d|>2" = "red", "v > 0.045" = "pink", "Cook's D>0.1" = "purple")
ggplot(df, aes(bmi, lowplasma_01)) +
  geom_point() +
  geom_smooth() +
  xlab("BMI (kg/m²)") +
  ylab("Low betaplasma") +
  labs(title = "Low betaplasma (=1) or Not low betaplasma (=0) vs BMI") +
  theme(
    plot.title = element_text(size = 10),        # Title size
    axis.title.x = element_text(size = 9),       # X-axis label size
    axis.title.y = element_text(size = 9),       # Y-axis label size
    axis.text.x = element_text(size = 8),        # X-axis tick values
    axis.text.y = element_text(size = 8)         # Y-axis tick values
  )


## ----------------------------------------------------------------------------------------------------------------
bmi_glm <- glm(lowplasma_01 ~ bmi, family = "binomial", data = df)
bmi_glm
bmi_sum <- summary(bmi_glm)
bmi_sum

bmi_sum$coefficients
bhat <- bmi_glm$coefficients
ci.beta <- confint(bmi_glm)
or = exp(bhat)
ci.or <- exp(ci.beta)
cbind('beta'=bhat,ci.beta, `exp(beta)` = or, ci.or) |> round(digits = 2)


## ----------------------------------------------------------------------------------------------------------------
head(df)

df |> mutate(phat = predict(bmi_glm, type = "response")) -> bmi_pred
bmi_pred <- cbind(
  bmi_pred,
  logit = predict(bmi_glm, se.fit = TRUE))
glimpse(bmi_pred)
bmi_pred |> mutate(logit.residual.scale = NULL) -> bmi_pred

lambda <- qnorm(1 - 0.05/2)
bmi_pred |> mutate(
  logit.lwr = logit.fit - lambda*logit.se.fit,
  logit.upr = logit.fit + lambda*logit.se.fit,
  
  odds.lwr = exp(logit.lwr),
  odds.upr = exp(logit.upr),
  
  p.lwr = odds.lwr/(1 + odds.lwr),
  p.upr = odds.upr/(1 + odds.upr)) -> bmi_pred
glimpse(bmi_pred)

ggplot(bmi_pred, aes(bmi, lowplasma_01)) +
  geom_point() +
  geom_smooth(se = FALSE, linetype = "dashed") +
  geom_line(aes(y = phat), color = "red", size = 1) +
  geom_ribbon(aes(ymin = p.lwr, ymax = p.upr), alpha = 0.2) +
  xlab("BMI (kg/m²)") +
  ylab("Low betaplasma") +
  labs(title = "Low betaplasma (=1) or Not low betaplasma (=0) vs BMI",
       caption = "red = fitted line, with 95% confidence interval, 
       blue dashed = moving average") +
  theme(
    plot.title = element_text(size = 10),        # Title size
    axis.title.x = element_text(size = 9),       # X-axis label size
    axis.title.y = element_text(size = 9),       # Y-axis label size
    axis.text.x = element_text(size = 8),        # X-axis tick values
    axis.text.y = element_text(size = 8),         # Y-axis tick values,        # Y-axis tick values
    plot.caption = element_text(size = 8)        # Caption size
  )


## ----------------------------------------------------------------------------------------------------------------
bmi_sum$coefficients[2]
bhat <- bmi_glm$coefficients[2]
bhat
ci.beta1 <- confint(bmi_glm)
ci.beta1 <- ci.beta1[2,]
#ci.beta1
delta_x <- c(+1,-1,-10)
or = exp(bhat*delta_x)
ci.or <- exp(outer(ci.beta1, delta_x, "*"))

#ci.or
cbind('delta_x' = delta_x,'bhat'=bhat, `OR` = or, t(ci.or)) |> round(digits = 4)




## ----------------------------------------------------------------------------------------------------------------
# Wald test
bmi_sum$coefficients

# LR-test
D_diff <- bmi_sum$null.deviance - bmi_sum$deviance
df_diff <- bmi_sum$df.null - bmi_sum$df.residual
chi2_alpha <- qchisq(p = 1 - 0.05, df = df_diff)
Pvalue <- pchisq(q = D_diff, df = df_diff, lower.tail = FALSE)

cbind(D_diff, df_diff, chi2_alpha, Pvalue)


## ----------------------------------------------------------------------------------------------------------------
bmi_infl <- influence(bmi_glm)
glimpse(bmi_infl)
bmi_lm <- lm(lowplasma_01 ~ bmi, df)

bmi_pred <- cbind(bmi_pred,
                   #xbeta = predict(bmi_glm),
                   v_logistic = bmi_infl$hat,
                  v_linear = hatvalues(bmi_lm))
glimpse(bmi_pred)

pplus1_bmi <- length(bmi_glm$coefficients)
n <- nobs(bmi_glm)

ggplot(bmi_pred, aes(x = bmi)) +
  geom_point(aes(y=v_logistic, color = "logistic")) +
  geom_point(aes(y=v_linear, color = "linear")) +
  geom_hline(yintercept = c(2*pplus1_bmi/n), linetype = "dashed") +
  scale_color_manual(values = c("logistic" = "red", "linear" = "blue")) +
  labs(
    title = "Leverage of BMI",
    x = "BMI (kg/m²)", 
    y = "Leverage",
    color = "Model"
  ) +
  theme(
    plot.title = element_text(size = 10),        # Title size
    axis.title.x = element_text(size = 9),       # X-axis label size
    axis.title.y = element_text(size = 9),       # Y-axis label size
    axis.text.x = element_text(size = 8),        # X-axis tick values
    axis.text.y = element_text(size = 8),         # Y-axis tick values,        # Y-axis tick values
    plot.caption = element_text(size = 8),        # Caption size
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )


## ----------------------------------------------------------------------------------------------------------------
bmi_pred |> mutate(devresid = bmi_infl$dev.res,
                    stddevresid = devresid/sqrt(1 - v_logistic)) -> bmi_pred

ggplot(bmi_pred, aes(x = bmi, 
                      y = stddevresid, 
                      color = as.factor(lowplasma_hl))) +
  geom_point() +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), 
             linetype = c("dotted", "dashed", "solid", "dashed", "dotted"),
             linewidth = 1) +
  labs(title = "Standardised deviance residuals vs BMI",
       x = "BMI (kg/m²)", y = "devstd",
       color = "Y") +
  theme(
    plot.title = element_text(size = 10),        # Title size
    axis.title.x = element_text(size = 9),       # X-axis label size
    axis.title.y = element_text(size = 9),       # Y-axis label size
    axis.text.x = element_text(size = 8),        # X-axis tick values
    axis.text.y = element_text(size = 8),         # Y-axis tick values,        # Y-axis tick values
    plot.caption = element_text(size = 8),        # Caption size
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

filter(bmi_pred, abs(stddevresid)>2)

ggplot(bmi_pred, aes(bmi, lowplasma_01)) +
  geom_point() +
  geom_smooth(se = FALSE, linetype = "dashed") +
  geom_line(aes(y = phat), color = "red", size = 1) +
  geom_point(data = filter(bmi_pred, abs(stddevresid)>2), 
             aes(color = "|d|>2"), size = 3) +
  geom_ribbon(aes(ymin = p.lwr, ymax = p.upr), alpha = 0.2) +
  xlab("BMI (kg/m²)") +
  ylab("Low betaplasma") +
  labs(title = "Low betaplasma (=1) or Not low betaplasma (=0) vs BMI",
       caption = "red = fitted line, with 95% confidence interval, 
       blue dashed = moving average") + 
  scale_color_manual(values = highlightcolors)+
  theme(legend.position = "bottom",
    plot.title = element_text(size = 10),        # Title size
    axis.title.x = element_text(size = 9),       # X-axis label size
    axis.title.y = element_text(size = 9),       # Y-axis label size
    axis.text.x = element_text(size = 8),        # X-axis tick values
    axis.text.y = element_text(size = 8),         # Y-axis tick values,        # Y-axis tick values
    plot.caption = element_text(size = 8),        # Caption size
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )


## ----------------------------------------------------------------------------------------------------------------
bmi_pred <- mutate(bmi_pred, 
                    Dcook = cooks.distance(bmi_glm))
bmi_pred
ggplot(bmi_pred, aes(x = bmi, y = Dcook, 
                      color = as.factor(lowplasma_hl))) +
  geom_point() +
  #geom_point(data = filter(bmi_pred, v_logistic > 0.045), aes(color = "v > 0.045"),
   #          size = 3) +
  #geom_point(data = filter(bmi_pred, Dcook > 0.1), aes(color = "Cook's D>0.1"),
   #          size = 3) +
  #geom_point(data = filter(bmi_pred, abs(stddevresid)>2), 
   #          aes(color = "|d|>2"), size = 3) +
  geom_hline(yintercept = 4/n, linewidth = 1, linetype = "dashed") +
  labs(title = "Cook's distance vs BMI",
       x = "BMI (kg/m²)", 
       color = "Highlight") +
  theme(legend.position = "bottom",
    plot.title = element_text(size = 10),        # Title size
    axis.title.x = element_text(size = 9),       # X-axis label size
    axis.title.y = element_text(size = 9),       # Y-axis label size
    axis.text.x = element_text(size = 8),        # X-axis tick values
    axis.text.y = element_text(size = 8),         # Y-axis tick values,        # Y-axis tick values
    plot.caption = element_text(size = 8),        # Caption size
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )


## ----------------------------------------------------------------------------------------------------------------
arrange(filter(bmi_pred, abs(stddevresid)>2), Dcook)
# largest D cook = 0.17647045, high betaplasma, bmi=45.9	
# Second largest 0.08332085, bmi=39.5


## ----echo=FALSE, warning=FALSE, paged.print=FALSE----------------------------------------------------------------
mutate(df,
       sex = factor(sex,
                         levels = c(1, 2),
                         labels = c("male", "female"))) -> df
mutate(df,
       smokstat = factor(smokstat,
                     levels = c(1, 2, 3),
                     labels = c("never", "former", "current"))) -> df

glimpse(df)
df <- mutate(df, smokstat = relevel(smokstat, "never"), 
             sex = relevel(sex, "female"))
all_glm <- glm(lowplasma_01 ~ bmi + age + calories + fat + cholesterol + fiber + alcohol + betadiet + smokstat + sex + vituse, family = "binomial", data = df)
vif(all_glm)
summary(all_glm)


## ----------------------------------------------------------------------------------------------------------------
excl_glm <- glm(lowplasma_01 ~ (bmi + age + fat + cholesterol + fiber + alcohol + betadiet + smokstat + sex + vituse)^2, family = "binomial", data = df)
excl_glm


## ----warning=FALSE, include=FALSE, paged.print=FALSE-------------------------------------------------------------
excl_sum <- summary(excl_glm)
#excl_sum
null_glm <- glm(lowplasma_01 ~ 1, family = "binomial", data = df)
backward_model <- step(excl_glm, 
     scope = list(lower = null_glm, upper = excl_glm),
     direction = "backward",
     k = log(nobs(excl_glm)))
#backward_model <- step(excl_glm, direction = "backward", k = log(nobs(excl_glm)))
# AIC unless we supply k=log(n)



## ----message=FALSE, warning=FALSE, include=FALSE-----------------------------------------------------------------
# BIC since we supplied k=log(n)

backward_stepwise_model <- step(backward_model, 
     scope = list(lower = null_glm, upper = excl_glm),
     direction = "both",
     k = log(nobs(excl_glm)))


## ----include=FALSE-----------------------------------------------------------------------------------------------
forward_model <- step(null_glm,
     k = log(nobs(excl_glm)),
     scope = list(lower = null_glm, upper = excl_glm),, 
     direction = "forward")


## ----message=FALSE, warning=FALSE, include=FALSE-----------------------------------------------------------------

forward_stepwise_model <- step(forward_model, 
     scope = list(lower = null_glm, upper = excl_glm),
     direction = "both",
     k = log(nobs(excl_glm)))


## ----------------------------------------------------------------------------------------------------------------
forward_stepwise_model
BIC(excl_glm, backward_model, backward_stepwise_model, forward_model, forward_stepwise_model)


## ----------------------------------------------------------------------------------------------------------------
backward_stepwise_sum <- summary(backward_stepwise_model)
backward_stepwise_sum

#glimpse(df)

forward_stepwise_sum <- summary(forward_stepwise_model)
forward_stepwise_sum

forward_stepwise_sum$coefficients
bhat <- forward_stepwise_model$coefficients
ci.beta <- confint(forward_stepwise_model)
or = exp(bhat)
ci.or <- exp(ci.beta)

cbind('beta'=bhat,ci.beta) |> round(digits = 4)


## ----------------------------------------------------------------------------------------------------------------
glimpse(df)
mutate(df, 
       vituse_new = factor(vituse == "no",
                         levels = c(FALSE, TRUE),
                         labels = c("yes", "no"))) -> df
glimpse(df)
count(df, vituse_new)
# not vituseoften = often + notoften as reference
df <- mutate(df, vituse_new = relevel(vituse_new, "yes"))
forward_stepwise_new <- glm(formula = lowplasma_01 ~ betadiet + bmi + vituse_new + age + 
    betadiet:bmi, family = "binomial", data = df)
forward_stepwise_new_sum <- summary(forward_stepwise_new)
forward_stepwise_new_sum


## ----------------------------------------------------------------------------------------------------------------
D_diff <- forward_stepwise_new_sum$deviance - forward_stepwise_sum$deviance
df_diff <- forward_stepwise_new_sum$df.residual - forward_stepwise_sum$df.residual
cbind(D_diff, df_diff)

chi2_alpha <- qchisq(1 - 0.05, df_diff)
Pvalue <- pchisq(D_diff, df_diff, lower.tail = FALSE)
cbind(D_diff, df_diff, chi2_alpha, Pvalue)
forward_stepwise_new_sum$deviance


## ----------------------------------------------------------------------------------------------------------------
anova_stepwise <- anova(forward_stepwise_new, forward_stepwise_model)
anova_stepwise


## ----------------------------------------------------------------------------------------------------------------

df1 <- data.frame(variable = names(null_glm$coefficients),
b_model1 = null_glm$coefficients, row.names = NULL)


df2 <- data.frame(variable = names(excl_glm$coefficients),
b_full = excl_glm$coefficients, row.names = NULL)

df3 <- data.frame(variable = names(backward_model$coefficients),
b_backward = backward_model$coefficients, row.names = NULL)

df4 <- data.frame(variable = names(backward_stepwise_model$coefficients),
b_backstep = backward_stepwise_model$coefficients, row.names = NULL)

df5 <- data.frame(variable = names(forward_model$coefficients),
b_forward = forward_model$coefficients, row.names = NULL)

df6 <- data.frame(variable = names(forward_stepwise_model$coefficients),
b_model6 = forward_stepwise_model$coefficients, row.names = NULL)

df7 <- data.frame(variable = names(forward_stepwise_new$coefficients),
b_redmodel = forward_stepwise_new$coefficients, row.names = NULL)

models <- full_join(df3,df4) |> full_join(df5) |> full_join(df7)
print(models)


## ----------------------------------------------------------------------------------------------------------------
aic <- AIC(backward_model, backward_stepwise_model, forward_model, forward_stepwise_new)
bic <- BIC(backward_model, backward_stepwise_model, forward_model, forward_stepwise_new)
collect.AICetc <- data.frame(aic, bic)
collect.AICetc


## ----------------------------------------------------------------------------------------------------------------

collect.AICetc |> mutate(df.1 = NULL) -> collect.AICetc
collect.AICetc

logLik(null_glm)
lnL0 <- logLik(null_glm)[1]
lnL0

collect.AICetc |> mutate(
  loglik =  c(logLik(backward_model)[1],
              logLik(backward_stepwise_model)[1],
              logLik(forward_model)[1],
              logLik(forward_stepwise_new)[1])) -> collect.AICetc
collect.AICetc

collect.AICetc |> mutate(
  #R2McF = 1 - loglik/lnL0,
  R2McF.adj = 1 - (loglik - (df - 1)/2)/lnL0) -> collect.AICetc
collect.AICetc




## ----------------------------------------------------------------------------------------------------------------
bestmodelaic_infl <- influence(backward_stepwise_model)
glimpse(bestmodelaic_infl)

bestmodelaic_pred <- cbind(df,
                   xbeta = predict(backward_stepwise_model),
                   v = bestmodelaic_infl$hat)
glimpse(bestmodelaic_pred)

bestmodelaic_pred |> mutate(devresid = bestmodelaic_infl$dev.res,
                    stddevresid = devresid/sqrt(1 - v)) -> bestmodelaic_pred
glimpse(bestmodelaic_pred)

pplus1_bestmodelaic <- length(backward_stepwise_model$coefficients)
n <- nobs(backward_stepwise_model)

bestmodelaic_pred |> mutate(pearson = bestmodelaic_infl$pear.res,
                    stdpearson = pearson/sqrt(1 - v)) -> bestmodelaic_pred
glimpse(bestmodelaic_pred)

ggplot(bestmodelaic_pred, aes(x = xbeta, 
                      y = stddevresid, 
                      color = as.factor(lowplasma_hl))) +
  geom_point() +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), 
             linetype = c("dotted", "dashed", "solid", "dashed", "dotted"),
             linewidth = 1) +
  labs(title = "Standardised deviance residuals vs linear predictor (Best AIC model)",
       x = "xb", y = "devstd",
       color = "Y") +
  theme(
    plot.title = element_text(size = 10),        
    axis.title.x = element_text(size = 9),       
    axis.title.y = element_text(size = 9),       
    axis.text.x = element_text(size = 8),        
    axis.text.y = element_text(size = 8),        
    legend.text = element_text(size = 8),        # Legend label size
    legend.title = element_text(size = 9)        # Legend title size
  )

## ----warning=FALSE-----------------------------------------------------------------------------------------------
bestmodelbic_infl <- influence(forward_stepwise_new)
glimpse(bestmodelbic_infl)

bestmodelbic_pred <- cbind(df,
                   xbeta = predict(forward_stepwise_new),
                   v = bestmodelbic_infl$hat)
glimpse(bestmodelbic_pred)

bestmodelbic_pred |> mutate(devresid = bestmodelbic_infl$dev.res,
                    stddevresid = devresid/sqrt(1 - v)) -> bestmodelbic_pred
glimpse(bestmodelbic_pred)

pplus1_bestmodelbic <- length(forward_stepwise_new$coefficients)
n <- nobs(forward_stepwise_new)

bestmodelbic_pred |> mutate(pearson = bestmodelbic_infl$pear.res,
                    stdpearson = pearson/sqrt(1 - v)) -> bestmodelbic_pred
glimpse(bestmodelbic_pred)

ggplot(bestmodelbic_pred, aes(x = xbeta, 
                      y = stddevresid, 
                      color = as.factor(lowplasma_hl))) +
  geom_point() +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), 
             linetype = c("dotted", "dashed", "solid", "dashed", "dotted"),
             linewidth = 1) +
  labs(title = "Standardised deviance residuals vs linear predictor (Best BIC model)",
       x = "xb", y = "devstd",
       color = "Y") +
  theme(
    plot.title = element_text(size = 10),        
    axis.title.x = element_text(size = 9),       
    axis.title.y = element_text(size = 9),       
    axis.text.x = element_text(size = 8),        
    axis.text.y = element_text(size = 8),        
    legend.text = element_text(size = 8),        # Legend label size
    legend.title = element_text(size = 9)        # Legend title size
  )
sum(abs(bestmodelaic_pred$stddevresid) > 2)

sum(abs(bestmodelbic_pred$stddevresid) > 2)


## ----------------------------------------------------------------------------------------------------------------
ggplot(bestmodelaic_pred, aes(sample = stdpearson)) +
  geom_qq() + geom_qq_line() +
  labs(title = "Q-Q-plot standardised residuals (best AIC model)",
       x = "theoretical", y = "sample")  +
  theme(
    plot.title = element_text(size = 10),        
    axis.title.x = element_text(size = 9),       
    axis.title.y = element_text(size = 9),       
    axis.text.x = element_text(size = 8),        
    axis.text.y = element_text(size = 8),        
    legend.text = element_text(size = 8),        # Legend label size
    legend.title = element_text(size = 9)        # Legend title size
  )

ggplot(bestmodelbic_pred, aes(sample = stdpearson)) +
  geom_qq() + geom_qq_line() +
  labs(title = "Q-Q-plot standardised residuals (best BIC model)",
       x = "theoretical", y = "sample") +
  theme(
    plot.title = element_text(size = 10),        
    axis.title.x = element_text(size = 9),       
    axis.title.y = element_text(size = 9),       
    axis.text.x = element_text(size = 8),        
    axis.text.y = element_text(size = 8),        
    legend.text = element_text(size = 8),        # Legend label size
    legend.title = element_text(size = 9)        # Legend title size
  )

sum(abs(bestmodelaic_pred$stdpearson) > 2)

sum(abs(bestmodelbic_pred$stdpearson) > 2)


## ----------------------------------------------------------------------------------------------------------------
df |> mutate(
  p_AIC = predict(backward_stepwise_model, type = "response"),
  p_BIC = predict(forward_stepwise_new, type = "response")
) -> pred_phat
glimpse(pred_phat)

pred_phat |> mutate(
  yhat_AIC = factor(p_AIC > 0.5,
                     levels = c(FALSE, TRUE),
                     labels = c("high", "low"))) -> pred_phat

pred_phat |> mutate(
  yhat_BIC = factor(p_BIC > 0.5,
                     levels = c(FALSE, TRUE),
                     labels = c("high", "low"))) -> pred_phat
#view(pred_phat)
cm_AIC <- confusionMatrix(
  data = pred_phat$yhat_AIC, 
  reference = pred_phat$lowplasma_hl,
  positive = "low")
cm_AIC

cm_BIC <- confusionMatrix(
  data = pred_phat$yhat_BIC, 
  reference = pred_phat$lowplasma_hl,
  positive = "low")
cm_BIC


## ----------------------------------------------------------------------------------------------------------------
roc_AIC <- roc(lowplasma_01 ~ p_AIC, data = pred_phat)
roc_AIC

coords(roc_AIC) |> head()

roc_BIC <- roc(lowplasma_01 ~ p_BIC, data = pred_phat)
roc_BIC

coords(roc_BIC) |> head()

ggroc(list(AIC = roc_AIC, BIC = roc_BIC)) +
  coord_fixed() +
  labs(title = "ROC-curves for model AIC and the BIC model") +
  theme(
    plot.title = element_text(size = 10),        
    axis.title.x = element_text(size = 9),       
    axis.title.y = element_text(size = 9),       
    axis.text.x = element_text(size = 8),        
    axis.text.y = element_text(size = 8),        
    legend.text = element_text(size = 8),        # Legend label size
    legend.title = element_text(size = 9)        # Legend title size
  )

## ----------------------------------------------------------------------------------------------------------------
aucs <- 
  data.frame(
    model = c("AIC", "BIC"),
    auc = c(auc(roc_AIC), auc(roc_BIC)),
    lwr = c(ci(roc_AIC)[1], ci(roc_BIC)[1]),
    upr = c(ci(auc(roc_AIC))[3], ci(auc(roc_BIC))[3]))
aucs

## ----------------------------------------------------------------------------------------------------------------
roc.test(roc_AIC, roc_BIC)


## ----------------------------------------------------------------------------------------------------------------
topleft_AIC <- coords(roc_AIC, "best", best.method = "closest.topleft")
topleft_AIC
topleft_AIC$threshold
topleft_BIC <- coords(roc_BIC, "best", best.method = "closest.topleft")
topleft_BIC


## ----------------------------------------------------------------------------------------------------------------


pred_phat |> mutate(
  yhat_new_AIC = factor(p_AIC > topleft_AIC$threshold,
                     levels = c(FALSE, TRUE),
                     labels = c("high", "low"))) -> pred_phat

pred_phat |> mutate(
  yhat_new_BIC = factor(p_BIC > topleft_BIC$threshold,
                     levels = c(FALSE, TRUE),
                     labels = c("high", "low"))) -> pred_phat

#view(pred_phat)
cm_new_AIC <- confusionMatrix(
  data = pred_phat$yhat_new_AIC, 
  reference = pred_phat$lowplasma_hl,
  positive = "low")
cm_new_AIC

cm_new_BIC <- confusionMatrix(
  data = pred_phat$yhat_new_BIC, 
  reference = pred_phat$lowplasma_hl,
  positive = "low")
cm_new_BIC


## ----------------------------------------------------------------------------------------------------------------
exp(best_BIC$coefficients)
exp(confint(best_BIC))

