
## 1. data input & look at a glance.
age <- c(46, 20, 52, 30, 57, 25, 28, 36, 22, 43, 57, 33, 22, 63, 40, 48, 28, 49, 52,
         58, 29, 34, 24, 50)
cholesterol <- c(3.5, 1.9, 4.0, 2.6, 4.5, 3.0, 2.9, 3.8, 2.1, 3.8, 4.1, 3.0, 2.5, 4.6,
                 3.2, 4.2, 2.3, 4.0, 4.3, 3.9, 3.3, 3.2, 2.5, 3.3)

## plot to estimated linear relation.
plot(age, cholesterol, xlab = "Patients age with hyperlipoproteinaemia", xlim = c(20, 70),
     ylab = "amount of cholesterol", panel.first = grid(col = "grey", lty = 2))
choldata <- data.frame(age = age,
                       cholesterol = cholesterol)
library(ggplot2)
ggplot(choldata, aes(x= age, y = cholesterol))+
  geom_point() +
  geom_smooth(method = "loess", color = "red")+
  geom_smooth(method = "lm")
### ------------------------------------------------------------
## Pearson correlation 
cor(age, cholesterol)
## testing against 0
cor.test(age, cholesterol)

## 2. Estimate B0 & B1 by least square methods.

fit <- lm(cholesterol ~ age, data = choldata)
fit
####------------------------
fit <- lm(cholesterol ~ age, data = choldata)
fit

## 3. The data & estimated linear curve
summary(fit)
###----------------------
qf(.975, df1 = 1, df2 = 22)
pf(.975, df1 = 1, df2 = 22)
qf(.975, df1 = 1, df2 = 22)-1


# Estimated linear curve
plot(age, cholesterol, xlab = "patients age", xlim = c(20, 70),
     ylab = "amount of cholesterol", panel.first = grid(col = "grey", lty = 2))
abline(fit, lwd = 2, col = "darkblue")

### adding confidence & prediction band with regression line
data <- data.frame(age = seq(20, 70, length = 50))
## confidence interval for estimated regression line
fit.c <- predict(fit, data, interval = "confidence")
## confidence interval for prediction
fit.p <- predict(fit, data, interval = "prediction")
plot(age, cholesterol, xlab = "Patients age",
     ylab = "amount of cholesterol", xlim = c(20, 70), ylim = c(1, 5),
     panel.first = grid(col = "grey", lty = 2))
abline(fit, lwd = 2, col = "darkblue")
matlines(data[["age"]], fit.c[,2:3],
         col = c("darkred", "darkred"), lwd = 2, lty = c(1, 1))
matlines(data[["age"]], fit.p[,2:3],
         col = c("orange", "orange"), lwd = 2, lty = c(1, 1))
legend("bottomright",
       legend = c("LS estimator", "confidence", "prediction"),
       fill = c("darkblue", "darkred", "orange"))



### 4. Diagnostic plots

par(mfrow = c(1, 2))
plot(fit, which = 1, pch = 19, main = "Diagnostic plot 1",
     panel.first = grid(col = "grey", lty = 2))
plot(fit, which = 2, pch = 19, main = "Diagostic plot 2",
     panel.first = grid(col = "grey", lty = 2))

## -----------------------------------------------
plot(fit, which = 3, pch = 19, main = "Diagnostic plot 3",
     panel.first = grid(col = "grey", lty = 2))
plot(fit, which = 4, pch = 19, main = "Diagnostic plot 4",
     panel.first = grid(col = "grey", lty = 2))

## ----------------------------------------------
plot(fit, which = 5, pch = 19, main = "Diagnostic plot 5",
     panel.first = grid(col = "grey", lty = 2))
plot(fit, which = 6, pch = 19, main = "Diagnostic plot 6",
     panel.first = grid(col = "grey", lty = 2))

par(mfcol = c(2,3))

## Cook's distance
round(cooks.distance(fit), 3)
## leverage
round(hatvalues(fit), 3)

### 5. --------------------------
confint(fit)
fit2 <- lm(cholesterol ~ I(age - 20))
confint(fit2)

####2###--------------------------------------

t <- (1:4)/30
h <- c(11.31, 17.25, 22.30, 24.85)
freefall <- data.frame(time = t, height = h)
ggplot(freefall, aes(x= time, y = height))+
  geom_point() +
  geom_smooth(method = "loess", color = "red")+
  geom_smooth(method = "lm")
fit3 <- lm(h ~ t)
fit3
summary(fit3)

### estimate of g

2*coef(fit3)[3]/100

###-------4=--------------------
plot(fit3, which = 1)
plot(fit3, which = 2)
plot(fit3, which = 3)
plot(fit3, which = 4)
plot(fit3, which = 5)
plot(fit3, which = 6)

###-------------------------
fit4 <- lm(h ~t + I(t^2),
           data = freefall)
fit4
summary(fit4)
####----------------------------------------
confint(fit3)

#### -------------------------------------3.-------------
########-------------------------------------------------

# Icudata <- read.xls(file = "diabetes.xls")
