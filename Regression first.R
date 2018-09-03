## ----Tribolium-----------------------------------------------------------
ICUData <- read.csv(file = "ICUData.csv")
#ICUData <- read.csv()
humidity <- c(0, 12, 29.5, 43, 53, 62.5, 75.5, 85, 93)
weight <- c(8.95, 8.14, 6.67, 6.08, 5.90, 5.83, 4.68, 4.20, 3.72)
plot(humidity, weight, xlab = "Percent rel. humidity", xlim = c(0, 100),
     ylab = "weight loss [mg]", panel.first = grid(col = "grey", lty = 2))

## ----Tribulium1----------------------------------------------------------
cor(humidity, weight)
## testing against 0
cor.test(humidity, weight)

## ----Tribulium2----------------------------------------------------------
fit <- lm(weight ~ humidity)
fit
plot(fit, which = 6)
plot(fit, which = 5)

## ----Tribulium3----------------------------------------------------------
summary(fit)

## ----Tribulium4----------------------------------------------------------
plot(humidity, weight, xlab = "Percent rel. humidity", xlim = c(0, 100),
     ylab = "weight loss [mg]", panel.first = grid(col = "grey", lty = 2))
abline(fit, lwd = 2, col = "darkblue")

## ----Tribulium5----------------------------------------------------------
pred.frame <- data.frame(humidity = seq(0, 100, length = 100))
## confidence interval for estimated regression line
fit.c <- predict(fit, pred.frame, interval = "confidence")
## confidence interval for prediction
fit.p <- predict(fit, pred.frame, interval = "prediction")
plot(humidity, weight, xlab = "Percent rel. humidity",
     ylab = "weight loss [mg]", xlim = c(0, 100), ylim = c(2.5, 10),
     panel.first = grid(col = "grey", lty = 2))
abline(fit, lwd = 2, col = "darkblue")
matlines(pred.frame[["humidity"]], fit.c[,2:3],
         col = c("darkred", "darkred"), lwd = 2, lty = c(1, 1))
matlines(pred.frame[["humidity"]], fit.p[,2:3],
         col = c("orange", "orange"), lwd = 2, lty = c(1, 1))
legend("topright",
       legend = c("LS estimator", "confidence", "prediction"),
       fill = c("darkblue", "darkred", "orange"))

## ----Tribulium6, fig.width = 9-------------------------------------------
par(mfrow = c(1, 2))
plot(fit, which = 1, pch = 19, main = "Diagnostic plot 1",
     panel.first = grid(col = "grey", lty = 2))
plot(fit, which = 2, pch = 19, main = "Diagostic plot 2",
     panel.first = grid(col = "grey", lty = 2))

## ----Tribulium7, fig.width = 9-------------------------------------------
par(mfrow = c(1, 2))
plot(fit, which = 3, pch = 19, main = "Diagnostic plot 3",
     panel.first = grid(col = "grey", lty = 2))
plot(fit, which = 4, pch = 19, main = "Diagnostic plot 4",
     panel.first = grid(col = "grey", lty = 2))

## ----Tribulium8, fig.width = 9-------------------------------------------
par(mfrow = c(1, 2))
plot(fit, which = 5, pch = 19, main = "Diagnostic plot 5",
     panel.first = grid(col = "grey", lty = 2))
plot(fit, which = 6, pch = 19, main = "Diagnostic plot 6",
     panel.first = grid(col = "grey", lty = 2))

## ----Tribulium9----------------------------------------------------------
## Cook's distance
round(cooks.distance(fit), 3)
## leverage
round(hatvalues(fit), 3)

## ----shortleaf-----------------------------------------------------------
## read in the dataset
shortleaf <- read.table(file = "shortleaf.txt", header = TRUE)
plot(shortleaf$Diam, shortleaf$Vol, xlab = "diameter [inch]",
     ylab = "volume [cubic feet]",
     panel.first = grid(col = "grey", lty = 2))

## ----shortleaf1----------------------------------------------------------
fit1 <- lm(Vol ~ Diam + I(Diam^2), data = shortleaf)
summary(fit1)

## ----shortleaf2----------------------------------------------------------
plot(shortleaf$Diam, shortleaf$Vol, xlab = "diameter [inch]",
     ylab = "volume [cubic feet]",
     panel.first = grid(col = "grey", lty = 2))
lines(shortleaf$Diam, predict(fit1), col = "darkblue", lwd = 2)

## ----shortleaf3, fig.height = 8.5, fig.width = 8.5-----------------------
par(mfrow = c(2,2))
plot(fit1)

## ----shortleaf4----------------------------------------------------------
fit2 <- lm(log(Vol) ~ I(log(Diam)), data = shortleaf)
summary(fit2)

## ----shortleaf5----------------------------------------------------------
plot(shortleaf$Diam, shortleaf$Vol, xlab = "diameter [inch]",
     ylab = "volume [cubic feet]",
     panel.first = grid(col = "grey", lty = 2))
lines(shortleaf$Diam, predict(fit1), col = "darkblue", lwd = 2)
lines(shortleaf$Diam, exp(predict(fit2)), col = "darkred", lwd = 2)
legend("topleft", fill = c("darkblue", "darkred"),
       legend = c("polynomial", "exponential"))
library(MASS)
boxcox(shortleaf$Vol ~ Diam, data = shortleaf)
boxcox(shortleaf$Vol ~ Diam,
       data = shortleaf,
       lambda = seq(from =0.1, to = 0.6, by = 0.2))

boxcox(Vol ~ Diam, data = shortleaf,
       lambda = seq(from = 0, 1, by = 0.05))
boxcox(Vol ~ Diam, data = shortleaf,
       lambda = seq(from = 0.2, 0.6, by = 0.01))




library(ggplot2)
Tribol <- data.frame(hum = humidity,
                     wt = weight)
ggplot(Tribol, aes(x= hum, y= wt)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  xlab("humidity") + ylab("weight")


## ----shortleaf6, fig.width=8.5, fig.height=8.5---------------------------
par(mfrow = c(2,2))
plot(fit2)

## ----FatherSon, message=FALSE--------------------------------------------
install.packages("UsingR", repos = "http://cran.r-project.org")
library(UsingR)
data(father.son)
plot(sheight ~ fheight, pch = 20, ylab = "Height of son [inch]",
     xlab = "Height of father [inch]", data = father.son)
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(lm(sheight ~ fheight, data = father.son), lwd=2, col = 'darkblue')

## ----RegressionToMean----------------------------------------------------
curve(dnorm(x, mean = 70, sd = 10), from = 35, to = 105, n = 501,
            lwd = 2, xlab = 'x', ylab = 'Density')
abline(h = 0, lwd = 2)
abline(v = 90, lty = 2, lwd = 2)

## ----ANOVA1way---------------
# --------------------------------------------
ICUData <- read.csv(file = "ICUData.csv")
plot(temperature ~ outcome,
     data = ICUData[-398,])

fit1 <- lm(temperature ~ outcome,
           data = ICUData[-398,])

summary(fit1)
anova(fit1)

fit2 <- aov(temperature ~ outcome,
            data = ICUData[-398,])

summary(fit2)
model.tables(fit2, type = "means")
model.tables(fit2, type = "effects")
#------------------------------------- try
fit3 <- aov(temperature ~ surgery,
            data = ICUData[-398,])

fit2 <- lm(temperature ~ outcome + surgery,
            data = ICUData[-398,])
summary(fit2)
anova(fit2)
#===============-------------
fit4 <- lm(temperature ~ outcome * surgery,
           data = ICUData[-398,])
summary(fit2)
anova(fit2)

plot(temperature ~ outcome + surgery,
           data = ICUData[-398,])
summary(fit2)
anova(fit2)
interaction.plot(ICUData$outcome, ICUData$surgery, ICUData$temperature)
interaction.plot(ICUData$outcome, ICUData$surgery, ICUData$temperature)
plot.design(temperature ~ outcome + surgery, data = ICUData[-398,])

fit5 <- aov(temperature ~ surgery + outcome,
           data = ICUData[-398,])
summary(fit5)

fit6 <- aov(temperature ~ outcome + surgery,
            data = ICUData[-398,])
summary(fit6)

library(car)
Anova(fit5)
Anova(fit6)

with(ICUData, tapply(temperatre, list(surgery, outcome, sex),
                     mean))

fit7 <- lm(temperature ~ outcome + surgery,
            data = ICUData[-398,])
summary(fit7)

fit9 <- lm(temperature ~ age + heart.rate + ICUData$SAPS.II + LOS, data = ICUData)
summary(fit9)
library(MASS)
stepAIC(fit9)



library(RColorBrewer)
mycol <- rev(brewer.pal(4, "Reds")[c(1,3,4)])
load("Sepsis.RData")
logPCT <- log(Sepsis$PCT.ng.ml)
plot(logPCT ~ Sepsis$Diagnosis, ylab = "log-PCT",
     xlab = "diagnosis", col = mycol) ## -> Boxplots

## ----ANOVA1way2----------------------------------------------------------
fit1 <- lm(logPCT ~ Diagnosis, data = Sepsis)
summary(fit1)

## ----ANOVA1way3----------------------------------------------------------
anova(fit1)

## ----ANOVA1way4----------------------------------------------------------
fit2 <- aov(logPCT ~ Diagnosis, data = Sepsis)
summary(fit2)

## ----ANOVA1way5----------------------------------------------------------
oneway.test(logPCT ~ Diagnosis, data = Sepsis, var.equal = TRUE)

## ----ANOVA1way6----------------------------------------------------------
oneway.test(logPCT ~ Diagnosis, data = Sepsis)

## ----ANOVA1way7----------------------------------------------------------
model.tables(fit2, type = "means")
model.tables(fit2, type = "effects")

## ----ANOVA2way, fig.width = 9--------------------------------------------
par(mfrow = c(1, 2))
plot(logPCT ~ Diagnosis + Gender, data = Sepsis) ## -> Boxplots

## ----ANOVA2way1----------------------------------------------------------
fit1 <- lm(logPCT ~ Diagnosis + Gender, data = Sepsis)
summary(fit1)
anova(fit1)

## ----ANOVA2way2----------------------------------------------------------
fit2 <- aov(logPCT ~ Diagnosis + Gender, data = Sepsis)
summary(fit2)
model.tables(fit2, type = "means")
model.tables(fit2, type = "effects")

## ----ANOVA2way3----------------------------------------------------------
fit3 <- lm(logPCT ~ Diagnosis + Gender + Diagnosis:Gender, data = Sepsis)
summary(fit3)
anova(fit3)

## ----ANOVA2way4----------------------------------------------------------
fit4 <- lm(logPCT ~ Diagnosis*Gender, data = Sepsis)
summary(fit4)
anova(fit4)

## ----ANOVA2way5----------------------------------------------------------
anova(fit2, fit4)

## ----ANOVA2way6----------------------------------------------------------
fit5 <- aov(logPCT ~ Diagnosis*Gender, data = Sepsis)
summary(fit5)
model.tables(fit5, type = "means")
model.tables(fit5, type = "effects")

## ----ANOVA2way7----------------------------------------------------------
interaction.plot(Sepsis$Diagnosis, Sepsis$Gender, logPCT)

## ----ANOVA2way8----------------------------------------------------------
interaction.plot(Sepsis$Gender, Sepsis$Diagnosis, logPCT)
###------------------------------------------------------------

fit2 <- lm(temperature ~ outcome + surgery,
           data = ICUData[-398,])
summary(fit2)
#### take all variables for test
fit3 <- lm(temperature ~ .,
           data = ICUData[-1,])
summary(fit3)


