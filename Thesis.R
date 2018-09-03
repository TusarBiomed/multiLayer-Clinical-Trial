
# Assume a clinical trial where several dose need to be compare.
# In case of type 2 diabetes trial we consider two endpoints.
# Those are 1. Endpoint X (Hemoglobin A1c or HbA1c) & 
# 2. Endpoint Y (Fasting serum glucose)
# We can use package "mvtnorm" for multi dimensional comparison
## Here we have to compare three different dose with placebo
# We labeled three dose as a, b & c and placebo as d.
# dose a = .25 mg, b = .50 mg & c = .80 mg
# We will look which dose gives statistically significant result in both endpoint.
library(mvtnorm)

# Here we can directly upload all primary parameters instead of doing it again & again.
load(file = "ParametersData.RData")
# number of simulations
M <- 1000

# sample size in each of four arm
n <- 80

# type I error
alpha <- 0.05

## Parameters of dose a 
# Standardized effect size(dosage) in two consecutive endpoint
mu1 <- c(.55, .63)
# Standard deviation sequentially 
SD1 <- c(.8, 1)
# for simplification we consider there are no any
# correlation between endpoints
rho1 <- 0
# Covariance matrix
S1 <- diag(SD1)
S1[1,2] <- rho1
S1[2,1] <- rho1


## Parameters of dose b
# Standardized effect size(dosage) in two consecutive endpoint
mu2 <- c(0.68, 0.78)
# Standard deviation sequentially 
SD2 <- c(.87, 1)
# for simplification we consider there are no any
# correlation between endpoints
rho2 <- 0
# Covariance matrix
S2 <- diag(SD2)
S2[1,2] <- rho2
S2[2,1] <- rho2


## Parameters of dose c
# Standardized effect size(dosage) in two consecutive endpoint
mu3 <- c(0.73, 0.81)
# Standard deviation sequentially 
SD3 <- c(1, .9)
# for simplification we consider there are no any
# correlation between endpoints
rho3 <- 0
# Covariance matrix
S3 <- diag(SD3)
S3[1,2] <- rho3
S3[2,1] <- rho3


## Parameters of placebo
# Standardized effect size(dosage) in two consecutive endpoint
mu4 <- c(0.20, 0.32) # Considering that, dosage less than .40 has 
# no any Clinical response. 
# Standard deviation sequentially
SD4 <- c(1, .95)
# for simplification we consider there are no any
# correlation between endpoints
rho4 <- 0
# Covariance matrix
S4 <- diag(SD4)
S4[1,2] <- rho4
S4[2,1] <- rho4

# Here we save all parameters value in a file named ParametersData.RData.
# save(file = "ParametersData.RData",
#   S1, S2, S3, S4, alpha, M, mu1, mu2, mu3, mu4, n, rho1,
#  rho2, rho3, rho4, SD1, SD2, SD3, SD4)

### Test for first group or group a
pvals1 <- numeric(M)
pvals2 <- numeric(M)
for(i in seq_len(M)){
  X1 <- rmvnorm(n, mean = mu1, sigma = S1)
  X2 <- rmvnorm(n, mean = mu4, sigma = S4)
  
  pvals1[i] <- t.test(X1[,1], X2[,1])$p.value
  pvals2[i] <- t.test(X1[,2], X2[,2])$p.value
}

### Test for second group or group b
pvals3 <- numeric(M)
pvals4 <- rep(1, M)
for(i in seq_len(M)){
  X3 <- rmvnorm(n, mean = mu2, sigma = S2)
  X4 <- rmvnorm(n, mean = mu4, sigma = S4)
  
  pvals3[i] <- t.test(X3[,1], X4[,1])$p.value
  pvals4[i] <- t.test(X3[,2], X4[,2])$p.value
}

### Test for Theird group or group c
pvals5 <- numeric(M)
pvals6 <- rep(1, M)
for(i in seq_len(M)){
  X5 <- rmvnorm(n, mean = mu3, sigma = S3)
  X6 <- rmvnorm(n, mean = mu4, sigma = S4)
  
  pvals5[i] <- t.test(X5[,1], X6[,1])$p.value
  pvals6[i] <- t.test(X5[,2], X6[,2])$p.value
}

# considering two equally important primary endpoints in all three doses
# and here the FWER should be half of alpha(0.25).
# for dose a
# empirical power at first endpoint
sum(pvals1 < 0.5*alpha)/M
# empirical power at second endpoint
sum(pvals2 < 0.5*alpha)/M
# empirical power at least one endpoint
sum((pvals1 < 0.5*alpha) | (pvals2 < 0.5*alpha))/M


# for dose b
# empirical power at first endpoint
sum(pvals3 < 0.5*alpha)/M
# empirical power at second endpoint
sum(pvals4 < 0.5*alpha)/M
# empirical power at least one endpoint
sum((pvals3 < 0.5*alpha) | (pvals4 < 0.5*alpha))/M


# for dose c
# empirical power at first endpoint
sum(pvals5 < 0.5*alpha)/M
# empirical power at second endpoint
sum(pvals6 < 0.5*alpha)/M
# empirical power at least one endpoint
sum((pvals5 < 0.5*alpha) | (pvals6 < 0.5*alpha))/M


############ -----------------------
##########--------------------------

## Lets Consider fixed sequence procedure in all three doses
## for dose a
# empirical power at first endpoint
sum(pvals1 < alpha)/M
first.a <- pvals1 < alpha
# empirical power at second endpoint
sum(pvals2[first.a] < alpha)/M

###### for dose b
# empirical power at first endpoint
sum(pvals3 < alpha)/M
first.b <- pvals3 < alpha
# empirical power at second endpoint
sum(pvals4[first.b] < alpha)/M

## for dose c
# empirical power at first endpoint
sum(pvals5 < alpha)/M
first.c <- pvals5 < alpha
# empirical power at second endpoint
sum(pvals6[first.c] < alpha)/M

# In all three cases or all two endpoint fixed sequence procedure cases-
# empirical power at least one endpoint are identical to power of first endpoint
# and
# empirical power at both endpoints are identical to power of second endpoint.



##################------------------------
## Applying Fallback procedure in all three groups/dose
# weights
w1 <- 0.6
w2 <- 0.4

# dose a
# empirical power at first endpoint
sum(pvals1 < w1*alpha)/M
first <- pvals1 < w1*alpha ### ***(should be sum & divided by M)
# empirical power: second endpoint
(sum(pvals2[first] < alpha) + sum(pvals2[!first] < w2*alpha))/M
second <- pvals2[!first] < w2*alpha
# empirical power at least one endpoints
(sum(first) + sum(second))/M


# dose b
## empirical power at first endpoint
sum(pvals3 < w1*alpha)/M
theird <- pvals3 < w1*alpha
# empirical power at second endpoint
(sum(pvals4[theird] < alpha) + sum(pvals4[!theird] < w2*alpha))/M
fourth <- pvals4[!theird] < w2*alpha
# empirical power at least one or both endpoints
(sum(theird) + sum(fourth))/M


# dose c
# empirical power at first endpoint
sum(pvals5 < w1*alpha)/M
fifth <- pvals5 < w1*alpha
# empirical power at second endpoint
(sum(pvals6[fifth] < alpha) + sum(pvals6[!fifth] < w2*alpha))/M
sixth <- pvals6[!fifth] < w2*alpha
# empirical power at least one or both endpoints
(sum(fifth) + sum(sixth))/M

#######################----------------------------------
## Which dose is statistically most significant?
# The treatment effect hugely depends on the dose. Here, we have to carefull 
# about probability of geting a significant result erroneously. Best decision 
# is to take account size of the treatment effect of other groups considering 
# there dose. In our case, dose b is most statistically significant. 


# Then effects of weights on fallback procedures.
# If we change weights of endpoints in Fallback procedure the result 
# will also be change. Generally while weight decrase chance of getting 
# a significnt result on that specific endpoint also decrase and 
# while increase the consequence also same
.
# Here we give two example using dose b(We consider dose b in all 
# analysis as it is our most statistically significant).
w1 <- 0.4
w2 <- 0.6
sum(pvals3 < w1*alpha)/M
# [1] 0.775
theird <- pvals3 < w1*alpha
(sum(pvals4[theird] < alpha) + sum(pvals4[!theird] < w2*alpha))/M
# [1] 0.8
fourth <- pvals4[!theird] < w2*alpha
(sum(theird) + sum(fourth))/M
# [1] 0.943
# -----------------------------------------
# example 2:
w1 <- 0.5
w2 <- 0.5
sum(pvals3 < w1*alpha)/M
#[1] 0.804
theird <- pvals3 < w1*alpha
(sum(pvals4[theird] < alpha) + sum(pvals4[!theird] < w2*alpha))/M
# [1] 0.798
fourth <- pvals4[!theird] < w2*alpha
(sum(theird) + sum(fourth))/M
# [1] 0.945
# Here in our evaluation when both endpoints are equally important or first 
# endpoint is little bit more weighted it seems good to get a significant result on both 
# endpoint.
################################################################
### effects of correlation(rho) on power
# write theory in script; correlation:
#####################################
# Sequence of correlation score
s <- seq(from = .01, to = 1.0, by = 0.1)
pvals3 <- lengths(s)
pvals4 <- lengths(s)
# define length of power as length of correlation
# first endpoint by Fixed-sequence
power3 <- numeric(length(s))
# second endpoint by Fixed-sequence
power4 <- numeric(length(s))
# first endpoint by Fallback
power.fall3 <- numeric(length(s))
# second endpoint by Fallback
power.fall4 <- numeric(length(s))
# number of simulation
M <- 100000
for(i in 1:length(s)){
  sigma <- lengths(s)
  S2 <- matrix(c(0.87, s[i]*sqrt(0.87), s[i]*sqrt(0.87), 1), ncol = 2)
  S4 <- matrix(c(1, s[i]*sqrt(0.95), s[i]*sqrt(0.95), 0.95), ncol = 2)
  for(j in 1:M){
    X3 <- rmvnorm(n, mean = mu2, sigma = S2)
    X4 <- rmvnorm(n, mean = mu4, sigma = S4)
    
    pvals3[j] <- t.test(X3[,1], X4[,1])$p.value
    pvals4[j] <- t.test(X3[,2], X4[,2])$p.value
  }
  power3[i] <- sum(pvals3 < alpha)/M
  first.fixc <- pvals3 < alpha
  power4[i] <- sum(pvals4[first.fixc] < alpha)/M
  power.fall3[i] <- sum(pvals3 < alpha*w1)/M
  first.fallc <- pvals3 < alpha*w1
  power.fall4[i] <- (sum(pvals4[first.fallc] < alpha) + 
                       sum(pvals4[!first.fallc] < w2*alpha))/M
  second.fixc <- pvals4[first.fixc] < alpha
  second.fallc1 <- pvals4[first.fallc] < alpha
  second.fallc2 <- pvals4[!first.fallc] < w2*alpha
}
# correlation between power at first endpoint & correlation(rho)
# in fixsequence procedure
cor.test(s, power3, method = "spearman")
# correlation between power at first endpoint & correlation(rho)
# in fallback procedure
cor.test(s, power.fall3, method = "spearman")
# correlation between power at second endpoint & correlation(rho)
# in fixsequence procedure
cor.test(s, power4, method = "spearman")
# correlation between power at second endpoint & correlation(rho)
# in fallback procedure
cor.test(s, power.fall4, method = "spearman")

#save(file = "CorrelationData.RData",
#     s, power3, power.fall3, power4, power.fall4, first.fixc,
#      first.fallc, second.fallc1, second.fallc2, second.fixc)
load(file = "CorrelationData.RData")
fit <- lm(power3 ~ s)
fit
fit1 <- lm(power4 ~ s)
fit1
summary(fit)
summary(fit1)
fit2 <- lm(power.fall4 ~ s)
fit2
summary(fit2)

# Visual representation of relation between correlation(rho) & power; for this 
# purpose we apply grammar of graphics plot2 or ggplot2 package.
library(ggplot2)
CorrData <- data.frame(correlation = s, power.fa3 = power.fall3,
                       power.fa4 = power.fall4,
                       power3 = power3, power4 = power4)
# At first endpoint considering fixed-sequence procedure;
ggplot(CorrData, aes(x=correlation, y=power3)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation of correlation & power by Fixed-sequence") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("correlation(rho)") + ylab("power at first endpoint") + ylim(0.8, 0.9) +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))


# At second endpoint, 
ggplot(CorrData, aes(x=correlation, y=power4)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation of correlation & power by Fixed-sequence") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("correlation(rho)") + ylab("power at second endpoint") + ylim(0.7, 0.85) + 
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
# In-case of fallback procedure
# At first endpoint
ggplot(CorrData, aes(x=correlation, y=power.fa3)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between correlation & power by Fallback") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("correlation(rho)") + ylab("power at first endpoint") + ylim(0.8, 0.9) +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
  
# At second endpoint
ggplot(CorrData, aes(x=correlation, y=power.fa4)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between correlation & power by Fallback") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("correlation(rho)") + ylab("power at second endpoint") + ylim(0.7, 0.85) +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
#### ?geom_pointrange / ?geom_linerange.
# At second endpoint,  + + + ylim(0, 1)
ggplot(CorrData, aes(x=correlation, y=power4)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation of correlation & power by Fixed-sequence") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("correlation(rho)") + ylab("power at second endpoint") + ylim(0.5, 1) +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
#------------------------------
ggplot(CorrData, aes(x=correlation, y=power.fa4)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between correlation & power by Fallback") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("correlation(rho)") + ylab("power at second endpoint") + ylim(0.5, 1) +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))

# Here the result is unexpected. It's usual that correlation between 
# two endpoint has a reversible effect on power. But our graph tell us opposite direction.
# This result is almost unrealistic.

### (*** but in one simulation I got this result)
# cor(s,power3)
# [1] -0.8449522
# > cor(s,power4)
# [1] -0.7564213
# Optimistic.
################################################
#### effects of sample size on power
# Number of simulation
M <- 1000
# Sequence of sample size
n <- seq(from = 40, to = 150, by = 10)
power5 <- numeric(length(n))
power6 <- numeric(length(n))
power.fall5 <- numeric(length(n))
power.fall6 <- numeric(length(n))
# Covariance matrix considering no any correlation,
S2 <- matrix(c(0.87, 0, 0, 1), ncol = 2)
S4 <- matrix(c(1, 0, 0, 0.95), ncol = 2)
for(i in 1:length(n)){
  pvals5 <- numeric(M)
  pvals6 <- rep(1, M)
  for(j in 1:M){
    X3 <- rmvnorm(n[i], mean = mu2, sigma = S2)
    X4 <- rmvnorm(n[i], mean = mu4, sigma = S4)
    
    pvals5[j] <- t.test(X3[,1], X4[,1])$p.value
    pvals6[j] <- t.test(X3[,2], X4[,2])$p.value
  }
  power5[i] <- sum(pvals5 < alpha)/M
  first.fixs <- pvals5 < alpha
  power6[i] <- sum(pvals6[first.fixs] < alpha)/M
  power.fall5[i] <- sum(pvals5 < alpha*w1)/M
  first.ep <- pvals5 < alpha*w1
  power.fall6[i] <- (sum(pvals6[first.ep] < alpha) + sum(pvals6[!first.ep] < w2*alpha))/M
}
# At first endpoint by fixed-sequence
cor.test(n, power5, method = "spearman")
# At first endpoint by fallback
cor.test(n, power.fall5, method = "spearman")
# At second endpoint by fixed-sequence
cor.test(n, power6, method = "spearman")
# At second endpoint by fallback
cor.test(n, power.fall6, method = "spearman")

# Visual representation of relation between sample size & power
sampleData <- data.frame(samplesize = n, power.fa5 = power.fall5,
                         power.fa6 = power.fall6,
                         power5 = power5, power6 = power6)
# At first endpoint
ggplot(sampleData, aes(x=samplesize, y=power5)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation of sample size & power by fixed-sequence") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("sample size(n)") + ylab("power at first endpoint") + ylim(0.3, 1.0) +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
# At second endpoint
ggplot(sampleData, aes(x=samplesize, y=power6)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation of sample size & power by fixed-sequence") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("sample size(n)") + ylab("power at second endpoint") + ylim(0.3, 1.0) +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
# In-case of fallback procedure
# At first endpoint
ggplot(sampleData, aes(x=samplesize, y=power.fa5)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between sample size & power by fallback") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("sample size(n)") + ylab("power at first endpoint") + ylim(0.3, 1.0) +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
# At second endpoint
ggplot(sampleData, aes(x=samplesize, y=power.fa6)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between sample size & power by fallback") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("sample size(n)") + ylab("power at second endpoint") + ylim(0.3, 1.0) +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))

# Both pair sounds good. It's easy to understand positive effect of 
# sample size on power.  
###################################################
### effets of standardized effect size(mu) on power.
# At first we take consider standardized effect size of first endpoint. 
# Standardized effect size
mu2 <- cbind(seq(from = .25, to = 1, by = .05), 0.78)
# Number of sample size
n <- 80
power7 <- numeric(nrow(mu2))
power8 <- numeric(nrow(mu2))
power.fall7 <- numeric(nrow(mu2))
power.fall8 <- numeric(nrow(mu2))
for(i in 1:nrow(mu2)){
  pvals7 <- numeric(M)
  pvals8 <- rep(1, M)
  for(j in 1:M){
    X3 <- rmvnorm(n, mean = mu2[i,], sigma = S2)
    X4 <- rmvnorm(n, mean = mu4, sigma = S4)
    
    pvals7[j] <- t.test(X3[,1], X4[,1])$p.value
    pvals8[j] <- t.test(X3[,2], X4[,2])$p.value
  }
  power7[i] <- sum(pvals7 < alpha)/M
  first.fixe <- pvals7 < alpha
  power8[i] <- sum(pvals8[first.fixe] < alpha)/M
  power.fall7[i] <- sum(pvals7 < alpha*w1)/M
  first.epe <- pvals7 < alpha*w1
  power.fall8[i] <- (sum(pvals8[first.epe] < alpha) + sum(pvals8[!first.epe] < w2*alpha))/M
}

# At first endpoint by fixed-sequence
cor.test(mu2[1:16], power7, method = "spearman")
# At first endpoint by fallback
cor.test(mu2[1:16], power.fall7, method = "spearman")
# At second endpoint by fixed-sequence
cor.test(mu2[1:16], power8, method = "spearman")
# At second endpoint by fallback
cor.test(mu2[1:16], power.fall8, method = "spearman")
#**#(???) (here it's also a indicator that one-endpoint has a limited or almost no 
# effect on another).

# Visual representation of relation between standardized effect size & power
standrData <- data.frame(standardized = mu2[1:16],
                         power.fa7 = power.fall7, power.fa8 = power.fall8,
                         power7 = power7, power8 = power8)
# At first endpoint
ggplot(standrData, aes(x=standardized, y=power7)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between mu1 & power by fixed-sequence") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("First standardized effect size(mu1)") + 
  ylab("power at first endpoint") +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
# At second endpoint
ggplot(standrData, aes(x=standardized, y=power8)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between mu1 & power by fixed-sequence") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("First standardized effect size(mu1)") + 
  ylab("power at second endpoint") +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
# Incase of Fallback procedure,
# At first endpoint
ggplot(standrData, aes(x=standardized, y=power.fa7)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between mu1 & power by fallback") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("First standardized effect size(mu1)") + ylab("power at first endpoint") +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
# At second endpoint
ggplot(standrData, aes(x=standardized, y=power.fa8)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between mu1 & power by fallback") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("First standardized effect size(mu1)") + ylab("power at second endpoint") +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))

###(***) here one thing is noticable that,
# In case of fallback procedure *standardized effect size of first*
# endpoint has a clear effect on SDEs of second endpoint.
# At second endpoint
ggplot(standrData, aes(x=standardized, y=power8)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between mu1 & power by fixed-sequence") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("First standardized effect size(mu1)") + 
  ylab("power at second endpoint") + ylim(0, 1) +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
# Incase of Fallback procedure,
# At second endpoint
ggplot(standrData, aes(x=standardized, y=power.fa8)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between mu1 & power by fallback") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("First standardized effect size(mu1)") + 
  ylab("power at second endpoint") + ylim(0, 1)
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))



### Here we consider second standardized effect size.
# Standardized effect size  
mu2 <- cbind(0.68, seq(from = .25, to = 1, by = .05))
power9 <- numeric(nrow(mu2))
power10 <- numeric(nrow(mu2))
power.fall9 <- numeric(nrow(mu2))
power.fall10 <- numeric(nrow(mu2))
for(i in 1:nrow(mu2)){
  pvals9 <- numeric(M)
  pvals10 <- rep(1, M)
  for(j in 1:M){
    X3 <- rmvnorm(n, mean = mu2[i,], sigma = S2)
    X4 <- rmvnorm(n, mean = mu4, sigma = S4)
    
    pvals9[j] <- t.test(X3[,1], X4[,1])$p.value
    pvals10[j] <- t.test(X3[,2], X4[,2])$p.value
  }
  power9[i] <- sum(pvals9 < alpha)/M
  first.fixes <- pvals9 < alpha
  power10[i] <- sum(pvals10[first.fixes] < alpha)/M
  power.fall9[i] <- sum(pvals9 < alpha*w1)/M
  first.epe2 <- pvals9 < alpha*w1
  power.fall10[i] <- (sum(pvals10[first.epe2] < alpha) + sum(pvals10[!first.epe2] < w2*alpha))/M
}
cor.test(mu2[,2], power9, method = "spearman")
cor.test(mu2[,2], power.fall9, method = "spearman")
cor.test(mu2[,2], power10, method = "spearman")
cor.test(mu2[,2], power.fall10, method = "spearman")


# Visual representation of relation between standardized effect size & power
standrData2 <- data.frame(standardized2 = mu2[,2],
                          power.fa9 = power.fall9, power.fa10 = power.fall10,
                          power9 = power9, power10 = power10)
# At first endpoint
ggplot(standrData2, aes(x=standardized2, y=power9)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between mu2 & power by fixed-sequence") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("standardized effect size(mu2)") + ylab("power at first endpoint") + ylim(0.8, 0.9) +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
# At second endpoint
ggplot(standrData2, aes(x=standardized2, y=power10)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between mu2 & power by fixed-sequence") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("standardized effect size(mu2)") + ylab("power at second endpoint") + ylim(0.01, 1.0) +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
# In case of Fallback procedure
# At first endpoint
ggplot(standrData2, aes(x=standardized2, y=power.fa9)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between mu2 & power by fallback") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("standardized effect size(mu2)") + ylab("power at first endpoint") + ylim(0.8, 0.9) +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
# At second endpoint
ggplot(standrData2, aes(x=standardized2, y=power.fa10)) + 
  geom_point(shape=19, alpha=0.25) +
  geom_smooth(method = loess) +
  ggtitle("Relation between mu2 & power by fallback") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  xlab("standardized effect size(mu2)") + ylab("power at second endpoint") + ylim(0.01, 1.0) +
  theme(axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"))


# This standardized effect size & power graph is expected. Although it's a strange to 
# look that, standardized effect size of second endpoint has an influence on first endpoint. 

###################################
# sims <- rqmvnorm(100, mean = 1:2, sigma = diag(2))
# plot(sims)

######################################################
### three endpoint trial(dose vs placebo)
### For Gatekeeping procedure
## parameters for drug 
mu5 <- c(.65, .73, .62)
SD5 <- c(.8, .70, .75)
rho5 <- 0
S5 <- diag(SD5)
S5[1,2] <- rho5
S5[1,3] <- rho5
S5[2,1] <- rho5
S5[2,3] <- rho5
S5[3,1] <- rho5
S5[3,2] <- rho5
## parameters for placebo
mu6 <- c(.21, .3, .20)
SD6 <- c(.8, .70, .75)
rho6 <- 0
S6 <- diag(SD6)
S6[1,2] <- rho6
S6[1,3] <- rho6
S6[2,1] <- rho6
S6[2,3] <- rho6
S6[3,1] <- rho6
S6[3,2] <- rho6
###---------------------
pvals7 <- numeric(M)
pvals8 <- numeric(M)
pvals9 <- numeric(M)
for(i in seq_len(M)){
  X7 <- rmvnorm(n, mean = mu5, sigma = S5)
  X8 <- rmvnorm(n, mean = mu6, sigma = S6)
  
  pvals7[i] <- t.test(X7[,1], X8[,1])$p.value
  pvals8[i] <- t.test(X7[,2], X8[,2])$p.value
  pvals9[i] <- t.test(X7[,3], X8[,3])$p.value
}

# empirical power at first endpoint
## Applying Fixed sequence procedure
sum(pvals7 < alpha)/M
first <- pvals7 < alpha
# empirical power at second endpoint
sum(pvals8[first] < alpha)/M
second <- pvals8[first] < alpha
# empirical power at theird endpoint
sum(pvals9[second] < alpha)/M

##### Applying Fallback procedure
## Applying Fallback procedure
# weights
w1 <- 0.5
w2 <- 0.3
w3 <- 0.2

# empirical power at first endpoint
sum(pvals7 < w1*alpha)/M
first <- pvals7 < w1*alpha
# empirical power at second endpoint
(sum(pvals8[first] < (w1*alpha + w2*alpha)) + sum(pvals8[!first] < w2*alpha))/M
second <- pvals8[first] < (w1*alpha + w2*alpha)
theird <- pvals8[!first] < w2*alpha

# empirical power at theird endpoint
(sum((second < alpha) + (theird < (w2*alpha + w3*alpha)) + pvals9[!first] < w3*alpha))/M

#################################################################
second <- pvals2[!first] < w2*alpha
# empirical power at least one endpoints
(sum(first) + sum(second))/M


### By using gMCP package we can visualize different sorts of multiple endpoint 
# example. This package depends on Java.
library(gMCP)
graphGUI()
graph <- BonferroniHolm(3)
gMCP(graph, pvalues=c(0.01,0.07,0.02), alpha=0.05)

####################################################

