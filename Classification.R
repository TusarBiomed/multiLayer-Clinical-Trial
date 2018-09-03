## ----message=FALSE-------------------------------------------------------
library(ROCR)
data(ROCR.hiv)
attach(ROCR.hiv)
pred.svm <- prediction(hiv.svm$predictions, hiv.svm$labels)
perf.svm <- performance(pred.svm, 'tpr', 'fpr')
pred.nn <- prediction(hiv.nn$predictions, hiv.svm$labels)
perf.nn <- performance(pred.nn, 'tpr', 'fpr')
plot(perf.svm, lty=1, col="red",
     main="SVMs and NNs for prediction of HIV-1 coreceptor usage")
plot(perf.nn, lty=1, col="blue", add=TRUE)
plot(perf.svm, avg="vertical", lwd=2, col="red",
     spread.estimate="stderror",plotCI.lwd=2,add=TRUE)
plot(perf.nn, avg="vertical", lwd=2, col="blue",
     spread.estimate="stderror",plotCI.lwd=2,add=TRUE)
legend("bottomright",c('SVM','NN'),col=c('red','blue'),lwd=3)
detach(ROCR.hiv)

## ----fig.width = 10------------------------------------------------------
set.seed(13)
x <- matrix(rnorm(3*18), ncol = 3)
x[16:18,1] <- x[16:18,1] + c(2.5, 1.5, 2) + 0.1
x[13:15,2] <- x[13:15,2] + c(2.5, 1.5, 2.5) + 0.1
x[10:12,3] <- x[10:12,3] - c(1.5, 1.5, 0.5) - 0.1
x <- x+7
group <- factor(c(rep("control", 9), rep("patient", 9)))
par(mfrow = c(1, 3))
plot(x[,1], pch = 19, col = c(rep("blue", 9), rep("red", 9)),
     main = expression(paste("Feature ", x[1])),
     ylab = expression(x[1]))
plot(x[,2], pch = 19, col = c(rep("blue", 9), rep("red", 9)),
     main = expression(paste("Feature ", x[2])),
     ylab = expression(x[2]))
plot(x[,3], pch = 19, col = c(rep("blue", 9), rep("red", 9)),
     main = expression(paste("Feature ", x[3])),
     ylab = expression(x[3]))

## ----feature1------------------------------------------------------------
t.test(x[,1] ~ group)
t.test(x[,2] ~ group)
t.test(x[,3] ~ group)

## ----message=FALSE-------------------------------------------------------
library(mclust)
library(RColorBrewer)
data(banknote)
mypalette <- brewer.pal(3, "Set1")[c(1,3)]
klassen <- banknote$Status
pairs(banknote[,-1], col = mypalette[klassen],
      main = "banknote-dataset", pch = 20)

## ------------------------------------------------------------------------
banknote.sel <- banknote[,c("Bottom", "Top")]
plot(banknote.sel, pch = 20, col = mypalette[klassen],
     main = "banknote-dataset - selected variables")

## ------------------------------------------------------------------------
mu1 <- colMeans(banknote.sel[klassen == "genuine",])
mu2 <- colMeans(banknote.sel[klassen == "counterfeit",])
Sigma.inv <- solve(cov(banknote.sel))
fun <- function(x1, mu1, mu2, Sigma.inv, N1, N2){
    f <- function(x2, x1, mu1, mu2, Sigma.inv, N1, N2){
        (t(c(x1, x2)) %*% Sigma.inv %*% (mu2 - mu1)
         - 0.5*(t(mu2) %*% Sigma.inv %*% mu2
                - t(mu1) %*% Sigma.inv %*% mu1)
         - log(N1/N2))
    }
    uniroot(f, c(0, 15), x1 = x1, mu1 = mu1, mu2 = mu2,
            Sigma.inv = Sigma.inv, N1 = N1, N2 = N2)$root
}
bottom.vec.lda <- c(7, 13)
top.vec.lda <- sapply(bottom.vec.lda, fun, mu1 = mu1, mu2 = mu2,
                      Sigma.inv = Sigma.inv, N1 = 100, N2 = 100)

## ------------------------------------------------------------------------
plot(banknote.sel, col = mypalette[klassen], lwd = 2,
     main = "Discriminant bound for LDA", pch = 20)
lines(bottom.vec.lda, top.vec.lda, lwd = 2, lty = 1,
      col = mypalette[1])
lines(bottom.vec.lda, top.vec.lda, lwd = 2, lty = 2,
      col = mypalette[2])
text(8, 8, "genuine", font = 2, col = mypalette[2], cex = 2)
text(11.5, 12, "counterfeit", font = 2, col = mypalette[1], cex = 2)

## ------------------------------------------------------------------------
Sigma.inv <- matrix(0, ncol = 2, nrow = 2)
Sigma.inv[1,1] <- 1/var(banknote.sel[,"Bottom"])
Sigma.inv[2,2] <- 1/var(banknote.sel[,"Top"])
fun <- function(x1, mu1, mu2, Sigma.inv, N1, N2){
    f <- function(x2, x1, mu1, mu2, Sigma.inv, N1, N2){
        (t(c(x1, x2)) %*% Sigma.inv %*% (mu2 - mu1)
         - 0.5*(t(mu2) %*% Sigma.inv %*% mu2
                - t(mu1) %*% Sigma.inv %*% mu1)
         - log(N1/N2))
    }
    uniroot(f, c(0, 15), x1 = x1, mu1 = mu1, mu2 = mu2,
            Sigma.inv = Sigma.inv, N1 = N1, N2 = N2)$root
}
bottom.vec.dlda <- c(7, 13)
top.vec.dlda <- sapply(bottom.vec.dlda, fun, mu1 = mu1, mu2 = mu2,
                       Sigma.inv = Sigma.inv, N1 = 100, N2 = 100)

## ------------------------------------------------------------------------
plot(banknote.sel, col = mypalette[klassen], lwd = 2,
     main = "Discriminant bound for DLDA")
lines(bottom.vec.dlda, top.vec.dlda, lwd = 2, lty = 1,
      col = mypalette[1])
lines(bottom.vec.dlda, top.vec.dlda, lwd = 2, lty = 2,
      col = mypalette[2])
text(8, 8, "genuine", font = 2, col = mypalette[2], cex = 2)
text(11.5, 12, "counterfeit", font = 2, col = mypalette[1], cex = 2)

## ------------------------------------------------------------------------
Sigma1 <- cov(banknote.sel[klassen == "genuine",])
Sigma1.inv <- solve(Sigma1)
Sigma2 <- cov(banknote.sel[klassen == "counterfeit",])
Sigma2.inv <- solve(Sigma2)
fun <- function(x1, mu1, mu2, Sigma1, Sigma1.inv, Sigma2,
                Sigma2.inv, N1, N2, int){
    f <- function(x2, x1, mu1, mu2, Sigma1, Sigma1.inv, Sigma2,
                  Sigma2.inv, N1, N2){
        x <- c(x1, x2)
        D1 <- (- 0.5*log(det(Sigma1))
               - 0.5*t(x - mu1) %*% Sigma1.inv %*% (x - mu1)
               + log(N1/(N1 + N2)))
        D2 <- (- 0.5*log(det(Sigma2))
               - 0.5*t(x - mu2) %*% Sigma2.inv %*% (x - mu2)
               + log(N2/(N1 + N2)))
        return(D2 - D1)
    }
    uniroot(f, int, x1 = x1, mu1 = mu1, mu2 = mu2,
            Sigma1 = Sigma1, Sigma1.inv = Sigma1.inv,
            Sigma2 = Sigma2, Sigma2.inv = Sigma2.inv,
            N1 = 100, N2 = 100)$root
}
bottom.vec <- seq(7, 13, 0.1)
top.vec <- sapply(bottom.vec, fun, mu1 = mu1, mu2 = mu2,
                  Sigma1 = Sigma1, Sigma1.inv = Sigma1.inv,
                  Sigma2 = Sigma2, Sigma2.inv = Sigma2.inv,
                  N1 = 100, N2 = 100, int = c(5, 15))
top.fun.qda <- splinefun(bottom.vec, top.vec, method = "natural")
bottom.vec.qda <- seq(7, 13, 0.01)

## ------------------------------------------------------------------------
plot(banknote.sel, col = mypalette[klassen], lwd = 2,
     main = "Discriminant bound for QDA")
lines(bottom.vec.qda, top.fun.qda(bottom.vec.qda), lwd = 2, lty = 1,
      col = mypalette[1])
lines(bottom.vec.qda, top.fun.qda(bottom.vec.qda), lwd = 2, lty = 2,
      col = mypalette[2])
text(8, 8, "genuine", font = 2, col = mypalette[2], cex = 2)
text(11.5, 12, "counterfeit", font = 2, col = mypalette[1], cex = 2)

## ------------------------------------------------------------------------
Sigma1 <- diag(apply(banknote.sel[klassen == "genuine",], 2, var))
Sigma1.inv <- solve(Sigma1)
Sigma2 <- diag(apply(banknote.sel[klassen == "counterfeit",], 2, var))
Sigma2.inv <- solve(Sigma2)
fun <- function(x1, mu1, mu2, Sigma1, Sigma1.inv, Sigma2, Sigma2.inv,
                N1, N2, int){
    f <- function(x2, x1, mu1, mu2, Sigma1, Sigma1.inv,
                  Sigma2, Sigma2.inv, N1, N2){
        x <- c(x1, x2)
        D1 <- (- 0.5*log(det(Sigma1))
               - 0.5*t(x - mu1) %*% Sigma1.inv %*% (x - mu1)
               + log(N1/(N1 + N2)))
        D2 <- (- 0.5*log(det(Sigma2))
               - 0.5*t(x - mu2) %*% Sigma2.inv %*% (x - mu2)
               + log(N2/(N1 + N2)))
        return(D2 - D1)
    }
    uniroot(f, int, x1 = x1, mu1 = mu1, mu2 = mu2,
            Sigma1 = Sigma1, Sigma1.inv = Sigma1.inv,
            Sigma2 = Sigma2, Sigma2.inv = Sigma2.inv,
            N1 = N1, N2 = N2)$root
}
bottom.vec <- seq(7, 13, 0.1)
top.vec <- sapply(bottom.vec, fun, mu1 = mu1, mu2 = mu2,
                  Sigma1 = Sigma1, Sigma1.inv = Sigma1.inv,
                  Sigma2 = Sigma2, Sigma2.inv = Sigma2.inv,
                  N1 = 100, N2 = 100, int = c(0, 20))
top.fun.dqda <- splinefun(bottom.vec, top.vec, method = "natural")
bottom.vec.dqda <- seq(7, 13, 0.01)

## ------------------------------------------------------------------------
plot(banknote.sel, col = mypalette[klassen], lwd = 2,
     main = "Discriminant bound for DQDA")
lines(bottom.vec.dqda, top.fun.dqda(bottom.vec.dqda), lwd = 2,
      lty = 1, col = mypalette[1])
lines(bottom.vec.dqda, top.fun.dqda(bottom.vec.dqda), lwd = 2,
      lty = 2, col = mypalette[2])
text(8, 8, "genuine", font = 2, col = mypalette[2], cex = 2)
text(11.5, 12, "counterfeit", font = 2, col = mypalette[1], cex = 2)

## ------------------------------------------------------------------------
mypalette1 <- brewer.pal(4, "Dark2")
plot(banknote.sel, col = mypalette[klassen], lwd = 2,
     main = "Discriminant bounds")
lines(bottom.vec.lda, top.vec.lda, lwd = 2, lty = 1,
      col = mypalette1[1])
lines(bottom.vec.dlda, top.vec.dlda, lwd = 2, lty = 1,
      col = mypalette1[2])
lines(bottom.vec.qda, top.fun.qda(bottom.vec.qda), lwd = 2,
      lty = 1, col = mypalette1[3])
lines(bottom.vec.dqda, top.fun.dqda(bottom.vec.dqda), lwd = 2,
      lty = 1, col = mypalette1[4])
legend("topright", legend = c("LDA", "DLDA", "QDA", "DQDA"),
       fill = mypalette1)
text(8, 8, "genuine", font = 2, col = mypalette[2], cex = 2)
text(11.5, 9.75, "counterfeit", font = 2, col = mypalette[1], cex = 2)

## ------------------------------------------------------------------------
#source("http://bioconductor.org/biocLite.R")
#biocLite("gpls")
library(gpls)
fit <- gpls(as.integer(klassen)-1 ~ Bottom + Top,
            data = banknote.sel)
bottom.vec.gpls <- c(7, 13)
top.vec.gpls <- (0.5 - fit$coeff[1]
                 - fit$coeff[2]*bottom.vec.gpls)/fit$coeff[3]
plot(banknote.sel, col = mypalette[klassen], lwd = 2,
     main = "Discriminant bound for GPLS")
lines(bottom.vec.gpls, top.vec.gpls, lwd = 2, lty = 1,
      col = mypalette[1])
lines(bottom.vec.gpls, top.vec.gpls, lwd = 2, lty = 2,
      col = mypalette[2])
text(8, 8, "genuine", font = 2, col = mypalette[2], cex = 2)
text(11.5, 12, "counterfeit", font = 2, col = mypalette[1], cex = 2)

## ------------------------------------------------------------------------
library(class)
fit <- knn(banknote.sel, banknote.sel, cl = klassen, k = 5)
plot(banknote.sel, col = mypalette[klassen], type = "n", lwd = 2,
     main = "Classification by 5-NN")
text(banknote.sel, as.character(fit), cex = 0.7,
     col = mypalette[klassen])

## ------------------------------------------------------------------------
library(e1071)
knn.obj <- tune.knn(banknote.sel, klassen, k = 1:10,
                    tunecontrol = tune.control(sampling = "boot"))
summary(knn.obj)

## ------------------------------------------------------------------------
plot(knn.obj, main = "Optimal k for banknote dataset")

## ------------------------------------------------------------------------
library(rpart)
fit <- rpart(klassen ~ Bottom + Top, data = banknote.sel,
             method = "class", control = rpart.control(minsplit = 2))
plot(fit, uniform = TRUE, margin = 0.1, branch = 1,
     main = "Decision tree for true vs. false")
text(fit, use.n = TRUE, all = TRUE, fancy = TRUE, fwidth = 0.7,
     fheight = 0.9)

## ----message=FALSE-------------------------------------------------------
library(randomForest)
randomForest(klassen ~ Bottom + Top, data = banknote.sel)

## ------------------------------------------------------------------------
library(e1071)
data.svm <- data.frame(klassen = klassen, banknote.sel)
fit <- svm(klassen ~ Top + Bottom, data = data.svm)
plot(fit, data.svm, symbolPalette = mypalette[c(1,2)],
     svSymbol = "o", grid = 200,
     col = brewer.pal(3, "Pastel1")[c(1,3)])

## ------------------------------------------------------------------------
svm.obj <- tune.svm(banknote.sel, klassen,
                    gamma = c(0.2, 0.5, 1, 2, 5),
                    cost = 2^(1:4))
summary(svm.obj)

## ------------------------------------------------------------------------
plot(svm.obj, main = "Optimale parameters for radial kernel")

## ------------------------------------------------------------------------
fit <- svm(klassen ~ Top + Bottom, data = data.svm,
           gamma = 0.5, cost = 2)
plot(fit, data.svm, symbolPalette = mypalette[c(1,2)],
     svSymbol = "o", grid = 200,
     col = brewer.pal(3, "Pastel1")[c(1,3)])

