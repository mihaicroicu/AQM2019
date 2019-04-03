################################################################################
# Monte Carlo Simulation and Resampling Methods for Social Science             #
# Thomas M. Carsey and Jeffrey J. Harden                                       #
# Chapter 6 File                                                               #
# Last update: 2/28/13                                                         #
################################################################################
# Set the working directory
setwd("C:/Users/jjharden/Dropbox/Research/Sim Book")
setwd("C:/Users/hegha560.000/Dropbox/Undervisning/UU/AdvancedQuantitativeMethods/WD")
setwd("C:/Users/hhegre/Dropbox/Undervisning/UU/AdvancedQuantitativeMethods/WD")

rm(list=ls(all=TRUE)) # Clear the workspace
library(ROCR)


# Maximum likelihood estimation illustration 

# An even more generic binomial PDF
pdf("pdf-binomial-4c.pdf")
tr <- 4
p <- 0.25
par(mar = c(5, 5.25, .5, .5))
plot(0:tr, dbinom(0:tr, tr, p), type = "h", lwd = 2, main = "", xlab = "",
     ylab = "", axes = FALSE,  ylim=c(0,0.5))
points(0:tr, dbinom(0:tr, tr, p), pch = 19)
axis(1, 0:tr, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression("X"), cex.lab = 1.5)
title(ylab = expression("P(X)"), line = 3.75, cex.lab = 1.5)
box()

dev.off()

# For lecture: 
set.seed(123456) # Set the seed for reproducible results
# estimates 
b0 <- .2 # True value for the intercept
b1 <- 1 # True value for the slope
n <- 1000 # Sample size
X <- runif(n, -5, 5) # Create a sample of n observations on the 
# independent variable X


# Logit
# Inverse Logit Function
inv.logit <- function(Z){
  return(exp(Z)/(1 + exp(Z)))
}

Z <- b0 + b1*X
p_logit <- inv.logit(b0 + b1*X)
p_probit <- pnorm(b0 + b1*X)
  
Y_logit <- rbinom(n, 1, inv.logit(b0 + b1*X)) # The true DGP, Bernoulli trials
Y_probit <- rbinom(n, 1, pnorm(b0 + b1*X)) # The true DGP, Bernoulli trials

Data <- cbind(X,p_logit,p_probit,Y_logit, Y_probit)
Data[1:10,]


pdf("link_functions.pdf", width=6, height = 3)
par(mfrow=c(1,2))
plot(Z,p_logit, abline(v=0, h=0.5, lwd = 1))
plot(Z,p_probit, abline(v=0, h=0.5, lwd = 1))
dev.off()

# OLS as a Probability Model
set.seed(123456) # Set the seed for reproducible results

reps <- 1000 # Set the number of repetitions at the top of the script
par.est <- matrix(NA, nrow = reps, ncol = 4) # Empty matrix to store the
                                             # estimates 
b0 <- .2 # True value for the intercept
b1 <- .5 # True value for the slope
n <- 1000 # Sample size
X <- runif(n, -1, 1) # Create a sample of n observations on the 
                     # independent variable X

for(i in 1:reps){ # Start the loop
Y <- rnorm(n, b0 + b1*X, 1) # The true DGP, Y ~ N(mu, sigma)
model <- lm(Y ~ X) # Estimate OLS model
vcv <- vcov(model) # Variance-covariance matrix
par.est[i, 1] <- model$coef[1] # Put the estimate for the intercept
                               # in the first column
par.est[i, 2] <- model$coef[2] # Put the estimate for the coefficient on
                               # X in the second column
par.est[i, 3] <- sqrt(diag(vcv)[1]) # SE of the intercept
par.est[i, 4] <- sqrt(diag(vcv)[2]) # SE of the coefficient on X
} # End the loop

# Logit
# Inverse Logit Function
inv.logit <- function(p){
  return(exp(p)/(1 + exp(p)))
}

# CP Function
coverage <- function(b, se, true, level = .95, df = Inf){ # Estimate, 
                                                          # standard error,
                                                          # true parameter, 
                                                          # confidence level, 
                                                          # and df  
qtile <- level + (1 - level)/2 # Compute the proper quantile
lower.bound <- b - qt(qtile, df = df)*se # Lower bound
upper.bound <- b + qt(qtile, df = df)*se # Upper bound 
# Is the true parameter in the confidence interval? (yes = 1)
true.in.ci <- ifelse(true >= lower.bound & true <= upper.bound, 1, 0)
cp <- mean(true.in.ci) # The coverage probability
mc.lower.bound <- cp - 1.96*sqrt((cp*(1 - cp))/length(b)) # Monte Carlo error  
mc.upper.bound <- cp + 1.96*sqrt((cp*(1 - cp))/length(b))  
return(list(coverage.probability = cp, # Return results
            true.in.ci = true.in.ci,
            ci = cbind(lower.bound, upper.bound),
            mc.eb = c(mc.lower.bound, mc.upper.bound)))
}

set.seed(32945) # Set the seed for reproducible results

reps <- 1000 # Set the number of repetitions at the top of the script
par.est.logit <- matrix(NA, nrow = reps, ncol = 4) # Empty matrix to store
                                                   # the estimates 
b0 <- .2 # True value for the intercept
b1 <- .5 # True value for the slope
n <- 1000 # Sample size
X <- runif(n, -1, 1) # Create a sample of n observations on the 
                     # independent variable X

for(i in 1:reps){ # Start the loop
Y <- rbinom(n, 1, inv.logit(b0 + b1*X)) # The true DGP, Bernoulli trials
model <- glm(Y ~ X, family = binomial (link = logit)) # Estimate logit model
vcv <- vcov(model) # Variance-covariance matrix
par.est.logit[i, 1] <- model$coef[1] # Put the estimate for the 
                                     # intercept in the first column
par.est.logit[i, 2] <- model$coef[2] # Put the estimate for the coefficient 
                                     # on X in the second column
par.est.logit[i, 3] <- sqrt(diag(vcv)[1]) # SE of the intercept
par.est.logit[i, 4] <- sqrt(diag(vcv)[2]) # SE of the coefficient on X
} # End the loop

cp.beta0.logit <- coverage(par.est.logit[ , 1], par.est.logit[ , 3], b0,
 df = n - model$rank)
cp.beta1.logit <- coverage(par.est.logit[ , 2], par.est.logit[ , 4], b1,
 df = n - model$rank)

pdf("logit-hist1.pdf")

par(mar = c(5, 5.25, .5, .5))
hist(par.est.logit[ , 1], breaks = 25, col = "gray50", ylim = c(0, 150),
 xlab = "", ylab = "", main = "", axes = FALSE)
axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression(hat(beta[0])), cex.lab = 1.5)
title(ylab = expression("Frequency"), line = 3.75, cex.lab = 1.5)
abline(v = b0, lwd = 4)
text(.05, 70, expression("True"~beta[0]~"= 0.20"), cex = 1.5)
box()

dev.off()

pdf("logit-hist2.pdf")

par(mar = c(5, 5.25, .5, .5))
hist(par.est.logit[ , 2], breaks = 25, xlim = c(.1, .9), col = "gray50",
 ylim = c(0, 200), xlab = "", ylab = "", main = "", axes = FALSE)
axis(1, at = seq(.1, .9, .1), cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression(hat(beta)[1]), cex.lab = 1.5)
title(ylab = expression("Frequency"), line = 3.75, cex.lab = 1.5)
abline(v = b1, lwd = 4)
text(.75, 125, expression("True"~beta[1]~"= 0.50"), cex = 1.5)
box()

dev.off()


# A prediction model
set.seed(4711) # Set the seed for reproducible results

n = 750

X1 <- runif(n, -1, 1) 
X2 <- runif(n, -1, 1) # 
X3 <- runif(n, -1, 1) 

b0 = 0
b1 = -1
b2 = 1
b3 = .5

Y <- rbinom(n, 1, inv.logit(b0 + b1*X1 + b2*X2 + b3*X3)) # The true DGP, Bernoulli trials

df = data.frame(cbind(Y,X1,X2,X3))

df_train <- df[1:500,]

df_test <- df[741:750,]

logit_model <- glm(Y ~ X1 + X2 + X3, family = binomial (link = logit), data=df_train) # Estimate logit model

summary(logit_model)


df_test$predprob <- predict(logit_model,newdata=df_test, type="response")
df_test$pred <- 0
df_test$pred[df_test$predprob>0.5] <- 1

summary(df_test$Y)
summary(df_test$pred)

conftable <- table(df_test$pred,df_test$Y)
addmargins(conftable)

test_sorted <- df_test[order(-df_test$predprob),] 
head(test_sorted, 25)

tpr = cbind(0,0,0,0,0,0,0,0,0,0)
fpr = cbind(0,0,0,0,0,0,0,0,0,0)
plot(fpr, tpr)

# ROC and PR curves:

library(ROCR)
pred <- prediction(df_test$predprob, df_test$Y)
roc <- performance(pred,"tpr","fpr")
plot(roc,colorize=TRUE)
pred
pr <- performance(pred,"prec","rec")
plot(pr,colorize=TRUE)

# Top 5 
