################################################################################
# Advanced Quantitative Methods Course                                         #
# H??vard Hegre, Mihai Croicu and David Randahl                                 #
# Lecture 2 File 1                                                             #
# Simple Monte Carlo experiment based on Carsey & Harden                       #
# Last update: 22 March 2019                                                   #
################################################################################

rm(list=ls(all=TRUE)) # Clear the workspace
# Import packages (remember to install them if they are new to your system)
library(xtable)
library(stargazer)

# Set the working directory

setwd("~/Dropbox/Apps/ShareLatex/Advanced Quantitative Methods/WD")
getwd()




# Basic OLS Example from C&H Chapters 1 and 4, now with 1000 reps and SEs
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
  Y <- b0 + b1*X + rnorm(n, 0, 1) # The true DGP, with N(0, 1) error
  model <- lm(Y ~ X) # Estimate OLS model
  vcv <- vcov(model) # Variance-covariance matrix
  par.est[i, 1] <- model$coef[1] # Put the estimate for the intercept
  # in the first column
  par.est[i, 2] <- model$coef[2] # Put the estimate for the coefficient on
  # X in the second column
  par.est[i, 3] <- sqrt(diag(vcv)[1]) # SE of the intercept
  par.est[i, 4] <- sqrt(diag(vcv)[2]) # SE of the coefficient on X
} # End the loop

# Coefficients
# Absolute Bias
ab.beta0 <- mean(abs(par.est[ , 1] - b0))
ab.beta1 <- mean(abs(par.est[ , 2] - b1)) 

ab.beta1

# MSE
mse.beta0 <- mean((par.est[ , 1] - b0)^2)
mse.beta1 <- mean((par.est[ , 2] - b1)^2)

mse.beta1

# Standard Errors
# Standard Deviation
sd.beta0 <- sd(par.est[ , 1]) # SD of the intercept estimates
mean.se.beta0 <- mean(par.est[ , 3]) # Mean SE of the intercept

sd.beta1 <- sd(par.est[ , 2]) # SD of the coefficient on X estimates
mean.se.beta1 <- mean(par.est[ , 4]) # Mean SE of the coefficient on X

#pdf("se-est1.pdf")

par(mar = c(5, 5.25, .5, .5))
hist(par.est[ , 3], breaks = 25, col = "purple", xlab = "", ylab = "",
     main = "", axes = FALSE)
axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression("Standard Error Estimates for"~beta[0]),
      cex.lab = 1.5)
title(ylab = expression("Frequency"), line = 3.75, cex.lab = 1.5)
abline(v = sd.beta0, lwd = 4)
text(.03275, 115, expression("SD of"~hat(beta[0])~"="~"0.0313973"),
     cex = 1.5)
box()

#dev.off()
#pdf("se-est2.pdf")

par(mar = c(5, 5.25, .5, .5))
hist(par.est[ , 4], breaks = 25, col = "blue", xlab = "", ylab = "",
     main = "", axes = FALSE)
axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression("Standard Error Estimates for"~beta[1]),
      cex.lab = 1.5)
title(ylab = expression("Frequency"), line = 3.75, cex.lab = 1.5)
abline(v = sd.beta1, lwd = 4)
text(.057, 150, expression("SD of"~hat(beta[1])~"="~"0.05501627"),
     cex = 1.5)
box()

#dev.off()

# CP Function
coverage <- function(b, se, true, level = .95, df = Inf){
  # b: Estimate, 
  # se: standard error,
  # true: true parameter, 
  # level: confidence level, 
  # df: degrees of freedom  
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
} # end of coverage function

cp.beta0 <- coverage(par.est[ , 1], par.est[ , 3], b0, df = n - model$rank)
cp.beta1 <- coverage(par.est[ , 2], par.est[ , 4], b1, df = n - model$rank)

cp.beta1$coverage.probability

#pdf("l5_cp_beta0.pdf")

par(mar = c(5, 6, .5, .5))
plot(seq(1, 100, length = 100), seq(.05, .4, length = 100), type = "n",
     axes = FALSE, xlab = "", ylab = "")
title(xlab = expression("100 Simulated Samples"), cex.lab = 1.5)
title(ylab = expression(hat(beta[0])), line = 3.75, cex.lab = 1.5)
box()
axis(1, at = seq(0, 100, 10), cex.axis = 1.25)
axis(2, at = seq(.05, 4, .05), cex.axis = 1.25, las = 2)
abline(h = b0, lwd = 2)
for (i in 1:100){
  points(i, par.est[i, 1], lwd = 2, col = ifelse(cp.beta0$true.in.ci[i] == 1,
                                                 "gray70", "gray20"), pch = 19)
  segments(i, cp.beta0$ci[i, 1], i, cp.beta0$ci[i, 2], lwd = 2,
           col = ifelse(cp.beta0$true.in.ci[i] == 1, "gray70", "gray20"))
}
legend("topleft", bty = "n", c(expression("CI includes true"~beta[0]),x??
                               expression("CI does not include true"~beta[0])),
       fill = c("gray70", "gray20"), cex = 1.5) 

#dev.off()

#pdf("l5_cp_beta1.pdf", width=10, height=6)

par(mar = c(5, 6, .5, .5))
plot(seq(100, 200, length = 100), seq(.25, .8, length = 100), type = "n",
     axes = FALSE, xlab = "", ylab = "")
title(xlab = expression("100 Simulated Samples"), cex.lab = 1.5)
title(ylab = expression(hat(beta[1])), line = 3.75, cex.lab = 1.5)
box()
axis(1, at = seq(100, 200, 10), labels = seq(0, 100, 10), cex.axis = 1.25)
axis(2, at = seq(.25, .75, .05), cex.axis = 1.25, las = 2)
abline(h = b1, lwd = 2)
for (i in 101:200){
  points(i, par.est[i, 2], lwd = 2, col = ifelse(cp.beta1$true.in.ci[i] == 1,
                                                 "gray70", "gray20"), pch = 19)
  segments(i, cp.beta1$ci[i, 1], i, cp.beta1$ci[i, 2], lwd = 2,
           col = ifelse(cp.beta1$true.in.ci[i] == 1, "gray70", "gray20"))
}
legend("topleft", bty = "n", c(expression("CI includes true"~beta[1]),
                               expression("CI does not include true"~beta[1])),
       fill = c("gray70", "gray20"), cex = 1.5) 
