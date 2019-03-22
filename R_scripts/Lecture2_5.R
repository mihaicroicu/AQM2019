################################################################################
# Advanced Quantitative Methods Course                                         #
# H??vard Hegre, Mihai Croicu and David Randahl                                 #
# Lecture 2 File 4                                                             #
# Multicollinearity                                                            #
# Last update: 22 March 2019                                                   #
################################################################################


rm(list=ls(all=TRUE)) # Clear the workspace
# Import packages (remember to install them if they are new to your system)
library(xtable)
library(stargazer)
library(mvtnorm)


# Set the working directory

setwd("~/Dropbox/Apps/ShareLatex/Advanced Quantitative Methods/WD")
getwd()


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
} # end of coverage function



# Multicollinearity
library(mvtnorm)
set.seed(121402) # Set the seed for reproducible results

reps <- 5000 # Set the number of repetitions at the top of the script
par.est.mc <- matrix(NA, nrow = reps, ncol = 4) # Empty matrix to store the
# estimates 
b0 <- .2 # True value for the intercept
b1 <- .5 # True value for the slopes
b2 <- .75 
n <- 500 # Sample size

# Levels of multicollinearity
mc.level <- c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, .99) 
# Matrix to store SD of the coefficient estimates. 
sd.betas <- matrix(NA, nrow = length(mc.level), ncol = 2)

# Matrix to store mean  of the standard error estimates. 
mean.se <- matrix(NA, nrow = length(mc.level), ncol = 2)


for(j in 1:length(mc.level)){ # Start the j loop
  X.corr <- matrix(c(1, mc.level[j], mc.level[j], 1), nrow = 2, ncol = 2)
  X <- rmvnorm(n, mean = c(0, 0), sigma = X.corr) # Create two correlated 
  X1 <- X[ , 1]                                   # independent variables
  X2 <- X[ , 2]
  
  for(i in 1:reps){ # Start the i loop
    Y <- b0 + b1*X1 + b2*X2 + rnorm(n, 0, 1) # The true DGP, with N(0, 1) error
    model <- lm(Y ~ X1 + X2) # Estimate OLS model
    vcv <- vcov(model) # Variance-covariance matrix
    par.est.mc[i, 1] <- model$coef[2] # Put the estimate for the coefficient on
    # X1 in the first column
    par.est.mc[i, 2] <- model$coef[3] # Put the estimate for the coefficient on
    # X2 in the second column
    par.est.mc[i, 3] <- sqrt(diag(vcv)[2]) # SE of the coefficient on X1
    par.est.mc[i, 4] <- sqrt(diag(vcv)[3]) # SE of the coefficient on X2
  } # End the i loop
  
  # I have corrected what seems an error in the original in the statement below and rather collect sd.betas for both Xs. 
  sd.betas[j, ] <- c(sd(par.est.mc[ , 1]), sd(par.est.mc[ , 2])) 
  
  mean.se[j, ] <- c(mean(par.est.mc[ , 3]), mean(par.est.mc[ , 4]))
  cat("Just completed correlation =", mc.level[j], 
      "(", j, "of", length(mc.level), ")", "\n")
} # End the j loop


pdf("mc-example_5000.pdf")

par(mar = c(5, 5.5, .5, .5))
plot(mc.level, sd.betas[ , 1], lwd = 3, ylim = c(0, .25), type = "b",
     xlab = "", ylab = "", main = "", axes = FALSE)
axis(1, at = mc.level, cex.axis = 1)
axis(2, at = seq(0, .25, .05), cex.axis = 1.25, las = 2)
title(xlab = expression("Correlation between"~X[1]~"and"~X[2]),
      cex.lab = 1.5)
title(ylab = expression("SD of"~hat(beta[1])), line = 3.75, cex.lab = 1.5)
box()

dev.off()


pdf("mc-example_se_5000.pdf")

par(mar = c(5, 5.5, .5, .5))
plot(mc.level, mean.se[ , 1], lwd = 3, ylim = c(0, .25), type = "b",
     xlab = "", ylab = "", main = "", axes = FALSE)
axis(1, at = mc.level, cex.axis = 1)
axis(2, at = seq(0, .25, .05), cex.axis = 1.25, las = 2)
title(xlab = expression("Correlation between"~X[1]~"and"~X[2]),
      cex.lab = 1.5)
title(ylab = expression("mean of se for beta1"), line = 3.75, cex.lab = 1.5)
box()

dev.off()


