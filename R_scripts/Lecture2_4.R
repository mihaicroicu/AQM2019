################################################################################
# Advanced Quantitative Methods Course                                         #
# H??vard Hegre, Mihai Croicu and David Randahl                                 #
# Lecture 2 File 4                                                             #
# Omitted variable bias                                                        #
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


# Omitted Variable
set.seed(37943) # Set the seed for reproducible results
library(mvtnorm)


reps <- 1000 # Set the number of repetitions at the top of the script

b0 <- .2 # True value for the intercept
b1 <- .5 # True value for the slopes
b2 <- .75
n <- 1000 # Sample size

# Level of IV correlation
cor.level <- c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, .99) 
# Empty matrix to store the estimates 
par.est.ov <- matrix(NA, nrow = reps, ncol = length(cor.level)) 
# Empty matrix to store mean squared error
mse.ov <- matrix(NA, nrow = length(cor.level), ncol = 1) 

for(j in 1:length(cor.level)){ # Start the j loop
  for(i in 1:reps){ # Start the loop
    X.corr <- matrix(c(1, cor.level[j], cor.level[j], 1), nrow = 2, ncol = 2)
    X <- rmvnorm(n, mean = c(0, 0), sigma = X.corr) # Create two correlated 
    X1 <- X[ , 1]                                   # independent variables
    X2 <- X[ , 2]
    Y <- b0 + b1*X1 + b2*X2 + rnorm(n, 0, 1) # The true DGP, with N(0, 1) error
    model <- lm(Y ~ X1) # Estimate OLS model
    par.est.ov[i, j] <- model$coef[2] # Put the estimate for the coefficient on
    # X1 in column j
  } # End the i loop
  mse.ov[j] <- mean((par.est.ov[ , j] - b1)^2)
  cat("Completed cor =", cor.level[j], "\n")
} # End the j loop


pdf("ov-density.pdf")

par(mar = c(5, 5.25, .5, .5))
plot(density(par.est.ov[ , 1]), lty = 1, xlim = c(.25, 1.5), ylim = c(0, 16),
     lwd = 3, xlab = "", ylab = "", main = "", axes = FALSE)
lines(density(par.est.ov[ , 11]), lwd = 3, lty = 2)
axis(1, at = seq(.25, 1.5, .25), cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression(hat(beta[1])), cex.lab = 1.5)
title(ylab = expression("Density"), line = 3.75, cex.lab = 1.5)
abline(v = b1, lwd = 2)
text(.8, 10, expression("True"~beta[1]~"= 0.50"), cex = 1.5)
box()
legend("topright", bty = "n", c(expression(r[X[1]~X[2]]~"= 0"),
                                expression(r[X[1]~X[2]]~"= 0.99")), lty = c(1, 2), lwd = 2, cex = 1.5)

dev.off()

pdf("ov-mse.pdf")

par(mar = c(5, 5.5, .5, .5))
plot(cor.level, mse.ov, type = "b", xlim = c(0, 1), lwd = 3, xlab = "",
     ylab = "", axes = FALSE)
axis(1, at = cor.level, cex.axis = 1)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression("Correlation between"~X[1]~"and"~X[2]),
      cex.lab = 1.5)
title(ylab = expression(hat(beta[1])~"Mean Squared Error"), line = 3.5,
      cex.lab = 1.5)
box()

dev.off()




