################################################################################
# Advanced Quantitative Methods Course                                         #
# H??vard Hegre, Mihai Croicu and David Randahl                                 #
# Lecture 2 File 3                                                             #
# Example: Measurement error
# Last update: 22 March 2019                                                   #
################################################################################


rm(list=ls(all=TRUE)) # Clear the workspace
# Import packages (remember to install them if they are new to your system)
library(xtable)
library(stargazer)

# Set the working directory

setwd("~/Dropbox/Apps/ShareLatex/Advanced Quantitative Methods/WD")
getwd()



# Measurement Error
set.seed(385062) # Set the seed for reproducible results

reps <- 1000 # Set the number of repetitions at the top of the script
b0 <- .2 # True value for the intercept
b1 <- .5 # True value for the slope
n <- 1000 # Sample size
X <- runif(n, -1, 1) # Create a sample of n observations on the 
# independent variable X

# Level of measurement error (SD of random noise)
e.level <- c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1) 
# Empty matrix to store absolute bias
ab.merror <- matrix(NA, nrow = reps, ncol = length(e.level)) 

# Empty array to store the estimates                                                
par.est.merror <- array(NA, c(reps, 2, length(e.level)))

for(j in 1:length(e.level)){ # Start the j loop
  ab.1 <- numeric(reps)
  par.est <- matrix(NA, nrow = reps, ncol = 2) # Empty matrix to store the
  # estimates      
  Xp <- X + rnorm(n, 0, e.level[j]) # X measured with error
  
  for(i in 1:reps){ 
    Y <- b0 + b1*X + rnorm(n, 0, 1) 
    model <- lm(Y ~ Xp) 
    par.est[i, 1] <- model$coef[1] # Put the estimate for the intercept
    # in the first column
    par.est[i, 2] <- model$coef[2] # Put the estimate for the coefficient on
    # X in the second column
    ab.1[i] <- abs(model$coef[2] - b1)
  } # End of i loop
  
  par.est.merror[ , , j] <- par.est
  ab.merror[ , j] <- ab.1
  cat("Completed e =", e.level[j], "\n")
  gc() # Clear out RAM for better performance 
} # End of j loop

pdf("merror-density.pdf")

par(mar = c(5, 5.25, .5, .5))
plot(density(par.est.merror[ , 2, 1]), lty = 1, xlim = c(0, 1),
     ylim = c(0, 15), lwd = 3, xlab = "", ylab = "", main = "", axes = FALSE)
lines(density(par.est.merror[ , 2, 11]), lwd = 3, lty = 2)
axis(1, at = seq(0, 1, .1), cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression(hat(beta[1])), cex.lab = 1.5)
title(ylab = expression("Density"), line = 3.75, cex.lab = 1.5)
abline(v = b1, lwd = 2)
text(.75, 7, expression("True"~beta[1]~"= 0.50"), cex = 1.5)
box()
legend("topright", bty = "n", c(expression(sigma[ME]~"= 0"),
                                expression(sigma[ME]~"= 1")), lty = c(1, 2), lwd = 3, cex = 1.5)

dev.off()

pdf("merror-ab.pdf")

par(mar = c(5, 6, .5, .5))
plot(rep(e.level[1], times = reps), ab.merror[ , 1], xlim = c(0, 1),
     ylim = c(0, .5), col = "gray60", xlab = "", ylab = "", axes = FALSE)
for(i in 2:length(e.level)){
  points(rep(e.level[i], times = reps), ab.merror[ , i], col = "gray60")
  lines(lowess(e.level, apply(ab.merror, 2, mean), f = .2), lwd = 3)
}
axis(1, at = e.level, cex.axis = 1)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression("SD of Measurement Error"), cex.lab = 1.5)
title(ylab = expression(hat(beta[1])~"Absolute Bias"), line = 3.75,
      cex.lab = 1.5)
box()

dev.off()

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