################################################################################
# Monte Carlo Simulation and Resampling Methods for Social Science             #
# Based on Thomas M. Carsey and Jeffrey J. Harden                              #
# Ordered logistic regression                                                  #
# Last update: 4/10/19                                                         #
################################################################################
# Set the working directory
setwd("C:/Users/jjharden/Dropbox/Research/Sim Book")
setwd("C:/Users/hegha560.000/Dropbox/Undervisning/UU/AdvancedQuantitativeMethods/WD")
setwd("C:/Users/hhegre/Dropbox/Undervisning/UU/AdvancedQuantitativeMethods/WD")

rm(list=ls(all=TRUE)) # Clear the workspace


# Count Models
# Poisson 
set.seed(3759) # Set the seed for reproducible results

reps <- 1000 # Set the number of repetitions at the top of the script
par.est.pois <- matrix(NA, nrow = reps, ncol = 4) # Empty matrix to store the
# estimates 
b0 <- .2 # True value for the intercept
b1 <- .5 # True value for the slope
n <- 1000 # Sample size
X <- runif(n, -1, 1) # Create a sample of n observations on the 
# independent variable X

for(i in 1:reps){ # Start the loop
  Y <- rpois(n, exp(b0 + b1*X)) # The true DGP
  model <- glm(Y ~ X, family = "poisson") # Estimate Poisson model
  vcv <- vcov(model) # Variance-covariance matrix
  par.est.pois[i, 1] <- model$coef[1] # Put the estimate for the intercept
  # in the first column
  par.est.pois[i, 2] <- model$coef[2] # Put the estimate for the coefficient on
  # X in the second column
  par.est.pois[i, 3] <- sqrt(diag(vcv)[1]) # SE of the intercept
  par.est.pois[i, 4] <- sqrt(diag(vcv)[2]) # SE of the coefficient on X
} # End the loop

# Means of the coefficient estimates
mean(par.est.pois[ , 1]) # Intercept
mean(par.est.pois[ , 2]) # Coefficient on X

# Coverage probabilities 
# Intercept
coverage(par.est.pois[ , 1], par.est.pois[ , 3], b0,
         df = n - model$rank)$coverage.probability
# Coefficient on X
coverage(par.est.pois[ , 2], par.est.pois[ , 4], b1,
         df = n - model$rank)$coverage.probability

pdf("od1.pdf")

par(mar = c(5, 5.25, .5, .5)) 
plot(1:n, Y, ylim = c(0, 20), col = "gray50", xlab = "", ylab = "",
     main = "", axes = FALSE)
axis(1, cex.axis = 1.25)
axis(2, at = 0:20, cex.axis = 1.25, las = 2)
title(xlab = expression("Observations"), cex.lab = 1.5)
title(ylab = expression("Y"), line = 3.75, cex.lab = 1.5)
abline(h = mean(Y), lwd = 2)
abline(h = var(Y), lwd = 2, lty = 2)
box()
legend("topleft", bty = "n", c(expression("Mean of Y"),
                               expression("Variance of Y")), lty = c(1, 2), lwd = 2, cex = 1.5)

dev.off()

# Poisson vs. Negative Binomial
library(MASS)
set.seed(763759) # Set the seed for reproducible results

reps <- 1000 # Set the number of repetitions at the top of the script
par.est.pnb <- matrix(NA, nrow = reps, ncol = 4) # Empty matrix to store the
# estimates 
b0 <- .2 # True value for the intercept
b1 <- .5 # True value for the slope
n <- 1000 # Sample size
X <- runif(n, -1, 1) # Create a sample of n observations on the 
# independent variable X

for(i in 1:reps){
  Y <- rnbinom(n, size = .5, mu = exp(b0 + b1*X)) # Generate data with
  # overdispersion
  model.p <- glm(Y ~ X, family = "poisson") # Estimate Poisson model
  model.nb <- glm.nb(Y ~ X) # Estimate NB model
  vcv.p <- vcov(model.p) # Variance-covariance matrices
  vcv.nb <- vcov(model.nb)
  par.est.pnb[i, 1] <- model.p$coef[2]  # Store the results
  par.est.pnb[i, 2] <- model.nb$coef[2] 
  par.est.pnb[i, 3] <- sqrt(diag(vcv.p)[2])    
  par.est.pnb[i, 4] <- sqrt(diag(vcv.nb)[2])
  cat("Completed", i, "of", reps, "\n")
}

# Means of the coefficient on X estimates
mean(par.est.pnb[ , 1]) # Poisson estimates
mean(par.est.pnb[ , 2]) # NB estimates

# MSE
mean((par.est.pnb[ , 1])^2) # Poisson MSE
mean((par.est.pnb[ , 2])^2) # NB MSE

# Coverage probabilities 
# Poisson SEs
coverage(par.est.pnb[ , 1], par.est.pnb[ , 3], b1,
         df = n - model.p$rank)$coverage.probability
# NB SEs
coverage(par.est.pnb[ , 2], par.est.pnb[ , 4], b1,
         df = n - model.nb$rank)$coverage.probability

pdf("od2.pdf")

par(mar = c(5, 5.25, .5, .5)) 
plot(1:n, Y, ylim = c(0, 20), col = "gray50", xlab = "", ylab = "",
     main = "", axes = FALSE)
axis(1, cex.axis = 1.25)
axis(2, at = 0:20, cex.axis = 1.25, las = 2)
title(xlab = expression("Observations"), cex.lab = 1.5)
title(ylab = expression("Y"), line = 3.75, cex.lab = 1.5)
abline(h = mean(Y), lwd = 2)
abline(h = var(Y), lwd = 2, lty = 2)
box()
legend("topleft", bty = "n", c(expression("Mean of Y"),
                               expression("Variance of Y")), lty = c(1, 2), lwd = 2, cex = 1.5)

dev.off()

# Negative Binomial vs. Zero-Inflated Negative Binomial
library(pscl)
# Zero-inflated negative binomial random number generator
rzinbinom <- function(n, mu, size, zprob){ 
  ifelse(rbinom(n, 1, zprob) == 1, 0, rnbinom(n, size = size, mu = mu))
}

set.seed(2837) # Set the seed for reproducible results

reps <- 1000 # Set the number of repetitions at the top of the script
par.est.zinb <- matrix(NA, nrow = reps, ncol = 4) # Empty matrix to store
# the estimates 
b0z <- -.8 # True value for the inflation intercept
b1z <- .3 # True value for the inflation slope                                                  
b0c <- .2 # True value for the count intercept
b1c <- .5 # True value for the count slope
n <- 1000 # Sample size
X <- runif(n, -1, 1) # Create a sample of n observations on the 
# independent variable X
Z <- rnorm(n, X, 1) # Inflation independent variable

for(i in 1:reps){
  # Generate data with a zero-inflation component
  Y.zi <- rzinbinom(n, mu = exp(b0c + b1c*X), size = .5,
                    zprob = exp(b0z + b1z*Z)/(1 + exp(b0z + b1z*Z)))
  
  # Generate data with no zero-inflation component
  Y.nozi <- rzinbinom(n, mu = exp(b0c + b1c*X), size = .5, zprob = 0) 
  model.nb1 <- glm.nb(Y.zi ~ X) # Standard negative binomial
  model.nb2 <- glm.nb(Y.nozi ~ X)
  model.zinb1 <- zeroinfl(Y.zi ~ X | Z, dist = "negbin") # Zero-inflated model
  model.zinb2 <- zeroinfl(Y.nozi ~ X | Z, dist = "negbin") 
  
  # Store the estimates of the coefficient on X (count equation)
  par.est.zinb[i, 1] <- model.nb1$coef[2] # Standard NB, with ZI
  par.est.zinb[i, 2] <- model.nb2$coef[2] # Standard NB, no ZI
  par.est.zinb[i, 3] <- as.numeric(model.zinb1$coef$count[2]) # ZI NB, with ZI
  par.est.zinb[i, 4] <- as.numeric(model.zinb2$coef$count[2]) # ZI NB, no ZI
  cat("Completed", i, "of", reps, "\n")
}

pdf("zinb-zi.pdf")

par(mar = c(5, 5.25, .5, .5)) 
plot(density(par.est.zinb[ , 1]), xlim = c(-.5, 1.5), ylim = c(0, 3.5),
     lwd = 3, xlab = "", ylab = "", main = "", axes = FALSE)
lines(density(par.est.zinb[ , 3]), lwd = 3, lty = 2)
axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression(hat(beta[1])), cex.lab = 1.5)
title(ylab = expression("Density"), line = 3.75, cex.lab = 1.5)
abline(v = b1, lwd = 2)
text(1, 2, expression("True"~beta[1]~"= 0.50"), cex = 1.5)
box()
legend("topright", bty = "n", c(expression("Standard NB"),
                                expression("ZINB")), lty = c(1, 2), lwd = 3, cex = 1.5)

dev.off()

pdf("zinb-nozi.pdf")

par(mar = c(5, 5.25, .5, .5)) 
plot(density(par.est.zinb[ , 2]), xlim = c(-.5, 1.5), ylim = c(0, 5),
     lwd = 3, xlab = "", ylab = "", main = "", axes = FALSE)
lines(density(par.est.zinb[ , 4]), lwd = 3, lty = 2)
axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression(hat(beta[1])), cex.lab = 1.5)
title(ylab = expression("Density"), line = 3.75, cex.lab = 1.5)
abline(v = b1, lwd = 2)
text(1, 2, expression("True"~beta[1]~"= 0.50"), cex = 1.5)
box()
legend("topright", bty = "n", c(expression("Standard NB"),
                                expression("ZINB")), lty = c(1, 2), lwd = 3, cex = 1.5)

dev.off()
