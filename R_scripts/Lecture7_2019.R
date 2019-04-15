################################################################################
# Script for lecture 7, Advanced Quantitative Methods, Uppsala University      #
# H??vard Hegre, April 2019.                                                    #
# Based on Monte Carlo Simulation and Resampling Methods for Social Science,   #
# Thomas M. Carsey and Jeffrey J. Harden, Chapter 6 File                       #
# Last update: 4/14/19                                                         #
################################################################################

# Set the working directory
setwd("C:/Users/jjharden/Dropbox/Research/Sim Book")
setwd("C:/Users/hegha560.000/Dropbox/Undervisning/UU/AdvancedQuantitativeMethods/WD")
setwd("C:/Users/hhegre/Dropbox/Undervisning/UU/AdvancedQuantitativeMethods/WD")
#setwd("D:/Users/peaceful/Desktop/WD")

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

# Some example poisson distributions
Poi_2 <- rpois(10000, 2) 
Poi_10 <- rpois(10000, 10) 
hist(Poi_2,breaks=50)
describe(Poi_2)
summary(Poi_2)
var(Poi_2)
hist(Poi_10,breaks=50)
describe(Poi_10)
summary(Poi_10)
var(Poi_10)
# Negative binomial
nb_size = 10000
nb_mu = 10
Nbin_10 <- rnbinom(10000, size = nb_size, mu = nb_mu)
# Variance is mu + mu^2/size
hist(Nbin_10,breaks=50)
describe(Nbin_10)
summary(Nbin_10)
var(Nbin_10)
nb_mu + nb_mu^2/nb_size


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



# Clustered Data
library(mvtnorm)
library(lme4)
# Function to compute robust cluster standard errors (Arai 2011)
rcse <- function(model, cluster){
  require(sandwich) 
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- model$rank
  dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
  uj  <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum))
  rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
  return(rcse.cov)
}

set.seed(28704) # Set the seed for reproducible results

reps <- 1000 # Set the number of repetitions at the top of the script
par.est.cluster <- matrix(NA, nrow = reps, ncol = 3) # Empty matrix to store
# the estimates 
se.est.cluster <- matrix(NA, nrow = reps, ncol = 4) # Empty matrix to store
# the standard errors  
b0 <- .2 # True value for the intercept
b1 <- .5 # True value for the slope
n <- 1000 # Sample size
p <- .5 # Rho; specifying the variances of effects 1 and 2
nc <- 50 # Number of clusters
c.label <- rep(1:nc, each = n/nc) # Cluster label

for(i in 1:reps){ # Start the loop
  i.sigma <- matrix(c(1, 0, 0, 1 - p), ncol = 2) # Level 1 effects; correlated as specified by (1-p)
  i.values <- rmvnorm(n = n, sigma = i.sigma)
  effect1 <- i.values[ , 1] # individual-level term going into X
  effect2 <- i.values[ , 2] # individual-level term going into the error term
  
  c.sigma <- matrix(c(1, 0, 0, p), ncol = 2) # Level 2 effects; correlated as specified by p
  c.values <- rmvnorm(n = nc, sigma = c.sigma)
  effect3 <- rep(c.values[ , 1], each = n/nc) # group-level effect going into X
  effect4 <- rep(c.values[ , 2], each = n/nc) # group-level effect going into the error term
  
  X <- effect1 + effect3 # X values unique to level 1 observations
  error <- effect2 + effect4 # 
  
  Y <- b0 + b1*X + error # True model
  
  model.ols <- lm(Y ~ X) # Model estimation 
  model.fe <- lm(Y ~ X + factor(c.label))
  model.mlm <- lmer(Y ~ X + (1|c.label))
  
  par.est.cluster[i, 1] <- model.ols$coef[2] # Coefficients
  par.est.cluster[i, 2] <- model.fe$coef[2]
  par.est.cluster[i, 3] <- fixef(model.mlm)[2]
  
  vcv.ols <- vcov(model.ols) # Variance-covariance matrices
  vcv.rcse <- rcse(model.ols, c.label)
  vcv.fe <- vcov(model.fe)
  vcv.mlm <- vcov(model.mlm)
  
  se.est.cluster[i, 1] <- sqrt(diag(vcv.ols)[2]) # Standard errors
  se.est.cluster[i, 2] <- sqrt(diag(vcv.rcse)[2])
  se.est.cluster[i, 3] <- sqrt(diag(vcv.fe)[2])
  se.est.cluster[i, 4] <- sqrt(diag(vcv.mlm)[2])
} # End the loop

# Coverage probabilities
ols.cp <- coverage(par.est.cluster[ , 1], se.est.cluster[ , 1], b1,
                   df = n - model.ols$rank)
rcse.cp <- coverage(par.est.cluster[ , 1], se.est.cluster[ , 2], b1,
                    df = n - model.ols$rank)
fe.cp <- coverage(par.est.cluster[ , 2], se.est.cluster[ , 3], b1,
                  df = n - model.fe$rank)
mlm.cp <- coverage(par.est.cluster[ , 3], se.est.cluster[ , 4], b1,
                   df = n - length(fixef(model.mlm)))

pdf("clustered-density.pdf")

par(mar = c(5, 5.25, .5, .5))
plot(density(par.est.cluster[ , 1]), lty = 1, xlim = c(.2, .8),
     ylim = c(0, 20), lwd = 3, xlab = "", ylab = "", main = "", axes = FALSE)
lines(density(par.est.cluster[ , 2]), lwd = 3, lty = 2)
lines(density(par.est.cluster[ , 3]), lwd = 3, lty = 3)
axis(1, at = seq(0, 1, .1), cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression(hat(beta[1])), cex.lab = 1.5)
title(ylab = expression("Density"), line = 3.75, cex.lab = 1.5)
abline(v = b1, lwd = 2)
text(.7, 7, expression("True"~beta[1]~"= 0.50"), cex = 1.5)
box()
legend("topright", bty = "n", c(expression("OLS"), expression("OLS with FE"),
                                expression("MLM")), lty = c(1, 2, 3), lwd = 3, cex = 1.5)

dev.off()

pdf("clustered-cp.pdf")

par(mar = c(5, 5.25, .5, .5))
plot(1, ols.cp$coverage.probability, pch = 19, xlim = c(0, 8),
     ylim = c(.5, 1), lwd = 3, xlab = "", ylab = "", main = "", axes = FALSE)
segments(1, ols.cp$mc.eb[1], 1, ols.cp$mc.eb[2], lwd = 2)
points(3, rcse.cp$coverage.probability, pch = 1, lwd = 3)
segments(3, rcse.cp$mc.eb[1], 3, rcse.cp$mc.eb[2], lwd = 2) 
points(5, fe.cp$coverage.probability, pch = 19, lwd = 3) 
segments(5, fe.cp$mc.eb[1], 5, fe.cp$mc.eb[2], lwd = 2) 
points(7, mlm.cp$coverage.probability, pch = 19, lwd = 3) 
segments(7, mlm.cp$mc.eb[1], 7, mlm.cp$mc.eb[2], lwd = 2) 
axis(1, at = c(1, 3, 5, 7), labels = c(expression("OLS"), expression("RCSE"),
                                       expression("OLS with FE"), expression("MLM")), cex.axis = 1.25)
axis(2, at = seq(.5, 1, .05), cex.axis = 1.25, las = 2)
title(xlab = expression("Estimator"), cex.lab = 1.5)
title(ylab = expression("Coverage Probability"), line = 3.75, cex.lab = 1.5)
abline(h = .95, lwd = 2, lty = 2)
box()

dev.off()

# Alternative multi-level DGP (HH addition)

set.seed(28704) # Set the seed for reproducible results

reps <- 1000 # Set the number of repetitions at the top of the script
par.est.cluster <- matrix(NA, nrow = reps, ncol = 3) # Empty matrix to store
# the estimates 
se.est.cluster <- matrix(NA, nrow = reps, ncol = 4) # Empty matrix to store
# the standard errors  
b0 <- .2 # True value for the intercept
b1 <- .5 # True value for the slope
a1 <- 1 # Term setting how much  alpha_j influences X: i.e. correlation between them
n <- 1000 # Sample size
nc <- 100 # Number of clusters
alpha_j <- rnorm(nc,0,1) # Clusters
alpha_ji <- rep(alpha_j[], each = n/nc) # group-level effect going into X
c.label <- rep(1:nc, each = n/nc) # Cluster label
X <- rnorm(n,a1*alpha_j,1) # X is drawn from a normal distribution with mean specified by a1*alpha_j; inducing correlation between X and groups

describe(cluster)
head(cluster)

for(i in 1:reps){ # Start the loop
  Y <- alpha_ji + b0 + b1*X + rnorm(n,0,1) # True model
  
  model.ols <- lm(Y ~ X) # Model estimation 
  model.fe <- lm(Y ~ X + factor(c.label))
  model.mlm <- lmer(Y ~ X + (1|c.label))
  
  par.est.cluster[i, 1] <- model.ols$coef[2] # Coefficients
  par.est.cluster[i, 2] <- model.fe$coef[2]
  par.est.cluster[i, 3] <- fixef(model.mlm)[2]
  
  vcv.ols <- vcov(model.ols) # Variance-covariance matrices
  vcv.rcse <- rcse(model.ols, c.label)
  vcv.fe <- vcov(model.fe)
  vcv.mlm <- vcov(model.mlm)
  
  se.est.cluster[i, 1] <- sqrt(diag(vcv.ols)[2]) # Standard errors
  se.est.cluster[i, 2] <- sqrt(diag(vcv.rcse)[2])
  se.est.cluster[i, 3] <- sqrt(diag(vcv.fe)[2])
  se.est.cluster[i, 4] <- sqrt(diag(vcv.mlm)[2])
} # End the loop


# Coverage probabilities
ols.cp <- coverage(par.est.cluster[ , 1], se.est.cluster[ , 1], b1,
                   df = n - model.ols$rank)
rcse.cp <- coverage(par.est.cluster[ , 1], se.est.cluster[ , 2], b1,
                    df = n - model.ols$rank)
fe.cp <- coverage(par.est.cluster[ , 2], se.est.cluster[ , 3], b1,
                  df = n - model.fe$rank)
mlm.cp <- coverage(par.est.cluster[ , 3], se.est.cluster[ , 4], b1,
                   df = n - length(fixef(model.mlm)))

pdf("clustered-density.pdf")

par(mar = c(5, 5.25, .5, .5))
plot(density(par.est.cluster[ , 1]), lty = 1, xlim = c(.2, .8),
     ylim = c(0, 20), lwd = 3, xlab = "", ylab = "", main = "", axes = FALSE)
lines(density(par.est.cluster[ , 2]), lwd = 3, lty = 2)
lines(density(par.est.cluster[ , 3]), lwd = 3, lty = 3)
axis(1, at = seq(0, 1, .1), cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression(hat(beta[1])), cex.lab = 1.5)
title(ylab = expression("Density"), line = 3.75, cex.lab = 1.5)
abline(v = b1, lwd = 2)
text(.7, 7, expression("True"~beta[1]~"= 0.50"), cex = 1.5)
box()
legend("topright", bty = "n", c(expression("OLS"), expression("OLS with FE"),
                                expression("MLM")), lty = c(1, 2, 3), lwd = 3, cex = 1.5)

dev.off()

pdf("clustered-cp.pdf")

par(mar = c(5, 5.25, .5, .5))
plot(1, ols.cp$coverage.probability, pch = 19, xlim = c(0, 8),
     ylim = c(.5, 1), lwd = 3, xlab = "", ylab = "", main = "", axes = FALSE)
segments(1, ols.cp$mc.eb[1], 1, ols.cp$mc.eb[2], lwd = 2)
points(3, rcse.cp$coverage.probability, pch = 1, lwd = 3)
segments(3, rcse.cp$mc.eb[1], 3, rcse.cp$mc.eb[2], lwd = 2) 
points(5, fe.cp$coverage.probability, pch = 19, lwd = 3) 
segments(5, fe.cp$mc.eb[1], 5, fe.cp$mc.eb[2], lwd = 2) 
points(7, mlm.cp$coverage.probability, pch = 19, lwd = 3) 
segments(7, mlm.cp$mc.eb[1], 7, mlm.cp$mc.eb[2], lwd = 2) 
axis(1, at = c(1, 3, 5, 7), labels = c(expression("OLS"), expression("RCSE"),
                                       expression("OLS with FE"), expression("MLM")), cex.axis = 1.25)
axis(2, at = seq(.5, 1, .05), cex.axis = 1.25, las = 2)
title(xlab = expression("Estimator"), cex.lab = 1.5)
title(ylab = expression("Coverage Probability"), line = 3.75, cex.lab = 1.5)
abline(h = .95, lwd = 2, lty = 2)
box()

dev.off()

summary(model.ols)
summary(model.fe)
summary(model.mlm)

