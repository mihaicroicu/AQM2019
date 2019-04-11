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
# Ordered Models
library(MASS)
set.seed(8732) # Set the seed for reproducible results


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




reps <- 1000 # Set the number of repetitions at the top of the script
par.est.oprobit <- matrix(NA, nrow = reps, ncol = 2) # Empty matrices to store
taus.oprobit <- matrix(NA, nrow = reps, ncol = 3)    # the estimates 
b0 <- 0 # True value for the intercept
b1 <- .5 # True value for the slope
n <- 1000 # Sample size
X <- rnorm(n, 0, 1) # Create a sample of n observations on the 
# independent variable X

XB <- b0 + b1*X # Systematic component
sd.error <- 1 # SD of the error of the unobserved Y*
# Define the true cutpoints
tau1 <- qnorm(.1, mean = mean(XB), sd = sqrt(var(XB) + sd.error^2)) 
tau2 <- qnorm(.5, mean = mean(XB), sd = sqrt(var(XB) + sd.error^2)) 
tau3 <- qnorm(.9, mean = mean(XB), sd = sqrt(var(XB) + sd.error^2)) 

for(i in 1:reps){ # Start the loop
  Y.star <- rnorm(n, XB, sd.error) # The unobserved Y*
  Y <- rep(NA, n) # Define Y as a vector of NAs with length n
  Y[Y.star < tau1] <- 1 # Set Y equal to a value according to Y.star
  Y[Y.star >= tau1 & Y.star < tau2] <- 2
  Y[Y.star >= tau2 & Y.star < tau3] <- 3
  Y[Y.star >= tau3] <- 4
  # Estimate ordered model
  model <- polr(as.ordered(Y) ~ X, method = "probit", Hess = TRUE) 
  vcv <- vcov(model) # Variance-covariance matrix
  par.est.oprobit[i, 1] <- model$coef[1] # Put the estimate for the coefficient
  # on X in the second column
  par.est.oprobit[i, 2] <- sqrt(diag(vcv)[1]) # SE of the coefficient on X
  taus.oprobit[i, ] <- model$zeta
  cat("Just completed iteration", i, "\n")
} # End the loop

pdf("ystar.pdf", width=5, height=2)

par(mar = c(5, 5.25, .5, .5))
plot(Y.star, Y, xlim = c(-3, 3), xlab = "", ylab = "", main = "",
     axes = FALSE, pch = 19)
axis(1, at = seq(-3, 3, 1), cex.axis = 1.25)
axis(2, at = 1:4, cex.axis = 1.25, las = 2)
title(xlab = expression("Y* (Unobserved/Continuous)"), cex.lab = 1.5)
title(ylab = expression("Y (Observed/Categorical)"), line = 3.75,
      cex.lab = 1.5)
abline(v = c(tau1, tau2, tau3), lwd = 2, lty = 3)
text(tau1 - .25, 3.5, expression(tau[1]), cex = 2)
text(tau2 - .25, 3.5, expression(tau[2]), cex = 2)
text(tau3 - .25, 3.5, expression(tau[3]), cex = 2)
box()

dev.off()

mean(par.est.oprobit[ , 1]) # Mean of coefficient on X estimates
# Compare the actual taus to the means of the tau estimates
data.frame(True = c(tau1, tau2, tau3), Estimated = apply(taus.oprobit, 2, mean))
# Coverage probability for the coefficient on X
coverage(par.est.oprobit[ , 1], par.est.oprobit[ , 2], b1,
         df = n - length(c(coef(model), model$zeta)))$coverage.probability

# Unordered Models
library(Zelig)
library(ZeligChoice)
library(nnet)

#zelig(as.factor(outcome) ~ bureau + democracy + army.1000 + duration + exports + border + lnpop + logincome, model = "mlogit", data = Table1)


set.seed(45262) # Set the seed for reproducible results

reps <- 1000 # Set the number of repetitions at the top of the script
par.est.mnl <- matrix(NA, nrow = reps, ncol = 4) # Empty matrix to store
# the estimates 
b0A <- .2 # True values for the intercepts
b0B <- -.2
b1A <- .5 # True values for the slopes
b1B <- .75
n <- 1000 # Sample size
X <- runif(n, -1, 1) # Create a sample of n observations on the 
# independent variable X
# Compute the probabilities of each outcome based on the DGP
pA <- exp(b0A + b1A*X)/(1 + exp(b0A + b1A*X) + exp(b0B + b1B*X))
pB <- exp(b0B + b1B*X)/(1 + exp(b0A + b1A*X) + exp(b0B + b1B*X))
pC <- 1 - pA - pB

pA[1:8]
pB[1:8]
pC[1:8]

for(i in 1:reps){ # Start the loop
  Y <- rep(NA, n) # Define Y as a vector of NAs with length n
  for(j in 1:n){ # Create the dependent variable in another loop
    Y[j] <- sample(c("A", "B", "C"), 1, replace = TRUE,
                   prob = c(pA[j], pB[j], pC[j]))
  }
  
# Alternative formulation to create dependent variable:
#  Y[j] <- rMultinom(rbind(c(pA[j], pB[j], pC[j])),1) # Multinomial draw over the probabilities for row i. 
  
  
  # Estimate a MNL model
  model <- zelig(as.factor(Y) ~ X, model = "mlogit",
                 data = data.frame(Y, X), cite = FALSE)
  vcv <- vcov(model) # Variance-covariance matrixlog()
  par.est.mnl[i, 1] <- coefficients(model)[3] # Coefficient on X, outcome 1 
  par.est.mnl[i, 2] <- coefficients(model)[4] # Coefficient on X, outcome 2
  par.est.mnl[i, 3] <- sqrt(diag(vcv)[3]) # SE of coefficient on X, outcome 1
  par.est.mnl[i, 4] <- sqrt(diag(vcv)[4]) # SE of coefficient on X, outcome 2
  cat("Just completed iteration", i, "\n")
} # End the loop

# Mean of coefficients on X estimates
mean(par.est.mnl[ , 1]) # Outcome 1
mean(par.est.mnl[ , 2]) # Outcome 2

# Coverage probabilities for the coefficients on X
# Outcome 1
coverage(par.est.mnl[ , 1], par.est.mnl[ , 3], b1A,
         df = n - length(coef(model)))$coverage.probability
# Outcome 2
coverage(par.est.mnl[ , 2], par.est.mnl[ , 4], b1B,
         df = n - length(coef(model)))$coverage.probability

# Ordered vs. MNL
library(Zelig)
set.seed(99999)

reps <- 1000 # Set the number of repetitions at the top of the script
d.pp <- array(NA, c(4, 3, reps)) # Empty array to store 
# simulated change in probabilities

# Ordered logit model DGP
b0 <- 0 # True value for the intercept
b1 <- .5 # True value for the slope
n <- 1000 # Sample size
X <- runif(n, -1, 1) # Create a sample of n observations on the 
# independent variable X

# MNL model DGP
b0A <- .2 # True values for the intercepts
b0B <- -.2
b1A <- .5 # True values for the slopes
b1B <- .75
n <- 1000 # Sample size

# Compute the probabilities of each outcome based on the DGP
pA <- exp(b0A + b1A*X)/(1 + exp(b0A + b1A*X) + exp(b0B + b1B*X))
pB <- exp(b0B + b1B*X)/(1 + exp(b0A + b1A*X) + exp(b0B + b1B*X))
pC <- 1 - pA - pB

for(i in 1:reps){
  # Ordered dependent variable
  Y.star <- rlogis(n, b0 + b1*X, 1) # The unobserved Y*
  # Define the true cutpoints
  tau1 <- quantile(Y.star, .25) 
  tau2 <- quantile(Y.star, .75) 
  Y.o <- rep(NA, n) # Define Y as a vector of NAs with length n
  Y.o[Y.star < tau1] <- 1 # Set Y equal to a value according to Y.star
  Y.o[Y.star >= tau1 & Y.star < tau2] <- 2
  Y.o[Y.star >= tau2] <- 3
  # Ordered data
  o.data <- data.frame(Y.o, X) # Put the data in a data frame
  
  # Unordered dependent variable
  Y.m <- rep(NA, n) # Define Y as a vector of NAs with length n
  for(j in 1:n){ # Create the dependent variable in another loop
    Y.m[j] <- sample(1:3, 1, replace = TRUE, prob = c(pA[j], pB[j], pC[j]))
  }
  # Unordered data
  m.data <- data.frame(Y.m, X) # Put the data in a data frame
  
  # Estimate the models with the ordered dependent variable
  o.correct <- zelig(as.ordered(Y.o) ~ X, model = "ologit",
                     data = o.data, cite = FALSE)
  m.incorrect <- zelig(as.factor(Y.o) ~ X, model = "mlogit",
                       data = o.data, cite = FALSE)
  
  # Estimate the models with the multinomial dependent variable
  m.correct <- zelig(as.factor(Y.m) ~ X, model = "mlogit",
                     data = m.data, cite = FALSE)
  o.incorrect <- zelig(as.ordered(Y.m) ~ X, model = "ologit",
                       data = m.data, cite = FALSE)
  
  # Set X to its minimum and maximum for each model
  x.oc <- setx(o.correct, X = min(X)) # For o.correct
  x.oc1 <- setx(o.correct, X = max(X)) 
  
  x.mi <- setx(m.incorrect, X = min(X)) # For m.incorrect
  x.mi1 <- setx(m.incorrect, X = max(X)) 
  
  x.mc <- setx(m.correct, X = min(X)) # For m.correct
  x.mc1 <- setx(m.correct, X = max(X))
  
  x.oi <- setx(o.incorrect, X = min(X)) # For o.incorrect
  x.oi1 <- setx(o.incorrect, X = max(X)) 
  
  # Compute the change in expected probabilities of falling in each category
  # when moving from the minimum to the maximum of X
  sim.oc <- sim(o.correct, x = x.oc, x1 = x.oc1)$qi$fd
  sim.mi <- sim(m.incorrect, x = x.mi, x1 = x.mi1)$qi$fd
  sim.mc <- sim(m.correct, x = x.mc, x1 = x.mc1)$qi$fd
  sim.oi <- sim(o.incorrect, x = x.oi, x1 = x.oi1)$qi$fd
  d.pp[1, , i] <- apply(sim.oc, 2, mean)
  d.pp[2, , i] <- apply(sim.mi, 2, mean)
  d.pp[3, , i] <- apply(sim.mc, 2, mean)
  d.pp[4, , i] <- apply(sim.oi, 2, mean)
  cat("Just completed iteration", i, "of", reps, "\n")
}

# Compute the average change in probability for each of the four models
dpp.means <- rbind(apply(d.pp[1, , ], 1, mean), apply(d.pp[2, , ], 1, mean), apply(d.pp[3, , ], 1, mean), apply(d.pp[4, , ], 1, mean))

# Compute the SD of the change in probability for each of the four models
dpp.sds <- rbind(apply(d.pp[1, , ], 1, sd), apply(d.pp[2, , ], 1, sd), apply(d.pp[3, , ], 1, sd), apply(d.pp[4, , ], 1, sd))

# LaTeX table
library(xtable)
results <- rbind(dpp.means, dpp.sds)

xtable(results, digits = 4, caption = "Means and Standard Deviations of the Change in expected probability for Each Category", label = "ord-unord-results", align = "lrrr")
