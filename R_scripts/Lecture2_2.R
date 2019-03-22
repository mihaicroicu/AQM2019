################################################################################
# Advanced Quantitative Methods Course                                         #
# H??vard Hegre, Mihai Croicu and David Randahl                                #
# Lecture 2 File 2                                                             #
# Experiment exploring violations of homoskedasticity                          #
# Last update: 22 March 2019                                                   #
################################################################################


rm(list=ls(all=TRUE)) # Clear the workspace
# Import packages (remember to install them if they are new to your system)
library(xtable)
library(stargazer)

# Set the working directory

setwd("~/Dropbox/Apps/ShareLatex/Advanced Quantitative Methods/WD")
getwd()


# Define required functions first

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


############################################################################
############################ Heteroskedasticity ############################
############################################################################
# Heteroskedasticity (Simulation 1 of 3)
set.seed(100484) # Set the seed for reproducible results

reps <- 1000 # Set the number of repetitions at the top of the script
par.est.ncv <- matrix(NA, nrow = reps, ncol = 4) # Empty matrix to store the
# estimates
sigma.est <- numeric(reps) # Empty vector to store sigma 
b0 <- .2 # True value for the intercept
b1 <- .5 # True value for the slope
n <- 1000 # Sample size
X <- runif(n, -1, 1) # Create a sample of n observations on the 
# independent variable X
gamma <- 1.5 # Heteroskedasticity parameter



for(i in 1:reps){ # Start the loop
  Y <- b0 + b1*X + rnorm(n, 0, exp(X*gamma)) # Now the error variance is a 
  # function of X plus random noise
  model <- lm(Y ~ X) # Estimate OLS model
  sigma.est[i] <- summary(model)$sigma # Store sigma
  vcv <- vcov(model) # Variance-covariance matrix
  par.est.ncv[i, 1] <- model$coef[1] # Put the estimate for the intercept
  # in the first column
  par.est.ncv[i, 2] <- model$coef[2] # Put the estimate for the coefficient on
  # X in the second column
  par.est.ncv[i, 3] <- sqrt(diag(vcv)[1]) # SE of the intercept
  par.est.ncv[i, 4] <- sqrt(diag(vcv)[2]) # SE of the coefficient on X
} # End the loop

#pdf("ncv-example.pdf")

par(mar = c(5, 5.25, .5, .5))
plot(X, Y, ylim = c(-10, 10), axes = FALSE, xlab = "",
     ylab = "", pch = 19)
title(xlab = expression("X"), cex.lab = 1.5)
title(ylab = expression("Y"), line = 3.75, cex.lab = 1.5)
box()
abline(lsfit(X, Y), lwd = 3)
axis(1, cex.axis = 1.25)
axis(2, at = seq(-10, 10, 2), cex.axis = 1.25, las = 2)

#dev.off()

# Report a table
samplemodel.ncv <- xtable(model)
print.xtable(samplemodel.ncv, type="latex", file="Samplemodel_ncv.tex")
summary(model)

# Plot a typical residual diagnostic plot
#pdf("ncv-diagnostics.pdf") 
plot(X, model$residuals)
#dev.off()



# Heteroskedasticity (Simulation 2 of 3)
# Compare to homoskedasticity with sigma set to the average value
# of the estimates of sigma from the last simulation
sigma <- mean(sigma.est)

set.seed(100484) # Set the seed for reproducible results

reps <- 1000 # Set the number of repetitions at the top of the script
par.est.ncv <- matrix(NA, nrow = reps, ncol = 6) # Empty matrix to store the
# estimates 
b0 <- .2 # True value for the intercept
b1 <- .5 # True value for the slope
n <- 1000 # Sample size
X <- runif(n, -1, 1) # Create a sample of n observations on the 
# independent variable X
gamma <- 1.5 # Heteroskedasticity parameter

for(i in 1:reps){ # Start the loop
  Y1 <- b0 + b1*X + rnorm(n, 0, exp(X*gamma)) # Y1: Heteroskedasticity
  Y2 <- b0 + b1*X + rnorm(n, 0, sigma) # Y2: Homoskedasticity, same average
  # sigma as Y1 
  model1 <- lm(Y1 ~ X) # Estimate OLS models
  model2 <- lm(Y2 ~ X)
  vcv <- vcov(model1) # Variance-covariance matrix (model 1)
  par.est.ncv[i, 1] <- model1$coef[1] # Put the estimate for the intercept
  # in the first column (model 1)
  par.est.ncv[i, 2] <- model1$coef[2] # Put the estimate for the coefficient on
  # X in the second column (model 1)
  par.est.ncv[i, 3] <- model2$coef[1] # Put the estimate for the intercept
  # in the first column (model 2)
  par.est.ncv[i, 4] <- model2$coef[2] # Put the estimate for the coefficient on
  # X in the second column (model 2) 
  par.est.ncv[i, 5] <- sqrt(diag(vcv)[1]) # SE of the intercept (model 1)
  par.est.ncv[i, 6] <- sqrt(diag(vcv)[2]) # SE of the coefficient on X (model 1)
} # End the loop

#pdf("ncv-coef1.pdf")

par(mar = c(5, 5.25, .5, .5))
plot(density(par.est.ncv[ , 1]), lty = 2, ylim = c(0, 12), lwd = 3,
     xlab = "", ylab = "", main = "", axes = FALSE)
lines(density(par.est.ncv[ , 3]), lwd = 3)
axis(1, cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression(hat(beta[0])), cex.lab = 1.5)
title(ylab = expression("Density"), line = 3.75, cex.lab = 1.5)
abline(v = b0, lwd = 2)
text(.1, 8, expression("True"~beta[0]~"= 0.20"), cex = 1.5)
box()
legend("topright", bty = "n", c(expression("Homoskedastic"),
                                expression("Heteroskedastic")), lty = c(1, 2), lwd = 3, cex = 1.5)

#dev.off()

#pdf("ncv-coef2.pdf", width=10, height=6)

par(mar = c(5, 5.25, .5, .5))
plot(density(par.est.ncv[ , 2]), lty = 2, ylim = c(0, 6), lwd = 3,
     xlab = "", ylab = "", main = "", axes = FALSE)
lines(density(par.est.ncv[ , 4]), lwd = 3)
axis(1, at = seq(0, 1, .1), cex.axis = 1.25)
axis(2, cex.axis = 1.25, las = 2)
title(xlab = expression(hat(beta[1])), cex.lab = 1.5)
title(ylab = expression("Density"), line = 3.75, cex.lab = 1.5)
abline(v = b1, lwd = 2)
text(.25, 4, expression("True"~beta[1]~"= 0.50"), cex = 1.5)
box()
legend("topright", bty = "n", c(expression("Homoskedastic"),
                                expression("Heteroskedastic")), lty = c(1, 2), lwd = 3, cex = 1.5)

#dev.off()

# Coverage Probabilities
cp.beta0.ncv <- coverage(par.est.ncv[ , 1], par.est.ncv[ , 5], b0,
                         df = n - model1$rank)
cp.beta0.ncv$coverage.probability
cp.beta0.ncv$mc.eb

cp.beta1.ncv <- coverage(par.est.ncv[ , 2], par.est.ncv[ , 6], b1,
                         df = n - model1$rank)
cp.beta1.ncv$coverage.probability
cp.beta1.ncv$mc.eb


#pdf("l5_cp-plot_heteroskedasticity.pdf", width=10, height=6)
# Changes from C&H original: Use the cp.beta1.ncv and par.est.ncv objects, 
# change the range for beta1 values on y axis

par(mar = c(5, 6, .5, .5))
plot(seq(100, 200, length = 100), seq(0, 1.5, length = 100), type = "n",
     axes = FALSE, xlab = "", ylab = "")
title(xlab = expression("100 Simulated Samples"), cex.lab = 1.5)
title(ylab = expression(hat(beta[1])), line = 3.75, cex.lab = 1.5)
box()
axis(1, at = seq(100, 200, 10), labels = seq(0, 100, 10), cex.axis = 1.25)
axis(2, at = seq(0, 1.5, .1), cex.axis = 1.25, las = 2)
abline(h = b1, lwd = 2)
for (i in 101:200){
  points(i, par.est.ncv[i, 2], lwd = 2, col = ifelse(cp.beta1.ncv$true.in.ci[i] == 1,
                                                     "gray70", "gray20"), pch = 19)
  segments(i, cp.beta1.ncv$ci[i, 1], i, cp.beta1.ncv$ci[i, 2], lwd = 2,
           col = ifelse(cp.beta1.ncv$true.in.ci[i] == 1, "gray70", "gray20"))
}
legend("topleft", bty = "n", c(expression("CI includes true"~beta[1]),
                               expression("CI does not include true"~beta[1])),
       fill = c("gray70", "gray20"), cex = 1.5) 

#dev.off()



# Heteroskedasticity: Robust SEs (Simulation 3 of 3)
# Just simulate the version with heteroskedasticity, assess performance of
# robust standard errors
library(sandwich)

set.seed(100484) # Set the seed for reproducible results

reps <- 1000 # Set the number of repetitions at the top of the script
par.est.ncv <- matrix(NA, nrow = reps, ncol = 4) # Empty matrix to store the
# estimates 
b0 <- .2 # True value for the intercept
b1 <- .5 # True value for the slope
n <- 1000 # Sample size
X <- runif(n, -1, 1) # Create a sample of n observations on the 
# independent variable X
gamma <- 1.5 # Heteroskedasticity parameter

for(i in 1:reps){ # Start the loop
  Y <- b0 + b1*X + rnorm(n, 0, exp(X*gamma)) # Now the error variance is a 
  # function of X plus random noise
  model <- lm(Y ~ X) # Estimate OLS model
  vcv <- vcovHC(model) # Robust variance-covariance matrix
  par.est.ncv[i, 1] <- model$coef[1] # Put the estimate for the intercept
  # in the first column
  par.est.ncv[i, 2] <- model$coef[2] # Put the estimate for the coefficient on
  # X in the second column
  par.est.ncv[i, 3] <- sqrt(diag(vcv)[1]) # SE of the intercept
  par.est.ncv[i, 4] <- sqrt(diag(vcv)[2]) # SE of the coefficient on X
} # End the loop

cp.beta1.ncv.robust <- coverage(par.est.ncv[ , 2], par.est.ncv[ , 4], b1,
                                df = n - model$rank)
cp.beta1.ncv.robust$coverage.probability
cp.beta1.ncv.robust$mc.eb # Simulation error
