
#Load packages
library(mvtnorm)
library(tidyverse)
library(sandwich)
## Inverse logit function
inv.logit <- function(Z){
  return(exp(Z)/(1 + exp(Z)))
}

## Bootstrapping
#Set seed
set.seed(3141)
#Set sample size
samp_size <- 1000

#Generate some independent data with some multicolinearity
xmat <- rmvnorm(samp_size,c(0,1),matrix(c(1,0.2,0.2,2),2,2,byrow=T))
x3 <- rpois(samp_size,apply(abs(xmat),1,sum))

# Make data.frame
df <- data.frame(x1 = xmat[,1],
                 x2 = xmat[,2],
                 x3 = x3)

# Create odds and probabilities
df <- df %>% mutate(y = 1.5*x1 + 0.7*x2 - 0.5*x3 +rnorm(samp_size,0,3))


# Estimate model
m1 <- lm(y~x1+x2+x3,data=df)
summary(m1)

##
b3_se <- summary(m1)$coef[4,2]
b3_est <- summary(m1)$coef[4,1]
b3_sim <- rnorm(1000,b3_est,b3_se)

## Setup bootstrap loop
b3_bs <- rep(NA,1000)
for(i in 1:1000){
  bs_index <- sample(1:nrow(df),replace=T)
  df_bs <- df[bs_index,]
  m1_bs <- lm(y~x1+x2+x3,data=df_bs)
  b3_bs[i] <- m1_bs$coefficients[4]
}

df_density <- data.frame(regular = b3_sim,
                          bootstrap = b3_bs) %>% gather(key="type",value="value")

ggplot(df_density,aes(x=value,col=type))+
  geom_density()


## OLS with heteroscedasticity
#Generate some independent data with some multicolinearity

# Create dependent variable with complex heteroskedasticity
df_het <- df %>% mutate(y_het = 1.5*x1 + 0.7*x2 - 0.5*x3 +rnorm(samp_size,0,apply(abs(df[,1:3]),1,sum)))


# Estimate model
m1_het <- lm(y_het~x1+x2+x3,data=df_het)
summary(m1_het)

##
b3_se_het <- summary(m1_het)$coef[4,2]
b3_est_het <- summary(m1_het)$coef[4,1]
b3_sim_het <- rnorm(1000,b3_est_het,b3_se_het)

## Setup bootstrap loop
b3_bs_het <- rep(NA,1000)
for(i in 1:1000){
  bs_index <- sample(1:nrow(df),replace=T)
  df_bs_het <- df_het[bs_index,]
  m1_bs_het <- lm(y_het~x1+x2+x3,data=df_bs_het)
  b3_bs_het[i] <- m1_bs_het$coefficients[4]
}

df_density_het <- data.frame(regular = b3_sim_het,
                             bootstrap = b3_bs_het) %>% gather(key="type",value="value")

ggplot(df_density_het,aes(x=value,col=type))+
  geom_density()

df_density_het %>% group_by(type) %>% summarize(mean=mean(value),
                                                sd = sd(value))

b3_se_hc <- sqrt(vcovHC(m1_het)[4,4])
b3_sim_hc <- rnorm(1000,b3_est_het,b3_se_hc)

df_density_het <- data.frame(regular = b3_sim_het,
                             bootstrap = b3_bs_het,
                             hc = b3_sim_hc) %>% gather(key="type",value="value")


ggplot(df_density_het,aes(x=value,col=type))+
  geom_density()

df_density_het %>% group_by(type) %>% summarize(mean=mean(value),
                                                sd = sd(value))

### Simulating QI

#Generate some independent data with some multicolinearity
xmat <- rmvnorm(samp.size,c(0,1),matrix(c(1,0.2,0.2,2),2,2,byrow=T))
x3 <- rpois(samp.size,apply(abs(xmat),1,sum))

# Make data.frame
df <- data.frame(x1 = xmat[,1],
                 x2 = xmat[,2],
                 x3 = x3)

# Create odds and probabilities
df <- df %>% mutate(z = 1.5*x1 + 0.7*x2 - 0.5*x3 +rnorm(samp.size,0,3),
                    y = rbinom(samp.size,1,inv.logit(z)))



# Create some objects related to x3 and create a DF with all values of X to simulate over
span_x3 <- seq(from=range(df$x3)[1],to=range(df$x3)[2],by=1)
length_x3 <- length(range(df$x3)[1]:range(df$x3)[2])
df_x_test <- data.frame(x1=mean(df$x1),x2=mean(df$x2),x3=span_x3)


# Extract coefs and vcov of model 1
coefs_m1 <- m1$coefficients
vcov_m1 <- vcov(m1)

#number of sims
n_sims <- 1000

# Draw coefficients
new_coefs <- rmvnorm(n_sims,coefs_m1,vcov_m1)

# Create out matrix
probs <- matrix(NA,n_sims,length_x3)
#simulate over all sims
for(i in 1:n_sims){
    m1_new <- m1 # copy model
    m1_new$coefficients <- new_coefs[i,] # enter new coefs
    probs[i,] <- predict(m1_new,newdata=df_x_test,type="response") # predict for all values of x3 
}
# Make dataframe and name columns
probs <- as.data.frame(probs)
colnames(probs) <- span_x3

#Tidy dataframe
probs_tidy <- gather(probs,key="x3",value="probability")

# Create lower and upper limits and summarize the simulations
lower <- 0.025
higher <- 0.975
df_summary <- probs_tidy %>% group_by(x3) %>% summarize(mean = mean(probability),
                                                        lo = quantile(probability,lower),
                                                        hi = quantile(probability,higher)) %>%   ungroup()

# Make plot
p1 <- ggplot(df_summary,aes(x=as.numeric(x3),y=mean))+
      geom_line() + scale_y_continuous(limits=c(0,1)) +
      labs(y = "probability",x="x3 value")+
      geom_line(aes(x=as.numeric(x3),y=lo),lty=2)+
      geom_line(aes(x=as.numeric(x3),y=hi),lty=2)

p1

### Observed value solution

# Create out matrix
probs_ovs <- matrix(NA,n_sims,length_x3)

#simulate over all sims
for(i in 1:n_sims){
  m1_new <- m1 # copy model
  m1_new$coefficients <- new_coefs[i,] # enter new coefs
  for(j in 1:length(span_x3)){
  df_new <- df %>% mutate(x3 = span_x3[j])
  probs_ovs[i,j] <- mean(predict(m1_new,newdata=df_new,type="response")) # predict for all values of x3 
  }
}

# Make dataframe and name columns
probs_ovs <- as.data.frame(probs_ovs)
colnames(probs_ovs) <- span_x3

#Tidy dataframe
probs_tidy_ovs <- gather(probs_ovs,key="x3",value="probability")

# Create lower and upper limits and summarize the simulations
lower <- 0.025
higher <- 0.975
df_summary_ovs <- probs_tidy_ovs %>% group_by(x3) %>% summarize(mean = mean(probability),
                                                        lo = quantile(probability,lower),
                                                        hi = quantile(probability,higher)) %>%   ungroup()


p2 <- p1 + geom_line(data=df_summary_ovs,aes(x=as.numeric(x3),y=mean),col="red")+
    geom_line(data=df_summary_ovs,aes(x=as.numeric(x3),y=lo),lty=2,col="red")+
  geom_line(data=df_summary_ovs,aes(x=as.numeric(x3),y=hi),lty=2,col="red")

p2
### Bootstrap solution (average case)

# Create out matrix
probs_bs <- matrix(NA,n_sims,length_x3)

#simulate over all sims
for(i in 1:n_sims){
  bs_index <- sample(1:nrow(df),replace=T) # Create bootstrap index
  df_bs <- df[bs_index,] # create bootstrapped data
  m1_bs <- glm(y~x1+x2+x3,family="binomial",data=df_bs) # create bootstrapped model
  probs_bs[i,] <- predict(m1_bs,newdata=df_x_test,type="response") # predict for all values of x3 
}


# Make dataframe and name columns
probs_bs <- as.data.frame(probs_bs)
colnames(probs_bs) <- span_x3

#Tidy dataframe
probs_tidy_bs <- gather(probs_bs,key="x3",value="probability")

# Create lower and upper limits and summarize the simulations
lower <- 0.025
higher <- 0.975
df_summary_bs <- probs_tidy_bs %>% group_by(x3) %>% summarize(mean = mean(probability),
                                                                lo = quantile(probability,lower),
                                                                hi = quantile(probability,higher)) %>%   ungroup()

p3 <- p1 + geom_line(data=df_summary_bs,aes(x=as.numeric(x3),y=mean),col="blue")+
  geom_line(data=df_summary_bs,aes(x=as.numeric(x3),y=lo),lty=2,col="blue")+
  geom_line(data=df_summary_bs,aes(x=as.numeric(x3),y=hi),lty=2,col="blue")

p3

p4 <- p2 + geom_line(data=df_summary_bs,aes(x=as.numeric(x3),y=mean),col="blue")+
  geom_line(data=df_summary_bs,aes(x=as.numeric(x3),y=lo),lty=2,col="blue")+
  geom_line(data=df_summary_bs,aes(x=as.numeric(x3),y=hi),lty=2,col="blue")

p4


### Bootstrap solution (bootstrapped average case)

# Create out matrix
probs_bs2 <- matrix(NA,n_sims,length_x3)

#simulate over all sims
for(i in 1:n_sims){
  bs_index <- sample(1:nrow(df),replace=T) # Create bootstrap index
  df_bs2 <- df[bs_index,] # create bootstrapped data
  m1_bs2 <- glm(y~x1+x2+x3,family="binomial",data=df_bs2) # create bootstrapped model
  df_x_test_bs2 <- data.frame(x1=mean(df_bs2$x1),x2=mean(df_bs2$x2),x3=span_x3) # Bootstrapped average case
  probs_bs2[i,] <- predict(m1_bs2,newdata=df_x_test_bs2,type="response") # predict for all values of x3 
}


# Make dataframe and name columns
probs_bs2 <- as.data.frame(probs_bs2)
colnames(probs_bs2) <- span_x3

#Tidy dataframe
probs_tidy_bs2 <- gather(probs_bs2,key="x3",value="probability")

# Create lower and upper limits and summarize the simulations
lower <- 0.025
higher <- 0.975
df_summary_bs2 <- probs_tidy_bs2 %>% group_by(x3) %>% summarize(mean = mean(probability),
                                                              lo = quantile(probability,lower),
                                                              hi = quantile(probability,higher)) %>%   ungroup()

p5 <- ggplot(df_summary_bs,aes(x=as.numeric(x3),y=mean))+
  geom_line(col="blue") + scale_y_continuous(limits=c(0,1)) +
  labs(y = "probability",x="x3 value")+
  geom_line(aes(x=as.numeric(x3),y=lo),lty=2,col="blue")+
  geom_line(aes(x=as.numeric(x3),y=hi),lty=2,col="blue")

p5 + geom_line(data=df_summary_bs2,aes(x=as.numeric(x3),y=mean),col="orange")+
  geom_line(data=df_summary_bs2,aes(x=as.numeric(x3),y=lo),lty=2,col="orange")+
  geom_line(data=df_summary_bs2,aes(x=as.numeric(x3),y=hi),lty=2,col="orange")








