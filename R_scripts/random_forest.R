 rm(list=ls())
library(ranger)
library(yardstick)
library(tidyverse)
library(pdp)

set.seed(1)

# AQM 2019 random forest example
# This example uses the iris dataset, loaded by default.
# This attempts to classify Iris flowers in three categories (virginica, versicolor, setosa)

# (C) Mihai Croicu 2019-04-25.

# Iris is a default R dataset based on Fisher's 1934 book. We'll copy it into flowers.
flowers <- iris
str(flowers)
head(flowers)

# Add a row number to the dataset to create our test and train partitions for cross-validation
flowers <- flowers %>% mutate(rowid = row_number())

# Sample 80% for training
train <- flowers %>% sample_frac(.80)

# And get the rest using anti_join - a join that returns what was NOT joined
test  <- flowers %>% anti_join(train, by = 'rowid')

# Remember to remove the row_id here!
train <- train %>% select(-rowid)
test <- test %>% select(-rowid)

#Fit our random forest. Default num.trees is 500, and min.node.size=10
#This can be good or bad for many tasks.
#Remember, if you have a classification (dichotomous or multinomial)
#You need to specify probability = T here to get probabilities. If you do not need probabilities, just 1/0 or classes, set to F
#And importance = 'permutation', or, if too slow, change to 'impurity'. If you do not need it, set to 'none'
rf_fit <- ranger(Species ~ ., data = train, 
                 num.trees = 2500, 
                 min.node.size = 1,
                 probability = TRUE,
                 importance = 'permutation')
rf_fit
summary(fit)

#Get the variable importance scores / feature importance scores (VIS/FIS)
rf_importance <- importance(rf_fit)
rf_importance

#Predict the test dataset from our training dataset
predict_obj <- predict(rf_fit, data=test)

#This returns an object containing probabilities. 
#If we want to inspect predictions we need to extract the data_frame
#Like this:
predict_df <- data.frame(predict_obj[['predictions']])

### ROC and PR curves. Here, the example is multiclass.
actuals <- test%>%select(Species)%>%rename(actual=Species)
full_predictions <- bind_cols(predict_df,actuals)

roc_data <- roc_curve(full_predictions,truth = actual, setosa, virginica, versicolor)
autoplot(roc_data)
roc_auc(full_predictions,truth = actual, setosa, virginica, versicolor)

pr_data <- pr_curve(full_predictions,truth = actual, setosa, virginica, versicolor)
autoplot(pr_data)
pr_auc(full_predictions,truth = actual, setosa, virginica, versicolor)


#### Partial dependency plots (PDP) (aka QI/Marginal Effects) plots:
# Compute PDPs.
# You can do it in the Carsey and Harden way, but there are tools for these for RF.
# If you do it the Carsey and Harden way, you don't sample from the covariance matrix. 
# Instead, you bootstrap the fitted trees a number of times (1..1000 times) 
# Making a new, bootstrapped random forest
# if you want uncertainty estimates.
# However they give us tools to do it automagically:

# For ?

#If you don't set the prob to T, you get a scale based on a centered logit, 
#i.e. how much the effect of the variable you choose changes at each level versus a level of 0

#Since both tidyverse and pdp contain a partial function
#We EXPLICITLY state what we want saying we want it from the pdp:: package

variable_of_interest = 'Petal.Width'
pdp_rf<-pdp::partial(rf_fit, pred.var=variable_of_interest, prob=T)
pdp_rf
plotPartial(pdp_rf)

# For multi-class classification, (multinomial dependent variable) this produces
# a prediction for the BASE CLASS (usually the first label, the factor with the lowest score)
# and the other classes together. Change the base class to see marginal effects for each class.
# e.g. in our case:
# setosa = 1 
# versicolor = 2    
# virginica = 3
# if we set variable_of_interest = 'Petal.Width' and which.class = 3 (virginica)
# it computes the probability of 
# petal width being being Iris Virginica across all petal widths. 
# holding all other predictors at mean.

variable_of_interest = 'Petal.Width'
pdp_rf<-pdp::partial(rf_fit, pred.var=variable_of_interest, prob=T, which.class=3)
pdp_rf
# Probability of a plant being Iris Virginica (3) given petal.width = 0...2.5
plotPartial(pdp_rf)
# Do this in a loop for all predicted classes!

# GOOD LUCK!
# /Mihai
