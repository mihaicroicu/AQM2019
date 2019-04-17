## Model refers to a model object, for instance an 'lm' or 'glm' object
## Data refers to the dataframe
## The function returns the bootstrapped coefficients from the model

bootstrap_function <- function(model,data){
  df_new <- data[sample(1:nrow(data),size=nrow(data),replace=T),]
  m2 <- update(model,data=df_new)
  if(class(model)=="zeroinfl")
    coefs <-as.data.frame(t(m2$coefficients$count))
  else{
    coefs <-as.data.frame(t(m2$coef))
  }
  return(coefs)
}