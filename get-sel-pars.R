# Get selectivity parameters
# Robyn Forrest Oct 21 2022
# Estimate the logistic selectivity parameters given selectivity at age

# FOR SIMPLICITY, RUN THIS ONCE TO OVERWRITE THE DOME-SHAPED SELECTIVITY FROM OPENMSE

library(tidyverse)
library(reshape2)

# data
data <- readRDS("Ex2dataORIG.rda")
dat <- data$BI
obs <- dat$sel
obs[8:length(obs)] <- 1
age <- 1:length(obs)

# parameters
ah <- 3
sig <- 0.2
pars <- c(ah,sig)

# model
get.sel <- function(pars){
  
  #model
  ah=pars[1]
  sig=pars[2]
  
  pred <- 1/(1+exp(-(age-ah)/sig)) # this is the same as pred <- plogis(age,ah,sig)
  
  #likelihood (just a simple sum of squares)
  like <- sum((obs-pred)^2)	
  
  #output
  output=list()
  output$like=like
  output$pred=pred
  return(output)
}

tmp <- get.sel(pars)
tmp$pred

fn=function(pars) {get.sel(pars)$like}

fit = optim(pars,fn, method="BFGS", hessian=F)
fit$par

Predicted <- get.sel(fit$par)$pred

# plot
D <- cbind(age,Predicted, obs) %>%
  as.data.frame() %>% 
  rename(Observed=obs) %>% 
  melt(id="age", value.name="Selectivity")
  
p <- ggplot(D)+
  geom_path(aes(x=age,y=Selectivity,colour=variable,linetype=variable), size=1.5)+
  theme_classic()
print(p)

## Overwrite data$BI$sel - ONLY DO THIS IF YOU WANT TO REPLACE THE SELECTIVITY
## FROM THE OPENMSE OM

# dat$sel <- Predicted
# data$BI <- dat
# saveRDS(data,file="Ex2data.rda")



