# Simple reversible jump model
# Extension of Green 2009 paper: Data modelled either by Normal or Location-Scale-t

library(LaplacesDemon) #for student-t distribution

rm(list = ls())
# Set seed
set.seed(1)

# Input data
# Decide to use either Normal or Location-Scale-t 
y=rnorm(1000,0,1) #Normal data
# y=rst(100,0,1,5) #Student-t data
plot(y)

# Set priors
# Assume equal chance of being in either model i.e. Prior set to be equal -> p(k=1)=p(k=2)=1/2
model_prior <- c(1/2,1/2)
mu_prior <- c(0,1) #prior parameters for mu
sigma_prior <- c(1) #prior parameters for sigma #Assume that sigma is distribution by half-normal distribution
df_prior <- c(3,30) #prior parameters for df

# Need to write a function to evaluate half-normal distribution
dhalfnorm <- function(y,scale_halfnorm,log){
  if(log=="TRUE"){
    den_halfnorm <- log(2*dnorm(y,0,scale_halfnorm))
  }
  else if(log=="FALSE"){
    den_halfnorm <- 2*dnorm(y,0,scale_halfnorm)
  }
  else {
    den_halfnorm <- 2*dnorm(y,0,scale_halfnorm)
  }
  return(den_halfnorm)
}

# proposal parameters
proposal_df <- c(3,30) #hyperparameter choices for df

# Begin sampler in k=1 i.e. normal model
sims <- 100000 #Iterations
mu_list <- rep(NA,sims) #Keep mu information
sigma_list <- rep(NA,sims) #Keep sigma information
df_list <- rep(NA,sims) #Keep df information

#Proposalsd for Metropolis Hastings
mu_proposal_mh <- c(0,1)
sigma_proposal_mh <- c(0,1)
df_proposal_mh <- c(0,1)
               
model_list <- rep(NA,sims) #Keep model information
para <- c(0,1) #Initialise mu and sigma
model <- 1 #Start in Normal model
for (s in 1:sims) {
  if(s%%10000==0){
    print(s)
  }
  if(model==1){#We propose jump to Student-t
    new_para <- c(para[1],para[2],runif(1,proposal_df[1],proposal_df[2]))
    
    old_like <- sum(dnorm(y,para[1],para[2],log = TRUE))
    new_like <- sum(dst(y,new_para[1],new_para[2],new_para[3],log=TRUE))
    old_prior <- dnorm(para[1],mu_prior[1],mu_prior[2],log = TRUE) + dhalfnorm(para[2],sigma_prior,log=TRUE)  
    new_prior <- dnorm(new_para[1],mu_prior[1],mu_prior[2],log = TRUE) + dunif(new_para[3],df_prior[1],df_prior[2],log = TRUE) +dhalfnorm(new_para[2],sigma_prior,log = TRUE)
      
    old_to_new <- dunif(new_para[3],proposal_df[1],proposal_df[2],log = TRUE)
    new_to_old <- log(1) #To reverse jump, we have a deterministic proposal
    
    if(runif(1)<exp(new_like+new_prior-old_like-old_prior+new_to_old-old_to_new)){
      para <- new_para
      model <- 2
    }
  }
  else{#We propose jump to normal
    new_para <- c(para[1],para[2])
    
    old_like <- sum(dst(y,para[1],para[2],para[3],log=TRUE))
    new_like <- sum(dnorm(y,new_para[1],new_para[2],log = TRUE))
    old_prior <- dnorm(para[1],mu_prior[1],mu_prior[2],log = TRUE) + dunif(para[3],df_prior[1],df_prior[2],log = TRUE) +dhalfnorm(para[2],sigma_prior,log = TRUE)
    new_prior <- dnorm(new_para[1],mu_prior[1],mu_prior[2],log = TRUE) +dhalfnorm(new_para[2],sigma_prior,log = TRUE)
      
    old_to_new <- log(1) #To reverse jump, we have a deterministic proposal
    new_to_old <- dunif(para[3],proposal_df[1],proposal_df[2],log = TRUE)
    
    if(runif(1)<exp(new_like+new_prior-old_like-old_prior+new_to_old-old_to_new)){
      para <- new_para
      model <- 1
    
  }
  }
  #resample parameters for either model
  #metropolis hasting sampling
  if(model==1){
    new_para <- c(para[1]+rnorm(1,mu_proposal_mh[1],mu_proposal_mh[2]),para[2])
    
    old_like <- sum(dnorm(y,para[1],para[2],log = TRUE))
    new_like <- sum(dnorm(y,new_para[1],new_para[2],log = TRUE))
    
    old_prior <- dnorm(para[1],mu_prior[1],mu_prior[2],log = TRUE) + dhalfnorm(para[2],sigma_prior,log = TRUE)
    new_prior <- dnorm(new_para[1],mu_prior[1],mu_prior[2],log = TRUE) + dhalfnorm(new_para[2],sigma_prior,log = TRUE)
    
    if(runif(1)<exp(new_like+new_prior-old_like-old_prior)){
      para <- new_para
    }
  }
  else {
    new_para <- c(para[1]+rnorm(1,mu_proposal_mh[1],mu_proposal_mh[2]),para[2],para[3]+rnorm(1,df_proposal_mh[1],df_proposal_mh[2]))
    
    if(new_para[3]<3){
      new_prior <- -Inf
      new_like <- -Inf
    }
    else{
      new_prior <- dnorm(new_para[1],mu_prior[1],mu_prior[2],log = TRUE) + dunif(new_para[3],df_prior[1],df_prior[2],log = TRUE) + dhalfnorm(new_para[2],sigma_prior,log = TRUE)
      new_like <- sum(dst(y,new_para[1],new_para[2],new_para[3],log = TRUE))
    }
    
    old_like <- sum(dst(y,para[1],para[2],para[3],log = TRUE))
    
    old_prior <- dnorm(para[1],mu_prior[1],mu_prior[2],log = TRUE) + dunif(para[3],df_prior[1],df_prior[2],log = TRUE) + dhalfnorm(para[2],sigma_prior,log = TRUE)
    
    if(runif(1)<exp(new_like+new_prior-old_like-old_prior)){
      para <- new_para
    }
  }
  
  model_list[s] <- model
  mu_list[s] <- para[1]
  sigma_list[s] <- para[2]
  df_list[s] <- para[3]
}

plot(factor(model_list))
hist(df_list,breaks = 300)