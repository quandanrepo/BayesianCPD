# Reduce the problem to two models within the one changepoint model:
# i.e. take tau being in two places and we can incorporate the changepoint prior as model prior
# Set up problem as before but we will only work with just unknown mean, known sigma

rm(list = ls())
set.seed(3)

# Data --------------------------------------------------------------------
# Currently use simulated data
cp <- 5
y <- c(rnorm(cp,0,1),rnorm(cp,1,1)) ; plot(y)


# Models to test ----------------------------------------------------------

num_models <- 2 #Just test two locations
tau_locations <- c(4,5) #Tau at 4 and 5

# Prior for parameters ----------------------------------------------------
mu_prior <- c(0,3) #Normal distribution 
known_sigma <- 1
model_prior <- rep(1,num_models)/num_models

# MCMC parameters ---------------------------------------------------------
proposal_sd_rj <- 0.3 #bigger jump for the reversible jump
proposal_sd_wm <- 0.1 #within model sample of parameters
sims <- 100000 
approx_samples <- 100000 #Number of samples for monte carlo estimate
perc_burnin <- 0.1
thin_every <- 5 # Add thinning as there might be autocorrelation since alot of proposals will get rejected due to out of bounds proposals for tau
proposal_tau <- function(tau){
  if (tau==4) {
    return(5)
  } else {
    return(4)
  }
}

# Set up parameters to store-----
num_of_cps_list <- rep(NA,sims)
mu_list <- list() #to assign to this list, use [[]]
tau_list <- list()


# Initialisation------
tau <- sample(4:5,1) #just run Gibbs on two models
mu <- mean(y) 
sigma <- sqrt(var(y))

# Analytical Solution ---------------------------------------------
marg_like=rep(NA,num_models)
# Write function for Marginal likelihood
marg_like_mu <- function(y,sigma,mu_0,sigma_0){
  return((sigma/(((sqrt(2*pi)*sigma)^length(y))*sqrt(length(y)*sigma_0^2+sigma^2)))*exp(-(sum(y^2)/(2*sigma^2))-(mu_0^2/(2*sigma_0^2)))*exp(((sigma_0^2*length(y)^2*mean(y)^2/sigma^2)+((sigma^2*mu_0^2)/sigma_0^2)+(2*length(y)*mean(y)*mu_0))/(2*(length(y)*sigma_0^2+sigma^2))))              
}
for (i in 1:num_models) {
  marg_like[i] <- marg_like_mu(y[1:tau_locations[i]],known_sigma,mu_prior[1],mu_prior[2])*marg_like_mu(y[(tau_locations[i]+1):(length(y))],known_sigma,mu_prior[1],mu_prior[2])
}
marg_like <- marg_like/sum(marg_like)

plot(marg_like)
lines(marg_like)
par(mfrow=c(2,1)) #have 2 plots in same space

results <- data.frame("Analytical"=marg_like)


# Approximation using Monte Carlo -----------------------------------------
# Define function for product of normals
prod_norm <- function(mu,sigma,y){prod(dnorm(y,mu,sigma))}

# Use approximation, i.e integrate over all mu
# Marginal likelihood P(Y|M_i)
pM1 <- rep(NA,num_models)
for (i in 1:num_models) {
  pM1[i] <- mean(sapply(rnorm(approx_samples,mu_prior[1],mu_prior[2]),prod_norm,known_sigma,y[1:tau_locations[i]]))*mean(sapply(rnorm(approx_samples,mu_prior[1],mu_prior[2]),prod_norm,known_sigma,y[(tau_locations[i]+1):(2*cp)]))
}
# Normalise
total <- sum(pM1)
pM1 <- pM1/total
# Store results in data frame
results <- cbind(results,"Approx"=c(pM1[1],pM1[2]))

# Gibbs + Metropolis Hasting ----------------------------------------------
# Gibbs Sampler on parameters + Metropolis Hasting on tau 
mu_posterior_dist <- function(y,known_sigma,mu_prior){
  sigma2_posterior <- 1/((1/(mu_prior[2]^2))+(length(y)/(known_sigma^2)))
  mu_posterior <- ((sum(y)/(known_sigma^2))+(mu_prior[1]/(mu_prior[2]^2)))*sigma2_posterior
  return(rnorm(1,mu_posterior,sqrt(sigma2_posterior)))
}
for (s in 1:sims) {
  if (s%%10000==0) {
    print(s)
  }
  #Resample parameters in segment 1 using gibbs sampler
  # Normal likelihood with normal prior
  new_mu_left <- mu_posterior_dist(y[1:tau],known_sigma,mu_prior)
  
  #Resample parameters in segment 2 using gibbs sampler
  # Normal likelihood with normal prior
  new_mu_right <- mu_posterior_dist(y[(tau+1):length(y)],known_sigma,mu_prior)
  
  #Resample tau using Metropolis hasting
  #Propose a new tau usings a symmetric proposal so no transition ratio
  pro_tau <- proposal_tau(tau)

  old_like <- sum(dnorm(y[1:tau],new_mu_left,known_sigma,log = TRUE)) + sum(dnorm(y[(tau+1):length(y)],new_mu_right,known_sigma,log = TRUE))
  new_like <- sum(dnorm(y[1:pro_tau],new_mu_left,known_sigma,log = TRUE)) + sum(dnorm(y[(pro_tau+1):length(y)],new_mu_right,known_sigma,log = TRUE))

  if (runif(1)<exp(new_like-old_like)) {
    tau <- pro_tau
  }

  tau_list[[s]] <- tau
  mu_list[[s]] <- c(new_mu_left,new_mu_right)
}
tau_list <- tau_list[-(1:(perc_burnin*sims))]
tau_freq <- as.data.frame(table(factor(unlist(tau_list[seq(1,length(tau_list),thin_every)]))))
tau_den <- tau_freq[,2]/sum(tau_freq[,2])

plot(marg_like)
lines(marg_like)
plot(tau_den)
lines(tau_den)

results <- cbind(results,"Gibbs+MH"=c(tau_den))

