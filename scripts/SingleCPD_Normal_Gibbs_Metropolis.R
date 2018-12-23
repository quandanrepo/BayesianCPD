# Metropolis Hasting for fixed model sampling on tau - compare with analytical and approximation
# Currently only dealing with either no changepoint or 1 changepoint at the correct position
# Assume that we know we are in a fixed model, say one changepoint.
# We derive the analytical solution using normal and use a Metropolis Hasting sampler to sample from the posterior distribution.
# Gibbs sampler on the parameters

rm(list = ls())

set.seed(3)

Data = ""

# Currently use simulated data

cp <- 5

y <- c(rnorm(cp,0,1),rnorm(cp,1,1)) ; plot(y)

#Use the same three methodolodies as before

#Unknown mean, known sigma
#Known mean, unknown sigma
#Unknown mean, Unknown sigma

# prior for parameter
mu_prior <- c(0,0.3) #Normal distribution #Mu is conditional on sigma i.e sigma^2/kappa #unless we choose unknown mean and known sigma (1)
sigma_prior <- c(5,0.6) #Use Scaled-Inverse-Chi-Squared for Conjugacy (df_para,scale_para) - mean of c(5,0.6) = 1
tau_prior <- function(tau,y,log=c(TRUE,FALSE)){ #Just assume uniform for now
  if (tau<1|tau>(length(y)-1)){
    return(0)
  } else {
    if (log) {
      return(log(1/(length(y)-1)))
    } else {
      return(1/(length(y)-1))
    }
  }
}

proposal_sd_rj <- 0.3 #bigger jump for the reversible jump
proposal_sd_wm <- 0.1 #within model sample of parameters
tau_lambda <- 2 # lambda for poisson jump 
sims <- 100000 
# approx_samples <- 1000 #Number of samples for monte carlo estimate, although we probably dont need it for this method
perc_burnin <- 0.1
thin_every <- 5 # Add thinning as there might be autocorrelation since alot of proposals will get rejected due to out of bounds proposals for tau
proposal_tau <- function(tau,tau_lambda){
  # return(tau+(((-1)^(sample(0:1,1)))*rpois(1,tau_lambda)))
  # return(sample(1:9,1))
  if (tau==4) {
    return(5)
  } else {
    return(4)
  }
}

# Set up parameters to store
num_of_cps_list <- rep(NA,sims)
# mu_list <- list() #to assign to this list, use [[]]
# sigma_list <- list() #to assign to this list, use [[]]
tau_list <- list()
segment_1 <- list()
segment_2 <- list()

# Initialisation
# tau <- sample(1:(length(y)-1),1) 
tau <- sample(4:5,1) #just run Gibbs on two models
mu <- mean(y)
sigma <- sqrt(var(y))

cpd_test <- as.numeric(readline(prompt = "Unknown Mu(1), Unknown Sigma(2), Unknown Mu and Sigma(3): ")) # Start with Unknown Mean, Known Sigma

if (cpd_test==1){ #Unknown Mu, known sigma^2
  known_sigma <- as.numeric(readline(prompt = "Known sigma: "))
  mu_prior[1] <- as.numeric(readline(prompt = "Prior Mean: "))
  mu_prior[2] <- as.numeric(readline(prompt = "Prior std: "))
  
  # Do analytical solution first
  marg_like=rep(NA,(length(y)-1))
  # Write function for Marginal likelihood
  marg_like_mu <- function(y,sigma,mu_0,sigma_0){
    return((sigma/(((sqrt(2*pi)*sigma)^length(y))*sqrt(length(y)*sigma_0^2+sigma^2)))*exp(-(sum(y^2)/(2*sigma^2))-(mu_0^2/(2*sigma_0^2)))*exp(((sigma_0^2*length(y)^2*mean(y)^2/sigma^2)+((sigma^2*mu_0^2)/sigma_0^2)+(2*length(y)*mean(y)*mu_0))/(2*(length(y)*sigma_0^2+sigma^2))))              
  }
  for (i in 1:(length(y)-1)) {
    marg_like[i] <- marg_like_mu(y[1:i],known_sigma,mu_prior[1],mu_prior[2])*marg_like_mu(y[(i+1):(length(y))],known_sigma,mu_prior[1],mu_prior[2])
  }
  # marg_like <- marg_like[-1]
  marg_like <- marg_like/sum(marg_like)
  
  plot(marg_like)
  lines(marg_like)
  par(mfrow=c(2,1)) #have 2 plots in same space
  
  results <- data.frame("Analytical"=marg_like)
  
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

    segment_1[[s]] <- mu_posterior_dist(y[1:tau],known_sigma,mu_prior)
    
    #Resample parameters in segment 2 using gibbs sampler
    # Normal likelihood with normal prior
    
    segment_2[[s]] <- mu_posterior_dist(y[(tau+1):length(y)],known_sigma,mu_prior)
    
    #Resample tau using Metropolis hasting
    #Propose a new tau usings a symmetric proposal so no transition ratio
    pro_tau <- proposal_tau(tau,tau_lambda)
    
    while (pro_tau<1|pro_tau>(length(y)-1)) {
      pro_tau <- proposal_tau(tau,tau_lambda)
    }
    
    old_like <- sum(dnorm(y[1:tau],segment_1[[s]],known_sigma,log = TRUE)) + sum(dnorm(y[(tau+1):length(y)],segment_2[[s]],known_sigma,log = TRUE))
    new_like <- sum(dnorm(y[1:pro_tau],segment_1[[s]],known_sigma,log = TRUE)) + sum(dnorm(y[(pro_tau+1):length(y)],segment_2[[s]],known_sigma,log = TRUE))
    
    old_prior <- tau_prior(tau,y,log = TRUE)
    new_prior <- tau_prior(pro_tau,y,log = TRUE)
    
    if (runif(1)<new_like+new_prior-old_like-old_prior) {
      tau <- pro_tau
    }
    
    tau_list[[s]] <- tau
  }
  tau_list <- tau_list[-(1:(perc_burnin*sims))]
  tau_freq <- as.data.frame(table(factor(unlist(tau_list[seq(1,length(tau_list),thin_every)]))))
  tau_den <- tau_freq[,2]/sum(tau_freq[,2])
  
  plot(marg_like)
  lines(marg_like)
  plot(tau_den)
  lines(tau_den)
  
  results <- cbind(results,"Gibbs+MH"=c(tau_freq[,2]))
  
} else if (cpd_test==2) {#Known mu, unknown sigma^2
  
} else {#Unknown Mu, Unknown Sigma^2
  
}





