# Metropolis Hasting for fixed model sampling on tau - compare with analytical and approximation
# Currently only dealing with either no changepoint or 1 changepoint at the correct position
# Assume that we know we are in a fixed model, say one changepoint.
# We derive the analytical solution using normal and use a Metropolis Hasting sampler to sample from the posterior distribution.
# Gibbs sampler on the parameters

# A note must be made regarding the proposals we use for changepoint
# Previously used a poisson changepoint jump on tau but this seems to make the algorithm oversample specific areas instead of traversing across all of the posterior space
# We need a more effective way to proposal changepoint apart from using a uniform proposal, which gives better results than poisson proposal but will probably do a bad job for a large number of data points
# try using tau_lambda as as.integer(n/2) since this will give the algorithm a mean jump of half the length of the data. i.e. from the middle location, you can jump to the ends

rm(list = ls())

set.seed(1)

# Data = "data/coal.csv"

# Currently use simulated data

cp <- 10

y <- c(rnorm(cp,1,1),rnorm(cp,1,2))
# y <- read.csv(Data)

par(mfrow=c(4,1)) #have 2 plots in same space
plot(y)

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

# proposal_sd_rj <- 0.3 #bigger jump for the reversible jump
# proposal_sd_wm <- 0.1 #within model sample of parameters
tau_lambda <- c(1,length(y)/2) # lambda for poisson jump 
large_jump <- 0.01 #prob of doing a large jump
sims <- 50000 
approx_samples <- 100 #Number of samples for monte carlo estimate
perc_burnin <- 0.1
thin_every <- 5 # Add thinning as there might be autocorrelation since alot of proposals will get rejected due to out of bounds proposals for tau
proposal_tau <- function(tau,tau_lambda,large_jump){
  if (runif(1)<large_jump) {
    return(tau+(((-1)^(sample(0:1,1)))*rpois(1,tau_lambda[2])))
  } else{
    return(tau+(((-1)^(sample(0:1,1)))*rpois(1,tau_lambda[1])))
  }
  # return(sample(1:(length(y)-1),1))
}

# Set up parameters to store
num_of_cps_list <- rep(NA,sims)
tau_list <- list()
segment_1 <- list() #to hold both mu and sigma
segment_2 <- list() #to hold both mu and sigma

# Initialisation
tau <- sample(1:(length(y)-1),1)
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
  
  results <- data.frame("Analytical"=marg_like)
  
  # Approximation using Monte Carlo -----------------------------------------
  # Define function for product of normals
  prod_norm <- function(mu,sigma,y){prod(dnorm(y,mu,sigma))}
  
  # Use approximation, i.e integrate over all mu
  # Marginal likelihood P(Y|M_i)
  pM1 <- rep(NA,length(y)-1)
  for (i in 1:(length(y)-1)) {
    pM1[i] <- mean(sapply(rnorm(approx_samples,mu_prior[1],mu_prior[2]),prod_norm,known_sigma,y[1:i]))*mean(sapply(rnorm(approx_samples,mu_prior[1],mu_prior[2]),prod_norm,known_sigma,y[(i+1):(length(y))]))
  }
  # Normalise
  total <- sum(pM1)
  pM1 <- pM1/total
  
  plot(pM1)
  lines(pM1)
  # Store results in data frame
  results <- cbind(results,"Approx"=c(pM1))
  
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
    pro_tau <- proposal_tau(tau,tau_lambda,large_jump)
    
    while (pro_tau<1|pro_tau>(length(y)-1)) {
      pro_tau <- proposal_tau(tau,tau_lambda,large_jump)
    }
    
    old_like <- sum(dnorm(y[1:tau],segment_1[[s]],known_sigma,log = TRUE)) + sum(dnorm(y[(tau+1):length(y)],segment_2[[s]],known_sigma,log = TRUE))
    new_like <- sum(dnorm(y[1:pro_tau],segment_1[[s]],known_sigma,log = TRUE)) + sum(dnorm(y[(pro_tau+1):length(y)],segment_2[[s]],known_sigma,log = TRUE))
    
    old_prior <- tau_prior(tau,y,log = TRUE)
    new_prior <- tau_prior(pro_tau,y,log = TRUE)
    
    if (runif(1)<exp(new_like+new_prior-old_like-old_prior)) {
      tau <- pro_tau
    }
    
    tau_list[[s]] <- tau
  }
  tau_list <- tau_list[-(1:(perc_burnin*sims))]
  tau_freq <- as.data.frame(table(factor(unlist(tau_list[seq(1,length(tau_list),thin_every)]))))
  tau_den <- tau_freq[,2]/sum(tau_freq[,2])
  
  # plot(marg_like)
  # lines(marg_like)
  plot(tau_den)
  lines(tau_den)
  
  results <- cbind(results,"Gibbs+MH"=c(tau_den))
  
} else if (cpd_test==2) {#Known mu, unknown sigma^2
  # Use approximation using samples from prior, reversible jump, and analytical
  known_mu <- as.numeric(readline(prompt = "Known mu: "))
  
  # Before anything, define function for product of normals
  prod_norm <- function(sigma,mu,y){prod(dnorm(y,mu,sigma))}
  
  # Define function to generate and evaluate Scaled-Inverse-Chi-Squared variables, parameters (scale,degree_of_freedom)
  rscaleinvchi <- function(n,df_para,scale_para){
    return((df_para*scale_para)/rchisq(n,df_para))
  }
  dscaleinvchi <- function(x,df_para,scale_para,log=c(TRUE,FALSE)){
    if (log) {
      return(log((scale_para*df_para/2)^(df_para/2)/(gamma(df_para/2))*exp((-df_para*scale_para)/(2*x))/(x^(1+(df_para/2)))))
    } else {
      return((scale_para*df_para/2)^(df_para/2)/(gamma(df_para/2))*exp((-df_para*scale_para)/(2*x))/(x^(1+(df_para/2))))
    } 
  }
  
  # Do analytical solution first
  marg_like=rep(NA,(length(y)-1))
  
  # Write function for Marginal likelihood
  marg_like_sigma <- function(y,mu,df_0,scale_0){
    return(((1/sqrt(2*pi))^(length(y)))*(((scale_0*df_0)/2)^(df_0/2))/gamma(df_0/2)*gamma((df_0+length(y))/2)*((2/((df_0*scale_0)+(sum((y-mu)^2))))^((df_0+length(y))/2)))
  }
  for (i in 1:(length(y)-1)) {
    marg_like[i] <- marg_like_sigma(y[1:i],known_mu,sigma_prior[1],sigma_prior[2])*marg_like_sigma(y[(i+1):(length(y))],known_mu,sigma_prior[1],sigma_prior[2])
  }
  # marg_like <- marg_like[-1]
  marg_like <- marg_like/sum(marg_like)
  
  plot(marg_like)
  lines(marg_like)

  results <- data.frame("Analytical"=marg_like)
  
  # Approximation using Monte Carlo -----------------------------------------
  # Define function for product of normals
  prod_norm <- function(sigma,mu,y){prod(dnorm(y,mu,sigma))}
  
  # Use approximation, i.e integrate over all mu
  # Marginal likelihood P(Y|M_i)
  pM1 <- rep(NA,length(y)-1)
  for (i in 1:(length(y)-1)) {
    pM1[i] <- mean(sapply(rnorm(approx_samples,sigma_prior[1],sigma_prior[2]),prod_norm,known_mu,y[1:i]))*mean(sapply(rnorm(approx_samples,sigma_prior[1],sigma_prior[2]),prod_norm,known_mu,y[(i+1):(length(y))]))
  }
  # Normalise
  total <- sum(pM1)
  pM1 <- pM1/total
  
  plot(pM1)
  lines(pM1)
  # Store results in data frame
  results <- cbind(results,"Approx"=c(pM1))
  
  # Gibbs Sampler on parameters + Metropolis Hasting on tau -------------------
  sigma_posterior_dist <- function(y,known_mu,sigma_prior){
    df_n <- sigma_prior[1]+length(y)
    scale2_n <-(1/df_n)*(sum((y-known_mu)^2)+(sigma_prior[1]*sigma_prior[2]))
    return(rscaleinvchi(1,df_n,scale2_n))
  }
  for (s in 1:sims) {
    if (s%%10000==0) {
      print(s)
    }
    #Resample parameters in segment 1 using gibbs sampler
    # Normal likelihood with normal prior
    
    segment_1[[s]] <- sigma_posterior_dist(y[1:tau],known_mu,sigma_prior)
    
    #Resample parameters in segment 2 using gibbs sampler
    # Normal likelihood with normal prior
    
    segment_2[[s]] <- sigma_posterior_dist(y[(tau+1):length(y)],known_mu,sigma_prior)
    
    #Resample tau using Metropolis hasting
    #Propose a new tau usings a symmetric proposal so no transition ratio
    pro_tau <- proposal_tau(tau,tau_lambda,large_jump)
    
    while (pro_tau<1|pro_tau>(length(y)-1)) {
      pro_tau <- proposal_tau(tau,tau_lambda,large_jump)
    }
    
    old_like <- sum(dnorm(y[1:tau],known_mu,segment_1[[s]],log = TRUE)) + sum(dnorm(y[(tau+1):length(y)],known_mu,segment_2[[s]],log = TRUE))
    new_like <- sum(dnorm(y[1:pro_tau],known_mu,segment_1[[s]],log = TRUE)) + sum(dnorm(y[(pro_tau+1):length(y)],known_mu,segment_2[[s]],log = TRUE))
    
    old_prior <- tau_prior(tau,y,log = TRUE)
    new_prior <- tau_prior(pro_tau,y,log = TRUE)
    
    if (runif(1)<exp(new_like+new_prior-old_like-old_prior)) {
      tau <- pro_tau
    }
    
    tau_list[[s]] <- tau
  }
  tau_list <- tau_list[-(1:(perc_burnin*sims))]
  tau_freq <- as.data.frame(table(factor(unlist(tau_list[seq(1,length(tau_list),thin_every)]))))
  tau_den <- tau_freq[,2]/sum(tau_freq[,2])
  
  # plot(marg_like)
  # lines(marg_like)
  plot(tau_den)
  lines(tau_den)
  
  results <- cbind(results,"Gibbs+MH"=c(tau_den))
  
} else {#Unknown Mu, Unknown Sigma^2
  # Do analytical solution first
  marg_like=rep(NA,(length(y)-1))
  
  # Before anything, define function for product of normals
  prod_norm <- function(mu,sigma,y){prod(dnorm(y,mu,sigma))}
  
  # Define function to generate and evaluate Scaled-Inverse-Chi-Squared variables, parameters (scale,degree_of_freedom)
  rscaleinvchi <- function(n,df_para,scale_para){
    return((df_para*scale_para)/rchisq(n,df_para))
  }
  dscaleinvchi <- function(x,df_para,scale_para,log=c(TRUE,FALSE)){
    if (log) {
      return(log((scale_para*df_para/2)^(df_para/2)/(gamma(df_para/2))*exp((-df_para*scale_para)/(2*x))/(x^(1+(df_para/2)))))
    } else {
      return((scale_para*df_para/2)^(df_para/2)/(gamma(df_para/2))*exp((-df_para*scale_para)/(2*x))/(x^(1+(df_para/2))))
    } 
  }
  
  # Write function for Marginal likelihood
  marg_like_mu_sigma <- function(y,mu_prior,sigma_prior){
    df_n <- sigma_prior[1]+length(y)
    kappa_n <- mu_prior[2]+length(y)
    mu_n <- (mu_prior[2]*mu_prior[1]+length(y)*mean(y))/kappa_n
    scale_n <- (1/df_n)*(sigma_prior[1]*sigma_prior[2]+sum((y-mean(y))^2)+((length(y)*mu_prior[2])/(mu_prior[2]+length(y)))*(mu_prior[1]-mean(y))^2)
    return((1/(pi^(length(y)/2)))*(sqrt(mu_prior[2]/kappa_n))*(gamma(df_n/2)/gamma(sigma_prior[1]/2))*((sigma_prior[1]*sigma_prior[2])^(sigma_prior[1]/2))/(df_n*scale_n)^(df_n/2))
  }
  
  # Do analytical solution first
  marg_like=rep(NA,(length(y)-1))
  
  for (i in 1:(length(y)-1)) {
    marg_like[i] <- marg_like_mu_sigma(y[1:i],mu_prior,sigma_prior)*marg_like_mu_sigma(y[(i+1):(length(y))],mu_prior,sigma_prior)
  }
  # marg_like <- marg_like[-1]
  marg_like <- marg_like/sum(marg_like)
  
  plot(marg_like)
  lines(marg_like)

  results <- data.frame("Analytical"=marg_like)
  

# Approximation Monte Carlo -----------------------------------------------
# Will be alot tricker to do a approximation since we need to do
# 1. Generate approx_sample sigma, and for each sigma, generate approx sample mu
# 2. Calculate the integration of both segments by taking mean

  # Marginal Likelihood P(Y|M_1)
  test_sigma_M1_1 <- sqrt(rscaleinvchi(approx_samples,sigma_prior[1],sigma_prior[2]))
  test_mu_M1_1 <- list()
  for (i in 1:length(test_sigma_M1_1)) {
    test_mu_M1_1[[i]] <- rnorm(approx_samples,mu_prior[1],sqrt(test_sigma_M1_1[i]^2/mu_prior[2]))
  }
  test_sigma_M1_2 <- sqrt(rscaleinvchi(approx_samples,sigma_prior[1],sigma_prior[2]))
  test_mu_M1_2 <- list()
  for (i in 1:length(test_sigma_M1_2)) {
    test_mu_M1_2[[i]] <- rnorm(approx_samples,mu_prior[1],sqrt(test_sigma_M1_2[i]^2/mu_prior[2]))
  }
  pM1 <- list()
  for (j in 1:(length(y)-1)) { #number of changepoint models
    temp_approx <- rep(NA,approx_samples) #initialise the storage
    for (i in 1:approx_samples) { # integrate over both parameters
      temp_approx[i] <- log(mean(sapply(test_mu_M1_1[[i]],prod_norm,test_sigma_M1_1[i],y[1:j]))*mean(sapply(test_mu_M1_2[[i]],prod_norm,test_sigma_M1_2[i],y[(j+1):(length(y))])))
    }
    pM1[j] <- mean(exp(temp_approx)) #need to take exp
  }
  approx_posterior <- unlist(pM1)
  # Normalise
  # total <- pM0+pM1
  # pM0 <- pM0/total #Not accurate due to the number of samples required for a good approximation
  approx_posterior <- approx_posterior/sum(approx_posterior)#
  
  plot(approx_posterior)
  lines(approx_posterior)
  # Store results in data frame
  results <- cbind(results,"Approx"=approx_posterior)


# Gibbs Sampler + MH on tau ------------------------------------------------------
  
  sigma_posterior_dist <- function(y,mu_prior,sigma_prior){
    df_n <- sigma_prior[1]+length(y)
    scale2_n <- (1/df_n)*(sigma_prior[1]*sigma_prior[2]+sum((y-mean(y))^2)+((length(y)*mu_prior[2])/(mu_prior[2]+length(y)))*(mu_prior[1]-mean(y))^2)
    return(rscaleinvchi(1,df_n,scale2_n))
  }
  
  mu_posterior_dist <- function(y,mu_prior,sigma){
    kappa_n <- mu_prior[2]+length(y)
    mu_n <- (mu_prior[2]*mu_prior[1]+length(y)*mean(y))/kappa_n
    return(rnorm(1,mu_n,sqrt((sigma^2)/(kappa_n))))
  }
  
  for (s in 1:sims) {
    if (s%%10000==0) {
      print(s)
    }
    #Resample parameters in segment 1 using gibbs sampler
    #sample sigma2 first, then mu conditional on sigma2
    temp_sigma2 <- sigma_posterior_dist(y[1:tau],mu_prior,sigma_prior)
    temp_mu <- mu_posterior_dist(y[1:tau],mu_prior,sqrt(temp_sigma2))
    segment_1[[s]] <- c(temp_mu,sqrt(temp_sigma2))
    
    #Resample parameters in segment 2 using gibbs sampler
    #sample sigma2 first, then mu conditional on sigma2
    temp_sigma2 <- sigma_posterior_dist(y[(tau+1):length(y)],mu_prior,sigma_prior)
    temp_mu <- mu_posterior_dist(y[(tau+1):length(y)],mu_prior,sqrt(temp_sigma2))
    segment_2[[s]] <- c(temp_mu,sqrt(temp_sigma2))
    
    #Resample tau using Metropolis hasting
    #Propose a new tau usings a symmetric proposal so no transition ratio
    pro_tau <- proposal_tau(tau,tau_lambda,large_jump)
    
    while (pro_tau<1|pro_tau>(length(y)-1)) {
      pro_tau <- proposal_tau(tau,tau_lambda,large_jump)
    }
    
    old_like <- sum(dnorm(y[1:tau],segment_1[[s]][1],segment_1[[s]][2],log = TRUE)) + sum(dnorm(y[(tau+1):length(y)],segment_2[[s]][1],segment_2[[s]][2],log = TRUE))
    new_like <- sum(dnorm(y[1:pro_tau],segment_1[[s]][1],segment_1[[s]][2],log = TRUE)) + sum(dnorm(y[(pro_tau+1):length(y)],segment_2[[s]][1],segment_2[[s]][2],log = TRUE))
    
    old_prior <- tau_prior(tau,y,log = TRUE)
    new_prior <- tau_prior(pro_tau,y,log = TRUE)
    
    if (runif(1)<exp(new_like+new_prior-old_like-old_prior)) {
      tau <- pro_tau
    }
    
    tau_list[[s]] <- tau
  }
  tau_list <- tau_list[-(1:(perc_burnin*sims))]
  tau_freq <- as.data.frame(table(factor(unlist(tau_list[seq(1,length(tau_list),thin_every)]))))
  tau_den <- tau_freq[,2]/sum(tau_freq[,2])

  # plot(marg_like)
  # lines(marg_like)
  plot(tau_den)
  lines(tau_den)
  
  results <- cbind(results,"Gibbs+MH"=c(tau_den))
  
  
}