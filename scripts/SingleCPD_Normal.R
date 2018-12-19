# Assume single changepoint model at a single point cp. Using normal distribution 

rm(list = ls())

set.seed(1)

# Load data in
Data=""
# y <- read.csv(Data)

# Or simulated data
cp <- 5
y1 <- rnorm(cp,0,1)
y2 <- rnorm(cp,0,2)

y <- c(y1,y2) ; plot(y)

approx_samples <- 10000

# prior for parameter
mu_prior <- c(0,3) #Normal distribution
sigma_prior <- c(5,0.6) #Use Scaled-Inverse-Chi-Squared for Conjugacy (df_para,scale_para) - mean of c(5,0.6) = 1

proposal_sd_rj <- 0.2 #bigger jump for the reversible jump
proposal_sd_wm <- 0.1 #within model sample of parameters
sims <- 100000
perc_burnin <- 0.1
num_of_cps_list <- rep(NA,sims)
mu_list <- list() #to assign to this list, use [[]]
sigma_list <- list() #to assign to this list, use [[]]
num_of_cps <- 0
mu <- mean(y)
sigma <- sqrt(var(y))
# Divide the problem into three problems
# 1) Unknown Mean, Known Sigma^2
# 2) Known Mean, Unknown Sigma^2
# 3) Unknown Mean, Unknown Sigma^2

cpd_test <- as.numeric(readline(prompt = "Unknown Mu(1), Unknown Sigma(2), Unknown Mu and Sigma(3): ")) # Start with Unknown Mean, Known Sigma

if (cpd_test==1){ #Unknown Mu, known sigma^2
  # Use approximation using samples from prior, reversible jump, and analytical
  known_sigma <- as.numeric(readline(prompt = "Known sigma: "))
  
  # Before anything, define function for product of normals
  prod_norm <- function(mu,sigma,y){prod(dnorm(y,mu,sigma))}
  
  # Use approximation, i.e integrate over all mu
  # Marginal likelihood P(Y|M_0)
  pM0 <- mean(sapply(rnorm(approx_samples,mu_prior[1],mu_prior[2]),prod_norm,known_sigma,y))
  # Marginal Likelihood P(Y|M_1)
  pM1 <- mean(sapply(rnorm(approx_samples,mu_prior[1],mu_prior[2]),prod_norm,known_sigma,y[1:cp]))*mean(sapply(rnorm(approx_samples,mu_prior[1],mu_prior[2]),prod_norm,known_sigma,y[(cp+1):(2*cp)]))
  # Normalise
  total <- pM0+pM1
  pM0 <- pM0/total
  pM1 <- pM1/total
  # Store results in data frame
  results <- data.frame("Approx"=c(pM0,pM1))
  
  # Now calculate using reversible jump
  # proposal_sd_rj <- 0.1 #bigger jump for the reversible jump
  # proposal_sd_wm <- 0.1 #within model sample of parameters
  # sims <- 100000
  # perc_burnin <- 0.1
  # num_of_cps_list <- rep(NA,sims)
  # mu_list <- list() #to assign to this list, use [[]]
  # num_of_cps <- 0
  # mu <- mean(y)
  
  for (s in 1:sims) {
    if (s%%10000==0) {
      print(s)
    }
    # If we have no changepoints, we jump to one changepoint at cp
    if (num_of_cps==0) {
      u1 <- rnorm(1,0,proposal_sd_rj)
      u2 <- rnorm(1,0,proposal_sd_rj)
      
      new_mu_left <- mean(y[1:cp])+u1
      new_mu_right <- mean(y[(cp+1):(2*cp)])+u2
      
      old_prior <- dnorm(mu,mu_prior[1],mu_prior[2],log = TRUE)
      new_prior <- dnorm(new_mu_left,mu_prior[1],mu_prior[2],log = TRUE) + dnorm(new_mu_right,mu_prior[1],mu_prior[2],log = TRUE)
      
      old_like <- sum(dnorm(y,mu,known_sigma,log = TRUE))
      new_like <- sum(dnorm(y[1:cp],new_mu_left,known_sigma,log = TRUE))+sum(dnorm(y[(cp+1):(2*cp)],new_mu_right,known_sigma,log = TRUE))
      
      old_to_new <- dnorm(u1,0,proposal_sd_rj,log = TRUE)+dnorm(u2,0,proposal_sd_rj,log = TRUE)
      new_to_old <- dnorm(mean(y)-mu,0,proposal_sd_rj,log = TRUE)
      if (runif(1)<exp(new_like+new_prior+new_to_old-old_like-old_prior-old_to_new)) {
        num_of_cps <- 1
        mu <- c(new_mu_left,new_mu_right)
      }
    } else { 
      # We are in the one changepoint model
      u <- rnorm(1,0,proposal_sd_rj)
      
      new_mu <- mean(y)+u
      
      old_prior <- dnorm(mu[1],mu_prior[1],mu_prior[2],log = TRUE)+dnorm(mu[2],mu_prior[1],mu_prior[2],log = TRUE)
      new_prior <- dnorm(new_mu,mu_prior[1],mu_prior[2],log = TRUE)
      
      old_like <- sum(dnorm(y[1:cp],mu[1],known_sigma,log = TRUE))+sum(dnorm(y[(cp+1):(2*cp)],mu[2],known_sigma,log = TRUE))
      new_like <- sum(dnorm(y,new_mu,known_sigma,log = TRUE))
      
      old_to_new <- dnorm(u,0,proposal_sd_rj,log = TRUE)
      new_to_old <- dnorm(mu[1]-mean(y[1:cp]),0,proposal_sd_rj,log = TRUE)+dnorm(mu[2]-mean(y[(cp+1):(2*cp)]),0,proposal_sd_rj,log = TRUE)
      
      if (runif(1)< exp(new_like+new_prior-old_like-old_prior+new_to_old-old_to_new)) {
        mu <- new_mu
        num_of_cps <- 0
      }
    }
    # Resample mu parameter from posterior, use Metropolis Hastings
    if (num_of_cps==0) {
      u <- rnorm(1,0,proposal_sd_wm)
      new_mu <- mu+u

      old_prior <- dnorm(mu,mu_prior[1],mu_prior[2],log = TRUE)
      new_prior <- dnorm(new_mu,mu_prior[1],mu_prior[2],log = TRUE)

      old_like <- sum(dnorm(y,mu,known_sigma,log = TRUE))
      new_like <- sum(dnorm(y,new_mu,known_sigma,log = TRUE))
      # Symmetric proposal so no transition probability
      if (runif(1)<exp(new_like+new_prior-old_prior-old_like)) {
        mu <- new_mu
      }
    } else{
      u1 <- rnorm(1,0,proposal_sd_wm)
      u2 <- rnorm(1,0,proposal_sd_wm)
      new_mu_left <- mu[1]+u1
      new_mu_right <- mu[2]+u2

      old_prior <- dnorm(mu[1],mu_prior[1],mu_prior[2],log = TRUE)+dnorm(mu[2],mu_prior[1],mu_prior[2],log = TRUE)
      new_prior <- dnorm(new_mu_left,mu_prior[1],mu_prior[2],log = TRUE)+dnorm(new_mu_right,mu_prior[1],mu_prior[2],log = TRUE)

      old_like <- sum(dnorm(y[1:cp],mu[1],known_sigma,log = TRUE))+sum(dnorm(y[(cp+1):(2*cp)],mu[2],known_sigma,log = TRUE))
      new_like <- sum(dnorm(y[1:cp],new_mu_left,known_sigma,log = TRUE))+sum(dnorm(y[(cp+1):(2*cp)],new_mu_right,known_sigma,log = TRUE))
      # Symmetric proposal so no transition probability
      if (runif(1)<exp(new_like+new_prior-old_prior-old_like)) {
        mu <- c(new_mu_left,new_mu_right)
      }
    }
    #now resample mu from its posterior
    # if (num_of_cps==0) {
    #   posteriorsd <- sqrt( ( 1/mu_prior[2]^2 + length(y) )^(-1) )
    #   posteriormu <- sum(y)/ ( 1/mu_prior[2]^2 + length(y)) 
    #   mu <- rnorm(1,posteriormu,posteriorsd)
    # } else {
    #   y1 <- y[1:cp]; y2 <- y[(cp+1):(2*cp)]
    #   posteriorsd <- sqrt( ( 1/mu_prior[2]^2 + length(y1) )^(-1) )
    #   posteriormu <- sum(y1)/ ( 1/mu_prior[2]^2 + length(y1)) 
    #   mu[1] <- rnorm(1,posteriormu,posteriorsd)
    #   
    #   posteriorsd <- sqrt( ( 1/mu_prior[2]^2 + length(y2) )^(-1) )
    #   posteriormu <- sum(y2)/ ( 1/mu_prior[2]^2 + length(y2)) 
    #   mu[2] <- rnorm(1,posteriormu,posteriorsd)
    # }
    num_of_cps_list[s] <- num_of_cps
    mu_list[[s]] <- mu
  }
  num_of_cps_list <- num_of_cps_list[-(1:(perc_burnin*sims))]
  mu_list <- mu_list[-(1:(perc_burnin*sims))]
  
  cps_0 <- sum(num_of_cps_list==0)
  cps_1 <- sum(num_of_cps_list==1)
  cps_0_perc <- cps_0/(cps_0+cps_1)
  cps_1_perc <- cps_1/(cps_0+cps_1)
  
  results <- cbind(results,"RJMCMC"=c(cps_0_perc,cps_1_perc))
  
  # Analytical results
  # We need to calculate marginal likelihood since p(M0|Y) proportional to p(Y|M0)p(M0) and p(M1|Y) proportional to p(Y|M1)p(M1)
  
  # Write function for Marginal likelihood
  marg_like_mu <- function(y,sigma,mu_0,sigma_0){
    return((sigma/(((sqrt(2*pi)*sigma)^length(y))*sqrt(length(y)*sigma_0^2+sigma^2)))*exp(-(sum(y^2)/(2*sigma^2))-(mu_0^2/(2*sigma_0^2)))*exp(((sigma_0^2*length(y)^2*mean(y)^2/sigma^2)+((sigma^2*mu_0^2)/sigma_0^2)+(2*length(y)*mean(y)*mu_0))/(2*(length(y)*sigma_0^2+sigma^2))))              
  }
  
  # For model M0, no changepoint
  # marg_like_M0 <- (known_sigma/((sqrt(2*pi)*known_sigma)*sqrt(length(y)*mu_prior[2]^2+known_sigma^2)))*exp(-(sum(y^2)/(2*known_sigma^2))-(mu_prior[1]^2/(2*mu_prior[2]^2)))*exp(((mu_prior[2]^2*length(y)^2*mean(y)^2/known_sigma^2)+((known_sigma^2*mu_prior[1]^2)/mu_prior[2]^2)+(2*length(y)*mean(y)*mu_prior[1]))/(2*(length(y)*mu_prior[2]^2+known_sigma^2)))              
  marg_like_M0 <- marg_like_mu(y,known_sigma,mu_prior[1],mu_prior[2])
    
  # For model M1, one changepoint
  marg_like_M1 <- marg_like_mu(y[1:cp],known_sigma,mu_prior[1],mu_prior[2])*marg_like_mu(y[(cp+1):(2*cp)],known_sigma,mu_prior[1],mu_prior[2])
  
  total_marg <- marg_like_M0+marg_like_M1
  
  marg_like_M0 <- marg_like_M0/total_marg
  marg_like_M1 <- marg_like_M1/total_marg
  
  results <- cbind(results,"Analytical"=c(marg_like_M0,marg_like_M1))
  
}else if(cpd_test==2){ #Known Mean, Unknown Sigma^2
  # Use approximation using samples from prior, reversible jump, and analytical
  known_mu <- as.numeric(readline(prompt = "Known mu: "))
  
  # Before anything, define function for product of normals
  prod_norm <- function(sigma,mu,y){prod(dnorm(y,mu,sigma))}
  
  # Use approximation, i.e integrate over all sigma
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
  
  # Marginal likelihood P(Y|M_0) 
  pM0 <- mean(sapply(sqrt(rscaleinvchi(approx_samples,sigma_prior[1],sigma_prior[2])),prod_norm,known_mu,y))
  # Marginal Likelihood P(Y|M_1) 
  pM1 <- mean(sapply(sqrt(rscaleinvchi(approx_samples,sigma_prior[1],sigma_prior[2])),prod_norm,known_mu,y[1:cp]))*mean(sapply(sqrt(rscaleinvchi(approx_samples,sigma_prior[1],sigma_prior[2])),prod_norm,known_mu,y[(cp+1):(2*cp)]))
  # Normalise
  total <- pM0+pM1
  pM0 <- pM0/total
  pM1 <- pM1/total
  # Store results in data frame
  results <- data.frame("Approx"=c(pM0,pM1))
  
  # Now calculate using reversible jump
  # proposal_sd_rj <- 0.1 #bigger jump for the reversible jump
  # proposal_sd_wm <- 0.1 #within model sample of parameters
  # sims <- 100000
  # perc_burnin <- 0.1
  # num_of_cps_list <- rep(NA,sims)
  # sigma_list <- list() #to assign to this list, use [[]]
  # num_of_cps <- 0
  # sigma <- sqrt(var(y)) #initialisation

  for (s in 1:sims) {
    if (s%%10000==0) {
      print(s)
    }
    # If we have no changepoints, we jump to one changepoint at cp
    if (num_of_cps==0) {
      u1 <- rnorm(1,0,proposal_sd_rj)
      u2 <- rnorm(1,0,proposal_sd_rj)

      new_sigma_left <- sqrt(var(y[1:cp])+u1)
      new_sigma_right <- sqrt(var(y[(cp+1):(2*cp)])+u2)

      old_prior <- dscaleinvchi(sigma^2,sigma_prior[1],sigma_prior[2],log = TRUE)
      new_prior <- dscaleinvchi(new_sigma_left^2,sigma_prior[1],sigma_prior[2],log = TRUE) + dscaleinvchi(new_sigma_right^2,sigma_prior[1],sigma_prior[2],log = TRUE)

      old_like <- sum(dnorm(y,known_mu,sigma,log = TRUE))
      new_like <- sum(dnorm(y[1:cp],known_mu,new_sigma_left,log = TRUE))+sum(dnorm(y[(cp+1):(2*cp)],known_mu,new_sigma_right,log = TRUE))

      old_to_new <- dnorm(u1,0,proposal_sd_rj,log = TRUE)+dnorm(u2,0,proposal_sd_rj,log = TRUE)
      new_to_old <- dnorm(sigma^2-var(y),0,proposal_sd_rj,log = TRUE)
      if (runif(1)<exp(new_like+new_prior+new_to_old-old_like-old_prior-old_to_new)) {
        num_of_cps <- 1
        sigma <- c(new_sigma_left,new_sigma_right)
      }
    } else {
      # We are in the one changepoint model
      u <- rnorm(1,0,proposal_sd_rj)

      new_sigma <- sqrt(var(y)+u)

      old_prior <- dscaleinvchi(sigma[1]^2,sigma_prior[1],sigma_prior[2],log = TRUE)+dscaleinvchi(sigma[2]^2,sigma_prior[1],sigma_prior[2],log = TRUE)
      new_prior <- dscaleinvchi(new_sigma^2,sigma_prior[1],sigma_prior[2],log = TRUE)

      old_like <- sum(dnorm(y[1:cp],known_mu,sigma[1],log = TRUE))+sum(dnorm(y[(cp+1):(2*cp)],known_mu,sigma[2],log = TRUE))
      new_like <- sum(dnorm(y,known_mu,new_sigma,log = TRUE))

      old_to_new <- dnorm(u,0,proposal_sd_rj,log = TRUE)
      new_to_old <- dnorm(sigma[1]^2-var(y[1:cp]),0,proposal_sd_rj,log = TRUE)+dnorm(sigma[2]^2-var(y[(cp+1):(2*cp)]),0,proposal_sd_rj,log = TRUE)

      if (runif(1)< exp(new_like+new_prior-old_like-old_prior+new_to_old-old_to_new)) {
        sigma <- new_sigma
        num_of_cps <- 0
      }
    }
    # Resample mu parameter from posterior, use Metropolis Hastings
    if (num_of_cps==0) {
      u <- rnorm(1,0,proposal_sd_wm)
      temp_new_sigma <- sigma^2+u
      if (temp_new_sigma<=0) {
      } else {
        new_sigma <- sqrt(sigma^2+u)

        old_prior <- dscaleinvchi(sigma^2,sigma_prior[1],sigma_prior[2],log = TRUE)
        new_prior <- dscaleinvchi(new_sigma^2,sigma_prior[1],sigma_prior[2],log = TRUE)
        
        old_like <- sum(dnorm(y,known_mu,sigma,log = TRUE))
        new_like <- sum(dnorm(y,known_mu,new_sigma,log = TRUE))
        # Symmetric proposal so no transition probability
        if (runif(1)<exp(new_like+new_prior-old_prior-old_like)) {
          sigma <- new_sigma
        }
      }
    } else{
      u1 <- rnorm(1,0,proposal_sd_wm)
      u2 <- rnorm(1,0,proposal_sd_wm)
      temp_new_sigma_left <- sigma[1]^2+u1
      temp_new_sigma_right <- sigma[2]^2+u2
      if (temp_new_sigma_left<=0|temp_new_sigma_right<=0) {
      } else{
        new_sigma_left <- sqrt(sigma[1]^2+u1)
        new_sigma_right <- sqrt(sigma[2]^2+u2)
        
        old_prior <- dscaleinvchi(sigma[1]^2,sigma_prior[1],sigma_prior[2],log = TRUE)+dscaleinvchi(sigma[2]^2,sigma_prior[1],sigma_prior[2],log = TRUE)
        new_prior <- dscaleinvchi(new_sigma_left^2,sigma_prior[1],sigma_prior[2],log = TRUE)+dscaleinvchi(new_sigma_right^2,sigma_prior[1],sigma_prior[2],log = TRUE)
        
        old_like <- sum(dnorm(y[1:cp],known_mu,sigma[1],log = TRUE))+sum(dnorm(y[(cp+1):(2*cp)],known_mu,sigma[2],log = TRUE))
        new_like <- sum(dnorm(y[1:cp],known_mu,new_sigma_left,log = TRUE))+sum(dnorm(y[(cp+1):(2*cp)],known_mu,new_sigma_right,log = TRUE))
        # Symmetric proposal so no transition probability
        if (runif(1)<exp(new_like+new_prior-old_prior-old_like)) {
          sigma <- c(new_sigma_left,new_sigma_right)
        }
      }
    }
    #now resample mu from its posterior
    # if (num_of_cps==0) {
    #   posteriorsd <- sqrt( ( 1/mu_prior[2]^2 + length(y) )^(-1) )
    #   posteriormu <- sum(y)/ ( 1/mu_prior[2]^2 + length(y))
    #   mu <- rnorm(1,posteriormu,posteriorsd)
    # } else {
    #   y1 <- y[1:cp]; y2 <- y[(cp+1):(2*cp)]
    #   posteriorsd <- sqrt( ( 1/mu_prior[2]^2 + length(y1) )^(-1) )
    #   posteriormu <- sum(y1)/ ( 1/mu_prior[2]^2 + length(y1))
    #   mu[1] <- rnorm(1,posteriormu,posteriorsd)
    #
    #   posteriorsd <- sqrt( ( 1/mu_prior[2]^2 + length(y2) )^(-1) )
    #   posteriormu <- sum(y2)/ ( 1/mu_prior[2]^2 + length(y2))
    #   mu[2] <- rnorm(1,posteriormu,posteriorsd)
    # }
    num_of_cps_list[s] <- num_of_cps
    sigma_list[[s]] <- sigma
  }
  num_of_cps_list <- num_of_cps_list[-(1:(perc_burnin*sims))]
  sigma_list <- sigma_list[-(1:(perc_burnin*sims))]

  cps_0 <- sum(num_of_cps_list==0)
  cps_1 <- sum(num_of_cps_list==1)
  cps_0_perc <- cps_0/(cps_0+cps_1)
  cps_1_perc <- cps_1/(cps_0+cps_1)

  results <- cbind(results,"RJMCMC"=c(cps_0_perc,cps_1_perc))
  
  # Analytical
  marg_like_sigma <- function(y,mu,df_0,scale_0){
    return(((1/sqrt(2*pi))^(length(y)))*(((scale_0*df_0)/2)^(df_0/2))/gamma(df_0/2)*gamma((df_0+length(y))/2)*((2/((df_0*scale_0)+(sum((y-mu)^2))))^((df_0+length(y))/2)))
  }
  
  # For Model M0
  marg_like_sigma_M0 <- marg_like_sigma(y,known_mu,sigma_prior[1],sigma_prior[2])
  
  # For Model M1
  marg_like_sigma_M1 <- marg_like_sigma(y[1:cp],known_mu,sigma_prior[1],sigma_prior[2])*marg_like_sigma(y[(cp+1):(2*cp)],known_mu,sigma_prior[1],sigma_prior[2])
  
  total_marg <- marg_like_sigma_M0+marg_like_sigma_M1
  
  marg_like_M0 <- marg_like_sigma_M0/total_marg
  marg_like_M1 <- marg_like_sigma_M1/total_marg
  
  results <- cbind(results,"Analytical"=c(marg_like_M0,marg_like_M1))

} else{ #Unknown mu and sigma^2
  
  # Use approximation using samples from prior, reversible jump, and analytical

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
  # Use approximation, i.e integrate over all mu and sigma
  # ######Need to add approximation for both mu and sigma
  # Marginal likelihood P(Y|M_0)
  test_sigma_M0 <- sqrt(rscaleinvchi(approx_samples,sigma_prior[1],sigma_prior[2]))
  test_sigma_M0_matrix <- matrix(test_sigma_M0,nrow = approx_samples,ncol = length(test_sigma_M0),byrow = TRUE)
  test_mu_M0 <- rnorm(approx_samples,mu_prior[1],mu_prior[2])
  pM0 <- list()
  for (i in 1:100) {
    pM0[i] <- mean(sapply(test_mu_M0,prod_norm,test_sigma_M0[i],y))
  }
  
  pM0 <- mean(sapply(rnorm(approx_samples,mu_prior[1],mu_prior[2]),prod_norm,known_sigma,y))
  # Marginal Likelihood P(Y|M_1)
  pM1 <- mean(sapply(rnorm(approx_samples,mu_prior[1],mu_prior[2]),prod_norm,known_sigma,y[1:cp]))*mean(sapply(rnorm(approx_samples,mu_prior[1],mu_prior[2]),prod_norm,known_sigma,y[(cp+1):(2*cp)]))
  # Normalise
  total <- pM0+pM1
  pM0 <- pM0/total
  pM1 <- pM1/total
  # Store results in data frame
  results <- data.frame("Approx"=c(pM0,pM1))
  
  # Now calculate using reversible jump
  # proposal_sd_rj <- 0.1 #bigger jump for the reversible jump
  # proposal_sd_wm <- 0.1 #within model sample of parameters
  # sims <- 100000
  # perc_burnin <- 0.1
  # num_of_cps_list <- rep(NA,sims)
  # mu_list <- list() #to assign to this list, use [[]]
  # num_of_cps <- 0
  # mu <- mean(y)
  
  for (s in 1:sims) {
    if (s%%10000==0) {
      print(s)
    }
    # If we have no changepoints, we jump to one changepoint at cp
    if (num_of_cps==0) {
      u1 <- rnorm(1,0,proposal_sd_rj)
      u2 <- rnorm(1,0,proposal_sd_rj)
      
      new_mu_left <- mean(y[1:cp])+u1
      new_mu_right <- mean(y[(cp+1):(2*cp)])+u2
      
      old_prior <- dnorm(mu,mu_prior[1],mu_prior[2],log = TRUE)
      new_prior <- dnorm(new_mu_left,mu_prior[1],mu_prior[2],log = TRUE) + dnorm(new_mu_right,mu_prior[1],mu_prior[2],log = TRUE)
      
      old_like <- sum(dnorm(y,mu,known_sigma,log = TRUE))
      new_like <- sum(dnorm(y[1:cp],new_mu_left,known_sigma,log = TRUE))+sum(dnorm(y[(cp+1):(2*cp)],new_mu_right,known_sigma,log = TRUE))
      
      old_to_new <- dnorm(u1,0,proposal_sd_rj,log = TRUE)+dnorm(u2,0,proposal_sd_rj,log = TRUE)
      new_to_old <- dnorm(mean(y)-mu,0,proposal_sd_rj,log = TRUE)
      if (runif(1)<exp(new_like+new_prior+new_to_old-old_like-old_prior-old_to_new)) {
        num_of_cps <- 1
        mu <- c(new_mu_left,new_mu_right)
      }
    } else { 
      # We are in the one changepoint model
      u <- rnorm(1,0,proposal_sd_rj)
      
      new_mu <- mean(y)+u
      
      old_prior <- dnorm(mu[1],mu_prior[1],mu_prior[2],log = TRUE)+dnorm(mu[2],mu_prior[1],mu_prior[2],log = TRUE)
      new_prior <- dnorm(new_mu,mu_prior[1],mu_prior[2],log = TRUE)
      
      old_like <- sum(dnorm(y[1:cp],mu[1],known_sigma,log = TRUE))+sum(dnorm(y[(cp+1):(2*cp)],mu[2],known_sigma,log = TRUE))
      new_like <- sum(dnorm(y,new_mu,known_sigma,log = TRUE))
      
      old_to_new <- dnorm(u,0,proposal_sd_rj,log = TRUE)
      new_to_old <- dnorm(mu[1]-mean(y[1:cp]),0,proposal_sd_rj,log = TRUE)+dnorm(mu[2]-mean(y[(cp+1):(2*cp)]),0,proposal_sd_rj,log = TRUE)
      
      if (runif(1)< exp(new_like+new_prior-old_like-old_prior+new_to_old-old_to_new)) {
        mu <- new_mu
        num_of_cps <- 0
      }
    }
    # Resample mu parameter from posterior, use Metropolis Hastings
    if (num_of_cps==0) {
      u <- rnorm(1,0,proposal_sd_wm)
      new_mu <- mu+u
      
      old_prior <- dnorm(mu,mu_prior[1],mu_prior[2],log = TRUE)
      new_prior <- dnorm(new_mu,mu_prior[1],mu_prior[2],log = TRUE)
      
      old_like <- sum(dnorm(y,mu,known_sigma,log = TRUE))
      new_like <- sum(dnorm(y,new_mu,known_sigma,log = TRUE))
      # Symmetric proposal so no transition probability
      if (runif(1)<exp(new_like+new_prior-old_prior-old_like)) {
        mu <- new_mu
      }
    } else{
      u1 <- rnorm(1,0,proposal_sd_wm)
      u2 <- rnorm(1,0,proposal_sd_wm)
      new_mu_left <- mu[1]+u1
      new_mu_right <- mu[2]+u2
      
      old_prior <- dnorm(mu[1],mu_prior[1],mu_prior[2],log = TRUE)+dnorm(mu[2],mu_prior[1],mu_prior[2],log = TRUE)
      new_prior <- dnorm(new_mu_left,mu_prior[1],mu_prior[2],log = TRUE)+dnorm(new_mu_right,mu_prior[1],mu_prior[2],log = TRUE)
      
      old_like <- sum(dnorm(y[1:cp],mu[1],known_sigma,log = TRUE))+sum(dnorm(y[(cp+1):(2*cp)],mu[2],known_sigma,log = TRUE))
      new_like <- sum(dnorm(y[1:cp],new_mu_left,known_sigma,log = TRUE))+sum(dnorm(y[(cp+1):(2*cp)],new_mu_right,known_sigma,log = TRUE))
      # Symmetric proposal so no transition probability
      if (runif(1)<exp(new_like+new_prior-old_prior-old_like)) {
        mu <- c(new_mu_left,new_mu_right)
      }
    }
    #now resample mu from its posterior
    # if (num_of_cps==0) {
    #   posteriorsd <- sqrt( ( 1/mu_prior[2]^2 + length(y) )^(-1) )
    #   posteriormu <- sum(y)/ ( 1/mu_prior[2]^2 + length(y)) 
    #   mu <- rnorm(1,posteriormu,posteriorsd)
    # } else {
    #   y1 <- y[1:cp]; y2 <- y[(cp+1):(2*cp)]
    #   posteriorsd <- sqrt( ( 1/mu_prior[2]^2 + length(y1) )^(-1) )
    #   posteriormu <- sum(y1)/ ( 1/mu_prior[2]^2 + length(y1)) 
    #   mu[1] <- rnorm(1,posteriormu,posteriorsd)
    #   
    #   posteriorsd <- sqrt( ( 1/mu_prior[2]^2 + length(y2) )^(-1) )
    #   posteriormu <- sum(y2)/ ( 1/mu_prior[2]^2 + length(y2)) 
    #   mu[2] <- rnorm(1,posteriormu,posteriorsd)
    # }
    num_of_cps_list[s] <- num_of_cps
    mu_list[[s]] <- mu
  }
  num_of_cps_list <- num_of_cps_list[-(1:(perc_burnin*sims))]
  mu_list <- mu_list[-(1:(perc_burnin*sims))]
  
  cps_0 <- sum(num_of_cps_list==0)
  cps_1 <- sum(num_of_cps_list==1)
  cps_0_perc <- cps_0/(cps_0+cps_1)
  cps_1_perc <- cps_1/(cps_0+cps_1)
  
  results <- cbind(results,"RJMCMC"=c(cps_0_perc,cps_1_perc))
  
  # Analytical results
  # We need to calculate marginal likelihood since p(M0|Y) proportional to p(Y|M0)p(M0) and p(M1|Y) proportional to p(Y|M1)p(M1)
  
  # Write function for Marginal likelihood
  marg_like_mu <- function(y,sigma,mu_0,sigma_0){
    return((sigma/(((sqrt(2*pi)*sigma)^length(y))*sqrt(length(y)*sigma_0^2+sigma^2)))*exp(-(sum(y^2)/(2*sigma^2))-(mu_0^2/(2*sigma_0^2)))*exp(((sigma_0^2*length(y)^2*mean(y)^2/sigma^2)+((sigma^2*mu_0^2)/sigma_0^2)+(2*length(y)*mean(y)*mu_0))/(2*(length(y)*sigma_0^2+sigma^2))))              
  }
  
  # For model M0, no changepoint
  # marg_like_M0 <- (known_sigma/((sqrt(2*pi)*known_sigma)*sqrt(length(y)*mu_prior[2]^2+known_sigma^2)))*exp(-(sum(y^2)/(2*known_sigma^2))-(mu_prior[1]^2/(2*mu_prior[2]^2)))*exp(((mu_prior[2]^2*length(y)^2*mean(y)^2/known_sigma^2)+((known_sigma^2*mu_prior[1]^2)/mu_prior[2]^2)+(2*length(y)*mean(y)*mu_prior[1]))/(2*(length(y)*mu_prior[2]^2+known_sigma^2)))              
  marg_like_M0 <- marg_like_mu(y,known_sigma,mu_prior[1],mu_prior[2])
  
  # For model M1, one changepoint
  marg_like_M1 <- marg_like_mu(y[1:cp],known_sigma,mu_prior[1],mu_prior[2])*marg_like_mu(y[(cp+1):(2*cp)],known_sigma,mu_prior[1],mu_prior[2])
  
  total_marg <- marg_like_M0+marg_like_M1
  
  marg_like_M0 <- marg_like_M0/total_marg
  marg_like_M1 <- marg_like_M1/total_marg
  
  results <- cbind(results,"Analytical"=c(marg_like_M0,marg_like_M1))
  

  
  
}
