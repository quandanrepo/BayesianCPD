# In this Gibbs sampler, they use a cumulative sum on the conditional distribution of the changepoint location and compare it to a random uniform number

MC<-1000#1000 # number of draws (chains)
N<-200 # length of chains
Y<-read.csv("data/coal.csv") # read in data
Y <- Y[,2]
n<-length(Y) # number of observations
m<-n # no change point
p<-rep(0,3*MC*N)#rep(0,3*MC*N) # array to store chains
dim(p)<-c(3,MC,N)#c(3,MC,N)
for (j in (1:MC)) {
  a1<-3 # parameter of priors
  a2<-1
  b1<-0.5
  b2<-0.5
  m<-as.integer(n*runif(1))+1
  for (i in (1:N)) {
    l1<-rgamma(1,a1+sum(Y[1:m]),m+b1) # step 1
    l2<-rgamma(1,a2+sum(Y)-sum(Y[1:m]),n-m+b2)
    log_pm<-(l2-l1)*(1:n)+cumsum(Y)*log(l1/l2)#exp((l2-l1)*(1:n))*(l1/l2)^cumsum(Y) # step 2 #take log as we get into infinity errors
    log_pm<-log_pm-max(log_pm) #as the numbers can get extremely small
    pm <- exp(log_pm) #back to normal space
    pm <- pm/sum(pm)
    m<-min(which(runif(1)<cumsum(pm)))#min((1:n)[runif(1)<cumsum(pm)])
    p[1,j,i]<-m # save result
    p[2,j,i]<-l1
    p[3,j,i]<-l2
  }
}

m_freq <- as.data.frame(table(factor(as.vector(p[1,,]))))
m_den <- m_freq[,2]/sum(m_freq[,2])
