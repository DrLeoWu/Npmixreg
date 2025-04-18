




library(mixtools)
library(parallel)

#######################set functional parameters
pi_1 <- function(x){
  return(exp(0.5*x)/(1+exp(0.5*x))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(3-sin(2*pi*x))
}

m_2 <- function(x){
  return(cos(3*pi*x))
}

sigma_1 <- function(x){
  return(0.6*exp(0.5*x))
}

sigma_2 <- function(x){
  return(0.5*exp(-0.2*x))
}
####################################################################################################
##### define a function to generate samples from the nonparametric normal mixture regression
#### x is a vector of univariate covariates 
#### p, m, sgm are n by c matrix of mix proportions, means, and standard deviations of normal distributions
rnpnormmix <- function(p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}
####################################################################################################
#### define the Epanechnikov kernel function with bandwith h
EpaKernel <- function(u,h){
  return(3*(1-(u/h)^2)/4*1*(abs(u)<=h)*(1/h))
}
####################################################################################################
##### define a function to implement the modified EM agorithm
##### N is the number of grid points, initl_m, initl_sgm,initl_pi are n by c (c=2)
npnormmixEM <- function(x,y,N=100,initl_m,initl_sgm,initl_pi,maxiter=1000,epsilon=0.00001,bw){
  u <- seq(min(x),max(x),l=N)
  l <- 0
  m_EM <- initl_m
  sgm_EM <- initl_sgm
  pi_EM <- initl_pi
  theta_old <- rep(0,6*N)
  theta_diff <- epsilon 
  
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    r <- pi_EM*dnorm(cbind(y,y),mean=m_EM,sd=sgm_EM)/rowSums(pi_EM*dnorm(cbind(y,y),mean=m_EM,sd=sgm_EM))
    pi1 <- rep(0,N)
    pi2 <- rep(0,N)
    m1 <- rep(0,N)
    m2 <- rep(0,N)
    sgm2_1 <- rep(0,N)
    sgm2_2 <- rep(0,N)
    
    for(j in 1:N){
      pi1[j] <- r[,1]%*%EpaKernel(x-u[j],bw)/sum(EpaKernel(x-u[j],bw))
      pi2[j] <- 1-pi1[j]
      #pi2[j] <- r[,2]%*%EpaKernel(x-u[j],bw)/sum(EpaKernel(x-u[j],bw))
      m1[j] <- (r[,1]*EpaKernel(x-u[j],bw))%*%y/(r[,1]%*%EpaKernel(x-u[j],bw))
      m2[j] <- (r[,2]*EpaKernel(x-u[j],bw))%*%y/(r[,2]%*%EpaKernel(x-u[j],bw))
      sgm2_1[j] <- (r[,1]*EpaKernel(x-u[j],bw))%*%((y-m1[j])^2)/(r[,1]%*%EpaKernel(x-u[j],bw))
      sgm2_2[j] <- (r[,2]*EpaKernel(x-u[j],bw))%*%((y-m2[j])^2)/(r[,2]%*%EpaKernel(x-u[j],bw))
    }
    
    if(min(c(sgm2_1,sgm2_2))==0){break}
    
    theta_new <- c(pi1,pi2,m1,m2,sgm2_1,sgm2_2)
    theta_diff <- max(abs(theta_new - theta_old))
    theta_old <- theta_new
    
    pi_EM[,1] <- approx(u,pi1,xout=x)$y
    pi_EM[,2] <- 1-pi_EM[,1]
    #pi_EM[,2] <- approx(u,pi2,xout=x)$y
    m_EM[,1] <- approx(u,m1,xout=x)$y
    m_EM[,2] <- approx(u,m2,xout=x)$y
    sgm_EM[,1] <- sqrt(approx(u,sgm2_1,xout=x)$y)
    sgm_EM[,2] <- sqrt(approx(u,sgm2_2,xout=x)$y)
  }
  if (l >= maxiter) {print("The algorithm does not converge")}
  if (min(c(sgm2_1,sgm2_2))==0) {
    print("Need larger bandwidth")
    return(list(pi1=rep(NA,length(y)),pi2=rep(NA,length(y)),m1=rep(NA,length(y)),m2=rep(NA,length(y)),sigma2_1=rep(NA,length(y)),sigma2_2=rep(NA,length(y)),iter=l))
  }else{return(list(pi1=pi_EM[,1],pi2=pi_EM[,2],m1=m_EM[,1],m2=m_EM[,2],sigma2_1=sgm_EM[,1]^2,sigma2_2=sgm_EM[,2]^2,iter=l))}
}
###########################################  n =200   ###############################################################


set.seed(293)
n <- 200
x <- runif(n)
p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))

kdde_0 <- kdde(x,deriv.order=0)
fx <- predict(kdde_0,x=x)
kdde_1 <- kdde(x,deriv.order=1)
grad_fx <- predict(kdde_1,x=x)

bs_AMISE <- function(j){
  library(mixtools)
  library(ks)
  library(fpc)
  library(sn)
  
  y <- rnpnormmix(p,mu,sgm)
  rmt <- regmix(cbind(x,x^2,x^3),y,nclust=2)
  
  coef_pi <- c(log(rmt$eps[1]/(1-rmt$eps[1])),0)
  coef_m1 <- rmt$coef[,1]
  coef_m2 <- rmt$coef[,2]
  coef_sigma21 <- c(log(rmt$var[1]),0)
  coef_sigma22 <- c(log(rmt$var[2]),0)
  
  X_pi <- X_sig2 <- cbind(rep(1,n),x) ## linear approximation for proportion functions and variance functions 
  X_m <- cbind(rep(1,n),x,x^2,x^3)
  parGMR <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi=coef_pi,initl_beta_m1=coef_m1,initl_beta_m2=coef_m2,initl_beta_sgm21=coef_sigma21,initl_beta_sgm22=coef_sigma22)
  coef_pi <- -parGMR$beta_pi
  coef_sigma22 <- parGMR$beta_sgm22
  coef_sigma21 <- parGMR$beta_sgm21
  coef_m1 <- parGMR$beta_m1
  coef_m2 <- parGMR$beta_m2
  
  est_unknown <- function(x,coef_pi,coef_m1,coef_m2,coef_sigma21,coef_sigma22){
    
    xd_pi <- xd_sigma <- c(1,x)
    xd_m <- c(1,x,x^2,x^3)
    p <- exp(c(coef_pi%*%xd_pi))/(1+exp(c(coef_pi%*%xd_pi)))
    m1 <- c(coef_m1%*%xd_m)
    m2 <- c(coef_m2%*%xd_m)
    sigma21 <- exp(c(coef_sigma21%*%xd_sigma))
    sigma22 <- exp(c(coef_sigma22%*%xd_sigma))
    
    deriv_theta_1 <- p^2*exp(-c(coef_pi%*%xd_pi))*(-coef_pi[2])
    deriv_theta_2 <- sigma21*coef_sigma21[2]
    deriv_theta_3 <- sigma22*coef_sigma22[2]
    deriv_theta_4 <- c(coef_m1[-1]%*%c(1,2*x,3*x^2))
    deriv_theta_5 <- c(coef_m2[-1]%*%c(1,2*x,3*x^2))
    deriv_theta <- c(deriv_theta_1,deriv_theta_2,deriv_theta_3,deriv_theta_4,deriv_theta_5)
    
    deriv2_theta_1 <- -exp(c(coef_pi%*%xd_pi))*(1-exp(c(coef_pi%*%xd_pi)))/(1+exp(c(coef_pi%*%xd_pi)))^3*coef_pi[2]^2
    deriv2_theta_2 <- sigma21*coef_sigma21[2]^2
    deriv2_theta_3 <- sigma22*coef_sigma22[2]^2
    deriv2_theta_4 <- c(coef_m1[-c(1,2)]%*%c(2,6*x))
    deriv2_theta_5 <- c(coef_m2[-c(1,2)]%*%c(2,6*x))
    deriv2_theta <- c(deriv2_theta_1,deriv2_theta_2,deriv2_theta_3,deriv2_theta_4,deriv2_theta_5)
    
    y_temp <- rnorm(10000)
    FI_integ <- function(y){
      eta <- p*dnorm(y,mean=m1,sd=sqrt(sigma21))+(1-p)*dnorm(y,mean=m2,sd=sqrt(sigma22))
      deriv_eta_1 <- dnorm(y,mean=m1,sd=sqrt(sigma21)) - dnorm(y,mean=m2,sd=sqrt(sigma22))
      deriv_eta_2 <- p*dnorm(y,mean=m1,sd=sqrt(sigma21))*(-1/(2*sigma21)+(y-m1)^2/(2*sigma21^2))
      deriv_eta_3 <- (1-p)*dnorm(y,mean=m2,sd=sqrt(sigma22))*(-1/(2*sigma22)+(y-m2)^2/(2*sigma22^2))
      deriv_eta_4 <- p*(y-m1)/sigma21*dnorm(y,mean=m1,sd=sqrt(sigma21))
      deriv_eta_5 <- (1-p)*(y-m2)/sigma22*dnorm(y,mean=m2,sd=sqrt(sigma22))
      deriv_eta <- c(deriv_eta_1,deriv_eta_2,deriv_eta_3,deriv_eta_4,deriv_eta_5)
      deriv2_eta_11 <- 0 
      deriv2_eta_12 <- dnorm(y,mean=m1,sd=sqrt(sigma21))*(-1/(2*sigma21)+(y-m1)^2/(2*sigma21^2))
      deriv2_eta_13 <- -dnorm(y,mean=m2,sd=sqrt(sigma22))*(-1/(2*sigma22)+(y-m2)^2/(2*sigma22^2))
      deriv2_eta_14 <- dnorm(y,mean=m1,sd=sqrt(sigma21))*(y-m1)/sigma21
      deriv2_eta_15 <- -dnorm(y,mean=m2,sd=sqrt(sigma22))*(y-m2)/sigma22
      deriv2_eta_21 <- deriv2_eta_12
      deriv2_eta_22 <- p*dnorm(y,mean=m1,sd=sqrt(sigma21))*(1/(2*sigma21^2)-(y-m1)^2/(sigma21^3)) + p*dnorm(y,mean=m1,sd=sqrt(sigma21))*(-1/(2*sigma21)+(y-m1)^2/(2*sigma21^2))^2
      deriv2_eta_23 <- 0
      deriv2_eta_24 <- p*(-y+m1)/sigma21^2*dnorm(y,mean=m1,sd=sqrt(sigma21))+p*(-1/(2*sigma21)+(y-m1)^2/(2*sigma21^2))*(y-m1)/sigma21*dnorm(y,mean=m1,sd=sqrt(sigma21))
      deriv2_eta_25 <- 0
      deriv2_eta_31 <- deriv2_eta_13
      deriv2_eta_32 <- 0
      deriv2_eta_33 <- (1-p)*dnorm(y,mean=m2,sd=sqrt(sigma22))*(1/(2*sigma22^2)-(y-m2)^2/(sigma22^3))+ (1-p)*dnorm(y,mean=m2,sd=sqrt(sigma22))*(-1/(2*sigma22)+(y-m2)^2/(2*sigma22^2))
      deriv2_eta_34 <- 0
      deriv2_eta_35 <- (1-p)*(-y+m2)/sigma22^2*dnorm(y,mean=m2,sd=sqrt(sigma22))+(1-p)*dnorm(y,mean=m2,sd=sqrt(sigma22))*(-1/(2*sigma22)+(y-m2)^2/(2*sigma22^2))*(y-m2)/sigma22
      deriv2_eta_41 <- deriv2_eta_14 
      deriv2_eta_42 <-  deriv2_eta_24
      deriv2_eta_43 <- 0
      deriv2_eta_44 <- -p/sigma21*dnorm(y,mean=m1,sd=sqrt(sigma21)) + p*(y-m1)^2/(sigma21)^2*dnorm(y,mean=m1,sd=sqrt(sigma21))
      deriv2_eta_45 <- 0
      deriv2_eta_51 <- deriv2_eta_15
      deriv2_eta_52 <- 0
      deriv2_eta_53 <- deriv2_eta_35
      deriv2_eta_54 <- 0
      deriv2_eta_55 <- -(1-p)/sigma22*dnorm(y,mean=m2,sd=sqrt(sigma22)) + (1-p)*(y-m2)^2/(sigma22)^2*dnorm(y,mean=m2,sd=sqrt(sigma22))
      deriv2_eta <- matrix(c(deriv2_eta_11,deriv2_eta_21,deriv2_eta_31,deriv2_eta_41,deriv2_eta_51,
                             deriv2_eta_12,deriv2_eta_22,deriv2_eta_32,deriv2_eta_42,deriv2_eta_52,
                             deriv2_eta_13,deriv2_eta_23,deriv2_eta_33,deriv2_eta_43,deriv2_eta_53,
                             deriv2_eta_14,deriv2_eta_24,deriv2_eta_34,deriv2_eta_44,deriv2_eta_54,
                             deriv2_eta_15,deriv2_eta_25,deriv2_eta_35,deriv2_eta_45,deriv2_eta_55),nrow=5)
      #q2 <- 1/eta*deriv2_eta - 1/eta^2*deriv_eta%*%t(deriv_eta)
      q2 <- (deriv2_eta - 1/eta*deriv_eta%*%t(deriv_eta))/dnorm(y)
      return(q2)
    }
    
    Lambda_integ <- function(y){
      eta <- p*dnorm(y,mean=m1,sd=sqrt(sigma21))+(1-p)*dnorm(y,mean=m2,sd=sqrt(sigma22))
      deriv_eta_1 <- dnorm(y,mean=m1,sd=sqrt(sigma21)) - dnorm(y,mean=m2,sd=sqrt(sigma22))
      deriv_eta_2 <- p*dnorm(y,mean=m1,sd=sqrt(sigma21))*(-1/(2*sigma21)+(y-m1)^2/(2*sigma21^2))
      deriv_eta_3 <- (1-p)*dnorm(y,mean=m2,sd=sqrt(sigma22))*(-1/(2*sigma22)+(y-m2)^2/(2*sigma22^2))
      deriv_eta_4 <- p*(y-m1)/sigma21*dnorm(y,mean=m1,sd=sqrt(sigma21))
      deriv_eta_5 <- (1-p)*(y-m2)/sigma22*dnorm(y,mean=m2,sd=sqrt(sigma22))
      deriv_eta <- c(deriv_eta_1,deriv_eta_2,deriv_eta_3,deriv_eta_4,deriv_eta_5)
      return(1/eta*deriv_eta%*%t(deriv_eta)/dnorm(y))
    }
    
    FI <- -matrix(apply(sapply(y_temp, FI_integ),1,mean),nrow=5)
    Lambda_core <- matrix(apply(sapply(y_temp, Lambda_integ),1,mean),nrow=5)
    deriv_lambda <- c(Lambda_core%*%deriv_theta)
    deriv2_lambda <- c(Lambda_core%*%deriv2_theta)
    return(list(FI=FI,deriv_lambda=deriv_lambda,deriv2_lambda=deriv2_lambda))
    
  }
  A <- B <- 0
  for(i in 1:n){
    quant_unknow <- est_unknown(x[i],coef_pi,coef_m1,coef_m2,coef_sigma21,coef_sigma22)
    A <- A + 3/5*tr(solve(quant_unknow$FI))/fx[i]/n
    B <- B + c((grad_fx[i]/fx[i]*quant_unknow$deriv_lambda+1/2*quant_unknow$deriv2_lambda)%*%solve(quant_unknow$FI)%*%solve(quant_unknow$FI)%*%(grad_fx[i]/fx[i]*quant_unknow$deriv_lambda+1/2*quant_unknow$deriv2_lambda)/25)
    print(i)
  }
  
  h_AMISE <- (A/(4*B))^{1/5}
  return(h_AMISE)
}

## generate 10 simulations 
the_cluster_bs_AMISE <- makeCluster(10)
clusterSetRNGStream(the_cluster_bs_AMISE,27)
clusterExport(the_cluster_bs_AMISE,c("rnpnormmix","emopt_pi","emopt_wnorm","normmixEM","n","x","p","mu","sgm","grad_fx","fx"))
h_AMISE <- clusterCall(cl = the_cluster_bs_AMISE, bs_AMISE,1:10)
stopCluster(the_cluster_bs_AMISE)

h_opt <- mean(unlist(h_AMISE),na.rm = TRUE) #0.08948493 0.07174553 0.07234931
h_opt <- 0.08948493
h_opt <- 0.11
##############################################################################################################################
############## Conduct simulation studies (500 samples, n=200)
set.seed(579)
initl_m <- cbind(m_1(x)+rnorm(1,sd=0.5),m_2(x)+rnorm(1,sd=0.5))  ### set initial values for the mean function 
initl_sgm <- cbind(rep(mean(sigma_1(x)),n),rep(mean(sigma_2(x)),n))   ## set the initial values for the standard deviation function
initl_pi <- cbind(rep(mean(pi_1(x)),n),rep(mean(pi_2(x)),n))  ## set the initial values for the mixing proportion function

## define a function for parallel simulation
sim_npmixreg <- function(j){
  # load libraries
  library(mixtools)
  library(missForest) ## library for imputation
  library(interp)
  y <- rnpnormmix(p,mu,sgm)
  output <- npnormmixEM(x,y,N=100,initl_m=initl_m,initl_sgm=initl_sgm,initl_pi=initl_pi,bw=h_opt) 
  RASE_pi <- sqrt((sum((output$pi1-pi_1(x))^2)+sum((output$pi2-pi_2(x))^2))/n)
  RASE_m  <- sqrt((sum((output$m1-m_1(x))^2)+sum((output$m2-m_2(x))^2))/n)
  RASE_sigma2 <- sqrt((sum((output$sigma2_1-sigma_1(x)^2)^2)+sum((output$sigma2_2-sigma_2(x)^2)^2))/n)
  return(list(RASE_pi=RASE_pi,RASE_m=RASE_m,RASE_sigma2=RASE_sigma2))
}

n_sim <- 500
the_cluster <- makeCluster(16)
clusterSetRNGStream(the_cluster,27)
clusterExport(the_cluster,c("npnormmixEM","EpaKernel","rnpnormmix","x","n","h_opt","initl_m","initl_sgm","initl_pi","pi_1",
                            "pi_2","m_1","m_2","sigma_1","sigma_2","p","mu","sgm"))
test_sim <- parLapply(the_cluster,1:n_sim,sim_npmixreg)
stopCluster(the_cluster)

RASEs <- matrix(unlist(test_sim),nrow=3,byrow = FALSE)
RASE_mean <- apply(RASEs,1,mean)
RASE_sd <- apply(RASEs,1,sd)




