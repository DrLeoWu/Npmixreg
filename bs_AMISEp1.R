
library(mixtools)
library(fpc)
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
##############################################################################################################################
##############  Define relavent functions for bandwidth selection based on AMISE
emopt_pi <- function(initl,X,w){ # initl=initial value for beta, X is n by p, w is n by 1
  
  emobj_pi <- function(beta_pi){
    nsamp <- nrow(X)
    return(-c(t(w)%*%X%*%beta_pi)/nsamp + sum(log(1+exp(X%*%beta_pi)))/nsamp)
  }
  
  emobj_pi_gr <- function(beta_pi){
    nsamp <- nrow(X)
    return(-t(X)%*%w/nsamp + t(X)%*%(1/(1+exp(-X%*%beta_pi)))/nsamp)
  }
  
  res <- optim(initl,emobj_pi,emobj_pi_gr,method = "BFGS")
  return(res$par)
}

emopt_wnorm <- function(initl_beta_m,initl_beta_sig2,X_m,X_sig2,y,w,maxiter=10000,epsilon=0.0001){  
  # this function solve the optimization problem of weighted log normal likelihood
  nsamp <- nrow(X_m)
  beta1_old <- initl_beta_m
  beta2_old <- initl_beta_sig2
  #w_new <- diag(w*exp(-X_sig2%*%beta2))
  l <- 0
  beta_diff <- epsilon 
  while(l < maxiter & beta_diff >= epsilon){
    l <- l + 1
    w_new <- diag(c(w*exp(-X_sig2%*%beta2_old)))
    beta1_new <- c(solve(t(X_m)%*%w_new%*%X_m)%*%t(X_m)%*%w_new%*%y)
    
    emobj_beta_sig2 <- function(beta_sig2){
      return(sum(w*(y-X_m%*%beta1_new)^2*exp(-X_sig2%*%beta_sig2)) + c(t(w)%*%(X_sig2%*%beta_sig2)))
    }
    
    emobj_beta_sig2_gr <- function(beta_sig2){
      return(-t(X_sig2)%*%(w*(y-X_m%*%beta1_new)^2*exp(-X_sig2%*%beta_sig2)) + t(X_sig2)%*%w)
    }
    
    rep2 <- optim(initl_beta_sig2,emobj_beta_sig2,emobj_beta_sig2_gr,method = "BFGS")
    beta2_new <- rep2$par
    beta_diff <- max(abs(c(beta1_new-beta1_old,beta2_new-beta2_old)))
    beta1_old <- beta1_new
    beta2_old <- beta2_new
  }
  
  return(list(beta_m=beta1_old,beta_sig2=beta2_old))
}

normmixEM <- function(x_pi,x_m,x_sgm2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22,maxiter=10000,epsilon=0.0001){ #x is n by p, and initl_beta is p by 5 (only for two clusters)
  beta_pi_old <- initl_beta_pi 
  beta_m1_old <- initl_beta_m1 
  beta_m2_old <- initl_beta_m2 
  beta_sgm21_old <- initl_beta_sgm21 
  beta_sgm22_old <- initl_beta_sgm22 
  w_base <- rep(0,length(y))
  l <- 0
  w_diff <- epsilon 
  
  while(l < maxiter & w_diff >= epsilon){
    l <- l + 1
    em_pi <- 1/(1+exp(-c(x_pi%*%beta_pi_old)))
    em_m1 <- c(x_m%*%beta_m1_old)
    em_m2 <- c(x_m%*%beta_m2_old)
    em_sgm21 <- exp(c(x_sgm2%*%beta_sgm21_old))
    em_sgm22 <- exp(c(x_sgm2%*%beta_sgm22_old))
    w_update <- (em_pi*dnorm(y,mean=em_m1,sd=sqrt(em_sgm21)))/(em_pi*dnorm(y,mean=em_m1,sd=sqrt(em_sgm21))+(1-em_pi)*dnorm(y,mean=em_m2,sd=sqrt(em_sgm22)))
    beta_pi_new <- emopt_pi(beta_pi_old,x_pi,w_update)
    normal1 <- emopt_wnorm(beta_m1_old,beta_sgm21_old,x_m,x_sgm2,y,w_update)
    normal2 <- emopt_wnorm(beta_m2_old,beta_sgm22_old,x_m,x_sgm2,y,1-w_update)
    beta_m1_new <- normal1$beta_m
    beta_m2_new <- normal2$beta_m
    beta_sgm21_new <- normal1$beta_sig2
    beta_sgm22_new <- normal2$beta_sig2
    
    beta_pi_old <- beta_pi_new
    beta_m1_old <- beta_m1_new
    beta_m2_old <- beta_m2_new
    beta_sgm21_old <- beta_sgm21_new
    beta_sgm22_old <- beta_sgm22_new
    w_diff <- max(abs(w_update-w_base))
    w_base <- w_update
    
    #print(l)
    #print(w_base)
    
  }
  if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(beta_pi=beta_pi_old,beta_m1=beta_m1_old,beta_m2=beta_m2_old,beta_sgm21=beta_sgm21_old,beta_sgm22=beta_sgm22_old,iter=l))
}

##############################################################################################################################

set.seed(893)
n <- 400
x <- runif(n)
p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))
y <- rnpnormmix(p,mu,sgm)
rmt <- regmix(cbind(x,x^2,x^3,x^4),y,nclust=2)

kdde_0 <- kdde(x,deriv.order=0)
fx <- predict(kdde_0,x=x)
kdde_1 <- kdde(x,deriv.order=1)
grad_fx <- predict(kdde_1,x=x)

coef_pi <- c(log(rmt$eps[1]/(1-rmt$eps[1])),0)
coef_m1 <- rmt$coef[,1]
coef_m2 <- rmt$coef[,2]
coef_sigma21 <- c(log(rmt$var[1]),0)
coef_sigma22 <- c(log(rmt$var[2]),0)

X_pi <- X_sig2 <- cbind(rep(1,n),x) ## linear approximation for proportion functions and variance functions 
X_m <- cbind(rep(1,n),x,x^2,x^3,x^4)
parGMR <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi=coef_pi,initl_beta_m1=coef_m1,initl_beta_m2=coef_m2,initl_beta_sgm21=coef_sigma21,initl_beta_sgm22=coef_sigma22)
coef_pi <- parGMR$beta_pi
coef_sigma22 <- parGMR$beta_sgm22
coef_sigma21 <- parGMR$beta_sgm21
coef_m1 <- parGMR$beta_m1
coef_m2 <- parGMR$beta_m2

est_unknown <- function(x,coef_pi,coef_m1,coef_m2,coef_sigma21,coef_sigma22){
  
  xd_pi <- xd_sigma <- c(1,x)
  xd_m <- c(1,x,x^2,x^3,x^4)
  p <- exp(c(coef_pi%*%xd_pi))/(1+exp(c(coef_pi%*%xd_pi)))
  m1 <- c(coef_m1%*%xd_m)
  m2 <- c(coef_m2%*%xd_m)
  sigma21 <- exp(c(coef_sigma21%*%xd_sigma))
  sigma22 <- exp(c(coef_sigma22%*%xd_sigma))
  
  deriv_theta_1 <- p^2*exp(-c(coef_pi%*%xd_pi))*(-coef_pi[2])
  deriv_theta_2 <- sigma21*coef_sigma21[2]
  deriv_theta_3 <- sigma22*coef_sigma22[2]
  deriv_theta_4 <- c(coef_m1[-1]%*%c(1,2*x,3*x^2,4*x^3))
  deriv_theta_5 <- c(coef_m2[-1]%*%c(1,2*x,3*x^2,4*x^3))
  deriv_theta <- c(deriv_theta_1,deriv_theta_2,deriv_theta_3,deriv_theta_4,deriv_theta_5)
  
  deriv2_theta_1 <- -exp(c(coef_pi%*%xd_pi))*(1-exp(c(coef_pi%*%xd_pi)))/(1+exp(c(coef_pi%*%xd_pi)))^3*coef_pi[2]^2
  deriv2_theta_2 <- sigma21*coef_sigma21[2]^2
  deriv2_theta_3 <- sigma22*coef_sigma22[2]^2
  deriv2_theta_4 <- c(coef_m1[-c(1,2)]%*%c(2,6*x,12*x^2))
  deriv2_theta_5 <- c(coef_m2[-c(1,2)]%*%c(2,6*x,12*x^2))
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
print(h_AMISE) ## n=200 h_AMISE = 0.08742973

## generate 10 simulations 
the_cluster_bs_AMISE <- makeCluster(10)
clusterSetRNGStream(the_cluster_bs_AMISE,27)
clusterExport(the_cluster_bs_AMISE,c("rnpnormmix","emopt_pi","emopt_wnorm","normmixEM","grad_wosort","deriv_2_ll","deriv_1_lambda"
                                     ,"deriv_2_lambda","mse_pt_normmix","n","x","X_pi","X_sig2","X_m","p","mu","sgm"))
h_AMISE <- clusterCall(cl = the_cluster_bs_AMISE, bs_AMISE,1:10)
stopCluster(the_cluster_bs_AMISE)

## find the average of the 10 selected bandwidth
h_sum <- matrix(rep(0,4),nrow=2)
for(i in 1:10){
  h_sum <- h_sum + h_AMISE[[i]]
}
h_AMISE <- h_sum/10

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

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

mean(unlist(h_AMISE),na.rm = TRUE) #0.08948493 0.07174553 0.07234931







