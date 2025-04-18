##########################################################################################################
##########################################################################################################
####### Simulation for nonparametric mixture of regression (Scenario 2, n=400, EX2 independent covariates)
##########################################################################################################
##########################################################################################################

## load libraries
library(mixtools)
library(sn)
library(matrixcalc)
library(Matrix)
library(ks)
library(pracma)
library(interp)
library(parallel)
library(dvmisc)
library(fpc)
library("scatterplot3d")
library(rgl)


##############################################################################################################################
##############  Define functions for data generation for a two component mixture of regression models

rnpnormmix <- function(p,m,sgm){ # p, m, sgm are n by c matrix of mix proportions, means, and standard deviations of normal distributions
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}
## define a function for setting correlations 
autocorr.mat <- function(p, rho = 0.5) {
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}

## define functions for simulation setup 
pi_1 <- function(x){
  n <- nrow(x)
  return(rep(0.6,n)) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(4+2*x[,1]+2*x[,2])
}

m_2 <- function(x){
  return(1-2*x[,1]+0.5*x[,2])
}

sigma_1 <- function(x){
  return(rep(0.8,nrow(x)))
}

sigma_2 <- function(x){
  return(rep(0.6,nrow(x)))
}

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

grad_wosort <- function(FF,x){
  temp_id <- 1:length(x)
  temp_id <- temp_id[order(x)]
  x_sort <- x[order(x)]
  F_sort <- FF[order(x)]
  temp_f <- gradient(F_sort, x_sort)
  return(temp_f[order(temp_id)])
}

deriv_2_ll <- function(y,x,thetas){
  ests <- thetas[thetas[,1]==x[1],]
  p1 <- ests[3]
  m1 <- ests[6]
  m2 <- ests[7]
  sig2_1 <- ests[4]
  sig2_2 <- ests[5]
  p2 <- 1-p1
  phi_1 <- dnorm(y,m1,sqrt(sig2_1))
  phi_2 <- dnorm(y,m2,sqrt(sig2_2))
  lh <- p1*phi_1+ p2*phi_2
  ss1 <- -1/(2*sig2_1)+(y-m1)^2/(2*sig2_1^2)
  ss2 <- -1/(2*sig2_2)+(y-m2)^2/(2*sig2_2^2)
  mm1 <- (y-m1)/sig2_1
  mm2 <- (y-m2)/sig2_2
  pf_1 <- c(phi_1-phi_2,p1*ss1*phi_1,p2*ss2*phi_2,p1*mm1*phi_1,p2*mm2*phi_2) ## first deriv
  pf_2 <- matrix(c(0,ss1*phi_1,-ss2*phi_2,mm1*phi_1,-mm2*phi_2,
                   ss1*phi_1,p1*(1/(2*sig2_1^2)-1/(2*sig2_1)-mm1^2/sig2_1+mm1^2/2)*phi_1,0,p1*(mm1^3/2-1.5*mm1/sig2_1)*phi_1,0,
                   -ss2*phi_2,0,p2*(1/(2*sig2_2^2)-1/(2*sig2_2)-mm2^2/sig2_2+mm2^2/2)*phi_2,0,p2*(mm2^3/2-1.5*mm2/sig2_2)*phi_2,
                   mm1*phi_1,p1*(mm1^3/2-1.5*mm1/sig2_1)*phi_1,0,-p1*phi_1/sig2_1+p1*phi_1*mm1^2,0,
                   -mm2*phi_2,0,p2*(mm2^3/2-1.5*mm2/sig2_2)*phi_2,0,-p2*phi_2/sig2_2+p2*phi_2*mm2^2),5,5,byrow=T)
  int_FIM <- 1/lh*pf_2-pf_1%*%t(pf_1)/(lh)^2
  
  return(int_FIM)
}

deriv_1_lambda <- function(y,x,thetas,ptheta_1,ptheta_2){
  ests <- thetas[thetas[,1]==x[1],]
  p1 <- ests[3]
  m1 <- ests[6]
  m2 <- ests[7]
  sig2_1 <- ests[4]
  sig2_2 <- ests[5]
  p2 <- 1-p1
  phi_1 <- dnorm(y,m1,sqrt(sig2_1))
  phi_2 <- dnorm(y,m2,sqrt(sig2_2))
  lh <- p1*phi_1+ p2*phi_2
  ss1 <- -1/(2*sig2_1)+(y-m1)^2/(2*sig2_1^2)
  ss2 <- -1/(2*sig2_2)+(y-m2)^2/(2*sig2_2^2)
  mm1 <- (y-m1)/sig2_1
  mm2 <- (y-m2)/sig2_2
  pf_1 <- c(phi_1-phi_2,p1*ss1*phi_1,p2*ss2*phi_2,p1*mm1*phi_1,p2*mm2*phi_2)
  eptheta_1 <- ptheta_1[ptheta_1[,1]==x[1],-1]
  eptheta_2 <- ptheta_2[ptheta_2[,1]==x[2],-1]
  E <- c(c(t(pf_1)%*%eptheta_1),c(t(pf_1)%*%eptheta_2))
  int_plambda_1 <- 1/lh^2*kronecker(t(pf_1),E)
  return(int_plambda_1)
}

deriv_2_lambda <- function(y,x,thetas,ptheta_1,ptheta_2,ptheta_11,ptheta_22,ptheta_12){
  ests <- thetas[thetas[,1]==x[1],]
  p1 <- ests[3]
  m1 <- ests[6]
  m2 <- ests[7]
  sig2_1 <- ests[4]
  sig2_2 <- ests[5]
  p2 <- 1-p1
  phi_1 <- dnorm(y,m1,sqrt(sig2_1))
  phi_2 <- dnorm(y,m2,sqrt(sig2_2))
  lh <- p1*phi_1+ p2*phi_2
  ss1 <- -1/(2*sig2_1)+(y-m1)^2/(2*sig2_1^2)
  ss2 <- -1/(2*sig2_2)+(y-m2)^2/(2*sig2_2^2)
  mm1 <- (y-m1)/sig2_1
  mm2 <- (y-m2)/sig2_2
  pf_1 <- c(phi_1-phi_2,p1*ss1*phi_1,p2*ss2*phi_2,p1*mm1*phi_1,p2*mm2*phi_2) ## first deriv
  pf_2 <- matrix(c(0,ss1*phi_1,-ss2*phi_2,mm1*phi_1,-mm2*phi_2,
                   ss1*phi_1,p1*(1/(2*sig2_1^2)-1/(2*sig2_1)-mm1^2/sig2_1+mm1^2/2)*phi_1,0,p1*(mm1^3/2-1.5*mm1/sig2_1)*phi_1,0,
                   -ss2*phi_2,0,p2*(1/(2*sig2_2^2)-1/(2*sig2_2)-mm2^2/sig2_2+mm2^2/2)*phi_2,0,p2*(mm2^3/2-1.5*mm2/sig2_2)*phi_2,
                   mm1*phi_1,p1*(mm1^3/2-1.5*mm1/sig2_1)*phi_1,0,-p1*phi_1/sig2_1+p1*phi_1*mm1^2,0,
                   -mm2*phi_2,0,p2*(mm2^3/2-1.5*mm2/sig2_2)*phi_2,0,-p2*phi_2/sig2_2+p2*phi_2*mm2^2),5,5,byrow=T)
  eptheta_1 <- ptheta_1[ptheta_1[,1]==x[1],-1]
  eptheta_11 <- ptheta_11[ptheta_11[,1]==x[1],-1]
  eptheta_2 <- ptheta_2[ptheta_2[,1]==x[2],-1]
  eptheta_22 <- ptheta_22[ptheta_22[,1]==x[2],-1]
  eptheta_12 <- ptheta_12[ptheta_12[,1]==x[1],-1]
  E1 <- c(t(eptheta_1)%*%pf_2%*%eptheta_1+c(t(pf_1)%*%eptheta_11))
  E2 <- c(t(eptheta_1)%*%pf_2%*%eptheta_2+c(t(pf_1)%*%eptheta_12))
  E3 <- c(t(eptheta_2)%*%pf_2%*%eptheta_2+c(t(pf_1)%*%eptheta_22))
  E <- matrix(c(E1,E2,E2,E3),2)
  
  int_plambda_2 <- t(1/lh^2*kronecker(pf_1,E))
  return(int_plambda_2)
}

mse_pt_normmix <- function(x0,Nsim=10000,thetas,ptheta_1,ptheta_2,ptheta_11,ptheta_22,ptheta_12,fx,grad_fx){
  set.seed(100)
  ests <- thetas[thetas[,1]==x0[1],]
  fxhat <- c(fx[fx[,1]==x0[1],-1])
  grad_fxhat <- grad_fx[grad_fx[,1]==x0[1],-1]
  p1 <- ests[3]
  m1 <- ests[6]
  m2 <- ests[7]
  sig2_1 <- ests[4]
  sig2_2 <- ests[5]
  y <- rnormmix(Nsim, c(p1,1-p1), c(m1,m2), c(sqrt(sig2_1),sqrt(sig2_2)))
  FIs<- sapply(y,deriv_2_ll,x=x0,thetas=thetas)
  FIM <- -matrix(apply(FIs,1,mean),5,5)
  plambda1s <- sapply(y,deriv_1_lambda,x=x0,thetas=thetas,ptheta_1=ptheta_1,ptheta_2=ptheta_2)
  plambda_1 <- matrix(apply(plambda1s,1,mean),nrow=2)
  plambda2s <- sapply(y,deriv_2_lambda,x=x0,thetas=thetas,ptheta_1=ptheta_1,ptheta_2=ptheta_2,ptheta_11=ptheta_11,ptheta_22=ptheta_22,ptheta_12=ptheta_12)
  plambda_2 <- matrix(apply(plambda2s,1,mean),nrow=4) 
  D_2 <- duplication.matrix(2)
  A1 <- 2/fxhat*vec(grad_fxhat%*%t(plambda_1[,1])) + plambda_2[,1]
  A2 <- 2/fxhat*vec(grad_fxhat%*%t(plambda_1[,2])) + plambda_2[,2]
  A3 <- 2/fxhat*vec(grad_fxhat%*%t(plambda_1[,3])) + plambda_2[,3]
  A4 <- 2/fxhat*vec(grad_fxhat%*%t(plambda_1[,4])) + plambda_2[,4]
  A5 <- 2/fxhat*vec(grad_fxhat%*%t(plambda_1[,5])) + plambda_2[,5]
  B <- bdiag(t(A1),t(A2),t(A3),t(A4),t(A5))
  #FIs <- sapply(y,deriv_2_ll,p1=p1,m1=m1,m2=m2,sig2_1=sig2_1,sig2_2=sig2_2)
  #return(-matrix(apply(FIs,1,mean),5,5))
  return(list(FIM=FIM,B=B,fxhat=fxhat))
}


##############################################################################################################################
############## Select bandwidth based on AMISE for 10 simulations (n=200)
set.seed(389)
n <- 400
x <- cbind(runif(n),runif(n))  ## generate x from uniform distributions 
# x <- rmvnorm(n,rep(0,2),autocorr.mat(2,0.5))  ## generate x from normal distribution 
X_pi <- X_sig2 <- cbind(rep(1,n),x) ## linear approximation for proportion functions and variance functions 
X_m <- cbind(rep(1,n),x,x^2,x[,1]*x[,2])
p <- cbind(pi_1(x),pi_2(x)) 
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))

## define a function to perform bandwidth selection based on AMISE
bs_AMISE <- function(j){ 
  ## load libraries
  library(mixtools)
  library(fpc)
  library(ks)
  library(matrixcalc)
  library(sn)
  library(Matrix)
  
  y <- rnpnormmix(p,mu,sgm)
  rmt <- regmix(x,y,nclust=2)
  initl_beta_pi <- c(rmt$eps[1],0,0)
  initl_beta_sgm21 <- c(rmt$var[1],0,0)
  initl_beta_sgm22 <- c(rmt$var[2],0,0)
  initl_beta_m1 <- c(rmt$coef[,1],0,0,0)
  initl_beta_m2 <- c(rmt$coef[,2],0,0,0)
  parGMR <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
  sw_beta_pi <- -parGMR$beta_pi
  sw_beta_sgm21 <- parGMR$beta_sgm22
  sw_beta_sgm22 <- parGMR$beta_sgm21
  sw_beta_m1 <- parGMR$beta_m2
  sw_beta_m2 <- parGMR$beta_m1
  ptheta1_1 <- exp(-X_pi%*%sw_beta_pi)/(1+exp(-X_pi%*%sw_beta_pi))^2*sw_beta_pi[2]
  ptheta2_1 <- exp(X_sig2%*%sw_beta_sgm21)*sw_beta_sgm21[2]
  ptheta3_1 <- exp(X_sig2%*%sw_beta_sgm22)*sw_beta_sgm22[2]
  ptheta4_1 <- sw_beta_m1[2] + 2*sw_beta_m1[4]*x[,1] + sw_beta_m1[6]*x[,2]
  ptheta5_1 <- sw_beta_m2[2] + 2*sw_beta_m2[4]*x[,1] + sw_beta_m2[6]*x[,2]
  ptheta_1 <- cbind(x[,1],ptheta1_1,ptheta2_1,ptheta3_1,ptheta4_1,ptheta5_1)
  
  ptheta1_2 <- exp(-X_pi%*%sw_beta_pi)/(1+exp(-X_pi%*%sw_beta_pi))^2*sw_beta_pi[3]
  ptheta2_2 <- exp(X_sig2%*%sw_beta_sgm21)*sw_beta_sgm21[3]
  ptheta3_2 <- exp(X_sig2%*%sw_beta_sgm22)*sw_beta_sgm22[3]
  ptheta4_2 <- sw_beta_m1[3] + 2*sw_beta_m1[5]*x[,1] + sw_beta_m1[6]*x[,1]
  ptheta5_2 <- sw_beta_m2[3] + 2*sw_beta_m2[5]*x[,1] + sw_beta_m2[6]*x[,1]
  ptheta_2 <- cbind(x[,2],ptheta1_2,ptheta2_2,ptheta3_2,ptheta4_2,ptheta5_2)
  
  pi_intm <- exp(-X_pi%*%sw_beta_pi)
  ptheta1_11 <- sw_beta_pi[2]^2*pi_intm*(-1+pi_intm)/(1+pi_intm)^3
  ptheta2_11 <- exp(X_sig2%*%sw_beta_sgm21)*sw_beta_sgm21[2]^2
  ptheta3_11 <- exp(X_sig2%*%sw_beta_sgm22)*sw_beta_sgm22[2]^2
  ptheta4_11 <- rep(2*sw_beta_m1[4],nrow(x))
  ptheta5_11 <- rep(2*sw_beta_m2[4],nrow(x))
  ptheta_11 <- cbind(x[,1],ptheta1_11,ptheta2_11,ptheta3_11,ptheta4_11,ptheta5_11)
  
  ptheta1_22 <- sw_beta_pi[3]^2*pi_intm*(-1+pi_intm)/(1+pi_intm)^3
  ptheta2_22 <- exp(X_sig2%*%sw_beta_sgm21)*sw_beta_sgm21[3]^2
  ptheta3_22 <- exp(X_sig2%*%sw_beta_sgm22)*sw_beta_sgm22[3]^2
  ptheta4_22 <- rep(2*sw_beta_m1[5],nrow(x))
  ptheta5_22 <- rep(2*sw_beta_m2[5],nrow(x))
  ptheta_22 <- cbind(x[,2],ptheta1_22,ptheta2_22,ptheta3_22,ptheta4_22,ptheta5_22)
  
  ptheta1_12 <- sw_beta_pi[2]*sw_beta_pi[3]*pi_intm*(-1+pi_intm)/(1+pi_intm)^3
  ptheta2_12 <- exp(X_sig2%*%sw_beta_sgm21)*sw_beta_sgm21[2]*sw_beta_sgm21[3]
  ptheta3_12 <- exp(X_sig2%*%sw_beta_sgm22)*sw_beta_sgm21[2]*sw_beta_sgm22[3]
  ptheta4_12 <- rep(2*sw_beta_m1[6],nrow(x))
  ptheta5_12 <- rep(2*sw_beta_m2[6],nrow(x))
  ptheta_12 <- cbind(x[,1],ptheta1_12,ptheta2_12,ptheta3_12,ptheta4_12,ptheta5_12)
  
  thetas <- cbind(x,1/(1+pi_intm),exp(X_sig2%*%sw_beta_sgm21),exp(X_sig2%*%sw_beta_sgm22),X_m%*%sw_beta_m1,X_m%*%sw_beta_m2)
  
  kdde_0 <- kdde(x,deriv.order=0)
  fx <- cbind(x[,1],predict(kdde_0,x=x))
  kdde_1 <- kdde(x,deriv.order=1)
  grad_fx <- cbind(x[,1],predict(kdde_1,x=x))
  
  D <- kronecker(rep(1,5),diag(4))
  D_2 <- duplication.matrix(2)
  M_1 <- kronecker(diag(2),matrix(c(1,0,0,0),2))%*%D_2
  M_2 <- kronecker(diag(2),matrix(c(0,1,1,0),2))%*%D_2
  M_3 <- kronecker(diag(2),matrix(c(0,0,0,1),2))%*%D_2
  M <- cbind(M_1,M_2,M_3)
  G_1 <- vec(M%*%kronecker(diag(3),c(1,0,0))+kronecker(diag(2),matrix(c(1,0,0,0),2))%*%D_2)
  G_2 <- vec(M%*%kronecker(diag(3),c(0,1,0))+kronecker(diag(2),matrix(c(0,1,1,0),2))%*%D_2)
  G_3 <- vec(M%*%kronecker(diag(3),c(0,0,1))+kronecker(diag(2),matrix(c(0,0,0,1),2))%*%D_2)
  G <- cbind(G_1,G_2,G_3)
  
  c1 <- matrix(rep(0,16),nrow=4)
  c2 <- 0
  for(i in 1:n){
    temp_para <- mse_pt_normmix(x0=x[i,],thetas=thetas,ptheta_1=ptheta_1,ptheta_2=ptheta_2,ptheta_11=ptheta_11,
                                ptheta_22=ptheta_22,ptheta_12=ptheta_12,fx=fx,grad_fx=grad_fx)
    c1 <- c1 + t(temp_para$B%*%D)%*%solve(temp_para$FIM)%*%solve(temp_para$FIM)%*%temp_para$B%*%D
    c2 <- c2 + tr(solve(temp_para$FIM))/temp_para$fxhat
  }
  c1 <- matrix(c1,4)
  newton.raphson_bw <- function(eta_int,c1,c2,tol=1e-6,nmax=1000){
    eta_old <- eta_int
    
    for(i in 1:nmax){
      vech_inv <- vech2mat(c(eta_old))
      par1 <- t(M%*%kronecker(diag(3),eta_old)+kronecker(diag(2),vech_inv)%*%D_2)%*%c1%*%kronecker(diag(2),vech_inv)%*%D_2%*%eta_old/2 -
        1*c2/(4*pi*n)*t(D_2)%*%vec(solve(vech_inv))/det(vech_inv)
      par2 <- kronecker(diag(3),t(c1%*%kronecker(diag(2),vech_inv)%*%D_2%*%eta_old))%*%G/2 +
        t(M%*%kronecker(diag(3),eta_old)+kronecker(diag(2),vech_inv)%*%D_2)%*%c1%*%(M%*%kronecker(diag(3),eta_old)+kronecker(diag(2),vech_inv)%*%D_2)/2+
        1*c2/(8*pi*n*det(vech_inv))*(t(D_2)%*%t(kronecker(solve(vech_inv),diag(2)))%*%(vec(diag(2))%*%t(vec(diag(2)))+2*diag(4))%*%(kronecker(diag(2),solve(vech_inv)))%*%D_2 +
                                       t(D_2)%*%vec(solve(vech_inv))%*%t(vec(solve(vech_inv)))%*%D_2)
      
      eta_new <- eta_old - solve(par2)%*%par1
      if(max(abs(eta_new-eta_old))<tol)
      {return(eta_new)}
      eta_old <- eta_new
      print(i)
      print(eta_old)
    }
    print("Need more iterations")
  }
  
  eta_int <- c(1,0,1)
  eta_opt <- newton.raphson_bw(eta_int,c1,c2)
  H_opt <- vech2mat(c(eta_opt))%*%vech2mat(c(eta_opt))
  return(H_opt)
}
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

#################################################################################################################################################
##### define a function to implement the modified EM agorithm for two component nonparametric mixture of regressions based on 100 grid points
##### N is the number of grid points, initl_m, initl_sgm,initl_pi are n by c (c=2)
npnormmixEM_grid <- function(x,y,N=100,initl_m,initl_sgm,initl_pi,maxiter=1000,epsilon=0.00001,bw){
  u <- as.matrix(expand.grid(seq(min(x[,1]),max(x[,1]),l=10),seq(min(x[,2]),max(x[,2]),l=10)))
  n <- length(y)
  l <- 0
  m_EM <- initl_m
  sgm_EM <- initl_sgm
  pi_EM <- initl_pi
  theta_old <- rep(0,6*N)
  theta_diff <- epsilon 
  
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    r <- pi_EM*dnorm(cbind(y,y),mean=m_EM,sd=sgm_EM)/rowSums(pi_EM*dnorm(cbind(y,y),mean=m_EM,sd=sgm_EM))
    pi1_grid <- rep(0,N)
    pi2_grid <- rep(0,N)
    m1_grid <- rep(0,N)
    m2_grid <- rep(0,N)
    sgm2_1_grid <- rep(0,N)
    sgm2_2_grid <- rep(0,N)
    
    for(j in 1:N){
      temp_u <- matrix(rep(u[j,],n),nrow=n,byrow=TRUE)
      #temp_EpaKernel <- EpaKernel(x-temp_u,bw) 
      temp_EpaKernel <- NormKernel(x-temp_u,bw)
      pi1_grid[j] <- r[,1]%*%temp_EpaKernel/sum(temp_EpaKernel)
      pi2_grid[j] <- 1-pi1_grid[j]
      m1_grid[j] <- (r[,1]*temp_EpaKernel)%*%y/(r[,1]%*%temp_EpaKernel)
      m2_grid[j] <- (r[,2]*temp_EpaKernel)%*%y/(r[,2]%*%temp_EpaKernel)
      sgm2_1_grid[j] <- (r[,1]*temp_EpaKernel)%*%((y-m1_grid[j])^2)/(r[,1]%*%temp_EpaKernel)
      sgm2_2_grid[j] <- (r[,2]*temp_EpaKernel)%*%((y-m2_grid[j])^2)/(r[,2]%*%temp_EpaKernel)
    }
    
    theta_new <- c(pi1_grid,pi2_grid,m1_grid,m2_grid,sgm2_1_grid,sgm2_2_grid)
    theta_diff <- max(abs(theta_new - theta_old))
    theta_old <- theta_new
    ### interpolation from grid points
    inter_pi1 <- interpp(x=u[,1],y=u[,2],z=pi1_grid,xo=x[,1],yo=x[,2],linear=FALSE,extrap=TRUE)$z
    if (anyNA(inter_pi1)){
      inter_pi1_narm <- missForest(cbind(x,inter_pi1))$ximp
      inter_pi1 <- inter_pi1_narm[,3]
    }
    inter_pi2 <- 1 - inter_pi1
    inter_m1 <- interpp(x=u[,1],y=u[,2],z=m1_grid,xo=x[,1],yo=x[,2],linear=FALSE,extrap=TRUE)$z
    if (anyNA(inter_m1)){
      inter_m1_narm <- missForest(cbind(x,inter_m1))$ximp
      inter_m1 <- inter_m1_narm[,3]
    }
    inter_m2 <- interpp(x=u[,1],y=u[,2],z=m2_grid,xo=x[,1],yo=x[,2],linear=FALSE,extrap=TRUE)$z
    if (anyNA(inter_m2)){
      inter_m2_narm <- missForest(cbind(x,inter_m2))$ximp
      inter_m2 <- inter_m2_narm[,3]
    }
    inter_sgm2_1 <- interpp(x=u[,1],y=u[,2],z=sgm2_1_grid,xo=x[,1],yo=x[,2],linear=FALSE,extrap=TRUE)$z
    if (anyNA(inter_sgm2_1)){
      inter_sgm2_1_narm <- missForest(cbind(x,inter_sgm2_1))$ximp
      inter_sgm2_1 <- inter_sgm2_1_narm[,3]
    }
    inter_sgm2_2 <- interpp(x=u[,1],y=u[,2],z=sgm2_2_grid,xo=x[,1],yo=x[,2],linear=FALSE,extrap=TRUE)$z
    if (anyNA(inter_sgm2_2)){
      inter_sgm2_2_narm <- missForest(cbind(x,inter_sgm2_2))$ximp
      inter_sgm2_2 <- inter_sgm2_2_narm[,3]
    }
    
    pi_EM[,1] <- inter_pi1
    pi_EM[,2] <- inter_pi2
    m_EM[,1] <- inter_m1
    m_EM[,2] <- inter_m2
    sgm_EM[,1] <- sqrt(inter_sgm2_1)
    sgm_EM[,2] <- sqrt(inter_sgm2_2)
  }
  if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(pi1=inter_pi1,pi2=inter_pi2,m1=inter_m1,m2=inter_m2,sigma2_1=inter_sgm2_1,sigma2_2=inter_sgm2_2,iter=l))
}

NormKernel <- function(u,h){
  n_u <- nrow(u)
  outcome <- rep(0,n_u)
  for (i in 1:n_u){
    outcome[i] <- dmvnorm(u[i,],sigma=h)
  }
  return(outcome)
}

##############################################################################################################################
############## Conduct simulation studies (500 samples, n=200)
set.seed(579)
initl_m <- cbind(m_1(x)+rnorm(1,sd=0.5),m_2(x)+rnorm(1,sd=0.5))  ### set initial values for the mean function 
initl_sgm <- cbind(rep(mean(sigma_1(x)),n),rep(mean(sigma_2(x)),n))   ## set the initial values for the standard deviation function
initl_pi <- cbind(rep(mean(pi_1(x)),n),rep(mean(pi_2(x)),n))  ## set the initial values for the mixing proportion function

## define a function for parallel simulation
sim_npmixreg <- function(i){
  # load libraries
  library(mixtools)
  library(missForest) ## library for imputation
  library(interp)
  library(fpc)
  y <- rnpnormmix(p,mu,sgm)
  nprmt <- npnormmixEM_grid(x,y,N=100,initl_m=initl_m,initl_sgm=initl_sgm,initl_pi=initl_pi,bw=h_AMISE) 
  RASE_pi_npmixreg <- sqrt((sum((nprmt$pi1-pi_1(x))^2)+sum((nprmt$pi2-pi_2(x))^2))/n)
  RASE_m_npmixreg  <- sqrt((sum((nprmt$m1-m_1(x))^2)+sum((nprmt$m2-m_2(x))^2))/n)
  RASE_sigma2_npmixreg <- sqrt((sum((nprmt$sigma2_1-sigma_1(x)^2)^2)+sum((nprmt$sigma2_2-sigma_2(x)^2)^2))/n)
  
  rmt <- regmix(x,y,nclust=1:2)
  xd <- cbind(rep(1,n),x)
  if(mean(xd%*%rmt$coef[,1]) > mean(xd%*%rmt$coef[,2])){
    RASE_pi_mixreg <- sqrt((sum((rep(rmt$eps[1],n)-pi_1(x))^2)+sum((rmt$pi2-pi_2(x))^2))/n)
    RASE_m_mixreg  <- sqrt((sum((xd%*%rmt$coef[,1]-m_1(x))^2)+sum((xd%*%rmt$coef[,2]-m_2(x))^2))/n)
    RASE_sigma2_mixreg <- sqrt((sum((rep(rmt$var[1],n)-sigma_1(x)^2)^2)+sum((rep(rmt$var[2],n)-sigma_2(x)^2)^2))/n)
  } else
  {
    RASE_pi_mixreg <- sqrt((sum((rep(rmt$eps[2],n)-pi_1(x))^2)+sum((rep(rmt$eps[1],n)-pi_2(x))^2))/n)
    RASE_m_mixreg  <- sqrt((sum((xd%*%rmt$coef[,2]-m_1(x))^2)+sum((xd%*%rmt$coef[,1]-m_2(x))^2))/n)
    RASE_sigma2_mixreg <- sqrt((sum((rep(rmt$var[2],n)-sigma_1(x)^2)^2)+sum((rep(rmt$var[1],n)-sigma_2(x)^2)^2))/n)
  }
  return(list(RASE_pi_npmixreg=RASE_pi_npmixreg,RASE_m_npmixreg=RASE_m_npmixreg, RASE_sigma2_npmixreg= RASE_sigma2_npmixreg,
              RASE_pi_mixreg=RASE_pi_mixreg,RASE_m_mixreg=RASE_m_mixreg,RASE_sigma2_mixreg=RASE_sigma2_mixreg))
}

n_sim <- 500
the_cluster <- makeCluster(16)
clusterSetRNGStream(the_cluster,27)
clusterExport(the_cluster,c("npnormmixEM_grid","NormKernel","rnpnormmix","x","n","h_AMISE","initl_m","initl_sgm","initl_pi","pi_1",
                            "pi_2","m_1","m_2","sigma_1","sigma_2","p","mu","sgm"))
test_sim <- parLapply(the_cluster,1:n_sim,sim_npmixreg)
stopCluster(the_cluster)

RASEs <- matrix(unlist(test_sim),nrow=6,byrow = FALSE)
RASE_mean <- apply(RASEs,1,mean)
RASE_sd <- apply(RASEs,1,sd)






