

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

n <- 200
h_grid <- c(0.05,0.08,0.11,0.14,0.17)
k_fold <- 5
h_opt <- rep(0,10)
set.seed(893)
for(q in 1:10){
  x <- runif(n)
  p <- cbind(pi_1(x),pi_2(x))
  mu <- cbind(m_1(x),m_2(x))
  sgm <- cbind(sigma_1(x),sigma_2(x))
  y <- rnpnormmix(p,mu,sgm)
  ## generate index for test and training sets
  ind_cv <- sample(x = rep(1:k_fold, each = n / k_fold), size = n)
  cv_error <- matrix(rep(0,length(h_grid)*k_fold),ncol=k_fold)
  
  for(i_para in 1:length(h_grid)){
    
    for(i_test in 1:k_fold){
      
      y_test <- y[ind_cv==i_test] ; y_train <- y[ind_cv!=i_test]
      x_test <- x[ind_cv==i_test] ; x_train <- x[ind_cv!=i_test]
      n_test <- length(y_test) ; n_train <- length(y_train)
      ## set initial values based on true functions
      initl_m_train <- cbind(m_1(x_train)+rnorm(1,sd=0.5),m_2(x_train)+rnorm(1,sd=0.5))  ### set initial values for the mean function 
      initl_sgm_train <- cbind(rep(mean(sigma_1(x_train)),n_train),rep(mean(sigma_2(x_train)),n_train))   ## set the initial values for the standard deviation function
      initl_pi_train <- cbind(rep(mean(pi_1(x_train)),n_train),rep(mean(pi_2(x_train)),n_train))  ## set the initial values for the mixing proportion function
      out <- npnormmixEM(x_train,y_train,initl_m=initl_m_train,initl_sgm=initl_sgm_train,initl_pi=initl_pi_train,bw=h_grid[i_para]) 
      ## interpolation from the training set 
      ### define functions for interpolation from training set in CV
      
      inter_m <- function(x_train,r_train,y_train,x,h){
        return(sum(r_train*EpaKernel(x_train-x,h)*y_train)/sum(r_train*EpaKernel(x_train-x,h)))
      }
      
      inter_sigma2 <- function(x_train,r_train,y_train,m_train,x,h){
        return(sum(r_train*EpaKernel(x_train-x,h)*(y_train-m_train)^2)/sum(r_train*EpaKernel(x_train-x,h)))
      }
      
      inter_pi <- function(x_train,r_train,x,h){
        return(sum(r_train*EpaKernel(x_train-x,h))/sum(EpaKernel(x_train-x,h)))
      }
      
      r1_train <- out$pi1*dnorm(y_train,mean=out$m1,sd=sqrt(out$sigma2_1))/(out$pi1*dnorm(y_train,mean=out$m1,sd=sqrt(out$sigma2_1))+
                                                                              out$pi2*dnorm(y_train,mean=out$m2,sd=sqrt(out$sigma2_2)))
      r2_train <- 1 - r1_train
      p1_test <- sapply(x_test,inter_pi,x_train=x_train,r_train=r1_train,h=h_grid[i_para])
      p2_test <- 1 - pi1_test 
      m1_test <- sapply(x_test,inter_m,x_train=x_train,r_train=r1_train,y_train=y_train,h=h_grid[i_para])
      m2_test <- sapply(x_test,inter_m,x_train=x_train,r_train=r1_train,y_train=y_train,h=h_grid[i_para])
      sigma2_1_test <- sapply(x_test,inter_sigma2,x_train=x_train,r_train=r1_train,y_train=y_train,m_train=out$m1,h=h_grid[i_para])
      sigma2_2_test <- sapply(x_test,inter_sigma2,x_train=x_train,r_train=r1_train,y_train=y_train,m_train=out$m2,h=h_grid[i_para])
      
      
      r1_test <- p1_test*dnorm(y_test,mean=m1_test,sd=sqrt(sigma2_1_test))/(p1_test*dnorm(y_test,mean=m1_test,sd=sqrt(sigma2_1_test))+p2_test*dnorm(y_test,mean=m2_test,sd=sqrt(sigma2_2_test)))
      r2_test <- 1 - r1_test
      y_pred <- r1_test*m1_test + r2_test*m2_test
      
      cv_error[i_para,i_test] <- sum((y_pred-y_test)^2)
    }
  }
  
  mincv <- min(rowMeans(cv_error,na.rm=TRUE),na.rm = TRUE)
  indexcv <- rowMeans(cv_error,na.rm=TRUE)== mincv
  indexcv[is.na(indexcv)] <- FALSE
  
  h_opt[q] <- h_grid[indexcv]
  print(q)
}
print(mean(h_opt))  ## h_opt = 0.11


###########################################  n =400   ###############################################################
n <- 400
h_grid <- c(0.05,0.08,0.11,0.14,0.17)
k_fold <- 5
h_opt <- rep(0,10)
set.seed(216)
for(q in 1:10){
  x <- runif(n)
  p <- cbind(pi_1(x),pi_2(x))
  mu <- cbind(m_1(x),m_2(x))
  sgm <- cbind(sigma_1(x),sigma_2(x))
  y <- rnpnormmix(p,mu,sgm)
  ## generate index for test and training sets
  ind_cv <- sample(x = rep(1:k_fold, each = n / k_fold), size = n)
  cv_error <- matrix(rep(0,length(h_grid)*k_fold),ncol=k_fold)
  
  for(i_para in 1:length(h_grid)){
    
    for(i_test in 1:k_fold){
      
      y_test <- y[ind_cv==i_test] ; y_train <- y[ind_cv!=i_test]
      x_test <- x[ind_cv==i_test] ; x_train <- x[ind_cv!=i_test]
      n_test <- length(y_test) ; n_train <- length(y_train)
      ## set initial values based on true functions
      initl_m_train <- cbind(m_1(x_train)+rnorm(1,sd=0.5),m_2(x_train)+rnorm(1,sd=0.5))  ### set initial values for the mean function 
      initl_sgm_train <- cbind(rep(mean(sigma_1(x_train)),n_train),rep(mean(sigma_2(x_train)),n_train))   ## set the initial values for the standard deviation function
      initl_pi_train <- cbind(rep(mean(pi_1(x_train)),n_train),rep(mean(pi_2(x_train)),n_train))  ## set the initial values for the mixing proportion function
      out <- npnormmixEM(x_train,y_train,initl_m=initl_m_train,initl_sgm=initl_sgm_train,initl_pi=initl_pi_train,bw=h_grid[i_para]) 
      ## interpolation from the training set 
      ### define functions for interpolation from training set in CV
      
      inter_m <- function(x_train,r_train,y_train,x,h){
        return(sum(r_train*EpaKernel(x_train-x,h)*y_train)/sum(r_train*EpaKernel(x_train-x,h)))
      }
      
      inter_sigma2 <- function(x_train,r_train,y_train,m_train,x,h){
        return(sum(r_train*EpaKernel(x_train-x,h)*(y_train-m_train)^2)/sum(r_train*EpaKernel(x_train-x,h)))
      }
      
      inter_pi <- function(x_train,r_train,x,h){
        return(sum(r_train*EpaKernel(x_train-x,h))/sum(EpaKernel(x_train-x,h)))
      }
      
      r1_train <- out$pi1*dnorm(y_train,mean=out$m1,sd=sqrt(out$sigma2_1))/(out$pi1*dnorm(y_train,mean=out$m1,sd=sqrt(out$sigma2_1))+
                                                                              out$pi2*dnorm(y_train,mean=out$m2,sd=sqrt(out$sigma2_2)))
      r2_train <- 1 - r1_train
      p1_test <- sapply(x_test,inter_pi,x_train=x_train,r_train=r1_train,h=h_grid[i_para])
      p2_test <- 1 - pi1_test 
      m1_test <- sapply(x_test,inter_m,x_train=x_train,r_train=r1_train,y_train=y_train,h=h_grid[i_para])
      m2_test <- sapply(x_test,inter_m,x_train=x_train,r_train=r1_train,y_train=y_train,h=h_grid[i_para])
      sigma2_1_test <- sapply(x_test,inter_sigma2,x_train=x_train,r_train=r1_train,y_train=y_train,m_train=out$m1,h=h_grid[i_para])
      sigma2_2_test <- sapply(x_test,inter_sigma2,x_train=x_train,r_train=r1_train,y_train=y_train,m_train=out$m2,h=h_grid[i_para])
      
      
      r1_test <- p1_test*dnorm(y_test,mean=m1_test,sd=sqrt(sigma2_1_test))/(p1_test*dnorm(y_test,mean=m1_test,sd=sqrt(sigma2_1_test))+p2_test*dnorm(y_test,mean=m2_test,sd=sqrt(sigma2_2_test)))
      r2_test <- 1 - r1_test
      y_pred <- r1_test*m1_test + r2_test*m2_test
      
      cv_error[i_para,i_test] <- sum((y_pred-y_test)^2)
    }
  }
  
  mincv <- min(rowMeans(cv_error,na.rm=TRUE),na.rm = TRUE)
  indexcv <- rowMeans(cv_error,na.rm=TRUE)== mincv
  indexcv[is.na(indexcv)] <- FALSE
  
  h_opt[q] <- h_grid[indexcv]
  print(q)
}
print(mean(h_opt))  ## h_opt = 0.083

###########################################  n =800   ###############################################################
n <- 800
h_grid <- c(0.05,0.08,0.11,0.14,0.17)
k_fold <- 5
h_opt <- rep(0,10)
set.seed(893)
for(q in 1:10){
  x <- runif(n)
  p <- cbind(pi_1(x),pi_2(x))
  mu <- cbind(m_1(x),m_2(x))
  sgm <- cbind(sigma_1(x),sigma_2(x))
  y <- rnpnormmix(p,mu,sgm)
  ## generate index for test and training sets
  ind_cv <- sample(x = rep(1:k_fold, each = n / k_fold), size = n)
  cv_error <- matrix(rep(0,length(h_grid)*k_fold),ncol=k_fold)
  
  for(i_para in 1:length(h_grid)){
    
    for(i_test in 1:k_fold){
      
      y_test <- y[ind_cv==i_test] ; y_train <- y[ind_cv!=i_test]
      x_test <- x[ind_cv==i_test] ; x_train <- x[ind_cv!=i_test]
      n_test <- length(y_test) ; n_train <- length(y_train)
      ## set initial values based on true functions
      initl_m_train <- cbind(m_1(x_train)+rnorm(1,sd=0.5),m_2(x_train)+rnorm(1,sd=0.5))  ### set initial values for the mean function 
      initl_sgm_train <- cbind(rep(mean(sigma_1(x_train)),n_train),rep(mean(sigma_2(x_train)),n_train))   ## set the initial values for the standard deviation function
      initl_pi_train <- cbind(rep(mean(pi_1(x_train)),n_train),rep(mean(pi_2(x_train)),n_train))  ## set the initial values for the mixing proportion function
      out <- npnormmixEM(x_train,y_train,initl_m=initl_m_train,initl_sgm=initl_sgm_train,initl_pi=initl_pi_train,bw=h_grid[i_para]) 
      ## interpolation from the training set 
      ### define functions for interpolation from training set in CV
      
      inter_m <- function(x_train,r_train,y_train,x,h){
        return(sum(r_train*EpaKernel(x_train-x,h)*y_train)/sum(r_train*EpaKernel(x_train-x,h)))
      }
      
      inter_sigma2 <- function(x_train,r_train,y_train,m_train,x,h){
        return(sum(r_train*EpaKernel(x_train-x,h)*(y_train-m_train)^2)/sum(r_train*EpaKernel(x_train-x,h)))
      }
      
      inter_pi <- function(x_train,r_train,x,h){
        return(sum(r_train*EpaKernel(x_train-x,h))/sum(EpaKernel(x_train-x,h)))
      }
      
      r1_train <- out$pi1*dnorm(y_train,mean=out$m1,sd=sqrt(out$sigma2_1))/(out$pi1*dnorm(y_train,mean=out$m1,sd=sqrt(out$sigma2_1))+
                                                                              out$pi2*dnorm(y_train,mean=out$m2,sd=sqrt(out$sigma2_2)))
      r2_train <- 1 - r1_train
      p1_test <- sapply(x_test,inter_pi,x_train=x_train,r_train=r1_train,h=h_grid[i_para])
      p2_test <- 1 - pi1_test 
      m1_test <- sapply(x_test,inter_m,x_train=x_train,r_train=r1_train,y_train=y_train,h=h_grid[i_para])
      m2_test <- sapply(x_test,inter_m,x_train=x_train,r_train=r1_train,y_train=y_train,h=h_grid[i_para])
      sigma2_1_test <- sapply(x_test,inter_sigma2,x_train=x_train,r_train=r1_train,y_train=y_train,m_train=out$m1,h=h_grid[i_para])
      sigma2_2_test <- sapply(x_test,inter_sigma2,x_train=x_train,r_train=r1_train,y_train=y_train,m_train=out$m2,h=h_grid[i_para])
      
      
      r1_test <- p1_test*dnorm(y_test,mean=m1_test,sd=sqrt(sigma2_1_test))/(p1_test*dnorm(y_test,mean=m1_test,sd=sqrt(sigma2_1_test))+p2_test*dnorm(y_test,mean=m2_test,sd=sqrt(sigma2_2_test)))
      r2_test <- 1 - r1_test
      y_pred <- r1_test*m1_test + r2_test*m2_test
      
      cv_error[i_para,i_test] <- sum((y_pred-y_test)^2)
    }
  }
  
  mincv <- min(rowMeans(cv_error,na.rm=TRUE),na.rm = TRUE)
  indexcv <- rowMeans(cv_error,na.rm=TRUE)== mincv
  indexcv[is.na(indexcv)] <- FALSE
  
  h_opt[q] <- h_grid[indexcv]
  print(q)
}
print(mean(h_opt))  ## h_opt = 0.077  

