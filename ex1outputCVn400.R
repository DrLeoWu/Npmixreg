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
  return(exp(0.5*x[,1])/(1+exp(0.5*x[,2]))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(4+7*dmvnorm(2*x-1,sigma = 0.3*diag(2)))
}
m_2 <- function(x){
  return(2-10*dmvnorm(2*x-1,sigma = 0.5*diag(2)))
}

sigma_1 <- function(x){
  return(0.6*exp(0.1*x[,1]+0.4*x[,2]))
}

sigma_2 <- function(x){
  return(0.5*exp(-0.1*x[,1]-0.2*x[,2]))
}

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
      x_na <- x[is.na(inter_pi1),]
      x_nna <- x[!is.na(inter_pi1),]
      pi1_nna <- inter_pi1[!is.na(inter_pi1)]
      if(is.matrix(x_na)){
        for(i_na in nrow(x_na)){
          dis <- rowSums((x_nna - x_na[i_na,])^2)
          if(x_na[i_na,1]==min(x[,1])||x_na[i_na,1]==max(x[,1])||x_na[i_na,2]==min(x[,2])||x_na[i_na,2]==max(x[,2])){
            inter_pi1[which(is.na(inter_pi1))[i_na]] <- pi1_nna[dis==min(dis)][1]
          }
        }
      } else {dis <- rowSums((x_nna - x_na)^2)
      if(x_na[1]==min(x[,1])||x_na[1]==max(x[,1])||x_na[2]==min(x[,2])||x_na[2]==max(x[,2])){
        inter_pi1[which(is.na(inter_pi1))] <- pi1_nna[dis==min(dis)][1]}}
      inter_pi1_narm <- missForest(cbind(x,inter_pi1))$ximp
      inter_pi1 <- inter_pi1_narm[,3]
    }
    inter_pi2 <- 1 - inter_pi1
    inter_m1 <- interpp(x=u[,1],y=u[,2],z=m1_grid,xo=x[,1],yo=x[,2],linear=FALSE,extrap=TRUE)$z
    if (anyNA(inter_m1)){
      x_na <- x[is.na(inter_m1),]
      x_nna <- x[!is.na(inter_m1),]
      m1_nna <- inter_m1[!is.na(inter_m1)]
      if(is.matrix(x_na)){
        for(i_na in nrow(x_na)){
          dis <- rowSums((x_nna - x_na[i_na,])^2)
          if(x_na[i_na,1]==min(x[,1])||x_na[i_na,1]==max(x[,1])||x_na[i_na,2]==min(x[,2])||x_na[i_na,2]==max(x[,2])){
            inter_m1[which(is.na(inter_m1))[i_na]] <- m1_nna[dis==min(dis)][1]
          }
        }
      } else {dis <- rowSums((x_nna - x_na)^2)
      if(x_na[1]==min(x[,1])||x_na[1]==max(x[,1])||x_na[2]==min(x[,2])||x_na[2]==max(x[,2])){
        inter_m1[which(is.na(inter_m1))] <- m1_nna[dis==min(dis)][1]}}
      
      inter_m1_narm <- missForest(cbind(x,inter_m1))$ximp
      inter_m1 <- inter_m1_narm[,3]
    }
    inter_m2 <- interpp(x=u[,1],y=u[,2],z=m2_grid,xo=x[,1],yo=x[,2],linear=FALSE,extrap=TRUE)$z
    if (anyNA(inter_m2)){
      x_na <- x[is.na(inter_m2),]
      x_nna <- x[!is.na(inter_m2),]
      m2_nna <- inter_m2[!is.na(inter_m2)]
      if(is.matrix(x_na)){
        for(i_na in nrow(x_na)){
          dis <- rowSums((x_nna - x_na[i_na,])^2)
          if(x_na[i_na,1]==min(x[,1])||x_na[i_na,1]==max(x[,1])||x_na[i_na,2]==min(x[,2])||x_na[i_na,2]==max(x[,2])){
            inter_m2[which(is.na(inter_m2))[i_na]] <- m2_nna[dis==min(dis)][1]
          }
        }
      } else {dis <- rowSums((x_nna - x_na)^2)
      if(x_na[1]==min(x[,1])||x_na[1]==max(x[,1])||x_na[2]==min(x[,2])||x_na[2]==max(x[,2])){
        inter_m2[which(is.na(inter_m2))] <- m2_nna[dis==min(dis)][1]}}
      
      inter_m2_narm <- missForest(cbind(x,inter_m2))$ximp
      inter_m2 <- inter_m2_narm[,3]
    }
    inter_sgm2_1 <- interpp(x=u[,1],y=u[,2],z=sgm2_1_grid,xo=x[,1],yo=x[,2],linear=FALSE,extrap=TRUE)$z
    if (anyNA(inter_sgm2_1)){
      x_na <- x[is.na(inter_sgm2_1),]
      x_nna <- x[!is.na(inter_sgm2_1),]
      sgm2_1_nna <- inter_sgm2_1[!is.na(inter_sgm2_1)]
      if(is.matrix(x_na)){
        for(i_na in nrow(x_na)){
          dis <- rowSums((x_nna - x_na[i_na,])^2)
          if(x_na[i_na,1]==min(x[,1])||x_na[i_na,1]==max(x[,1])||x_na[i_na,2]==min(x[,2])||x_na[i_na,2]==max(x[,2])){
            inter_sgm2_1[which(is.na(inter_sgm2_1))[i_na]] <- sgm2_1_nna[dis==min(dis)][1]
          }
        }
      } else {dis <- rowSums((x_nna - x_na)^2)
      if(x_na[1]==min(x[,1])||x_na[1]==max(x[,1])||x_na[2]==min(x[,2])||x_na[2]==max(x[,2])){
        inter_sgm2_1[which(is.na(inter_sgm2_1))] <- sgm2_1_nna[dis==min(dis)][1]}}
      
      inter_sgm2_1_narm <- missForest(cbind(x,inter_sgm2_1))$ximp
      inter_sgm2_1 <- inter_sgm2_1_narm[,3]
    }
    inter_sgm2_2 <- interpp(x=u[,1],y=u[,2],z=sgm2_2_grid,xo=x[,1],yo=x[,2],linear=FALSE,extrap=TRUE)$z
    if (anyNA(inter_sgm2_2)){
      x_na <- x[is.na(inter_sgm2_2),]
      x_nna <- x[!is.na(inter_sgm2_2),]
      sgm2_2_nna <- inter_sgm2_2[!is.na(inter_sgm2_2)]
      if(is.matrix(x_na)){
        for(i_na in nrow(x_na)){
          dis <- rowSums((x_nna - x_na[i_na,])^2)
          if(x_na[i_na,1]==min(x[,1])||x_na[i_na,1]==max(x[,1])||x_na[i_na,2]==min(x[,2])||x_na[i_na,2]==max(x[,2])){
            inter_sgm2_2[which(is.na(inter_sgm2_2))[i_na]] <- sgm2_2_nna[dis==min(dis)][1]
          }
        }
      } else {dis <- rowSums((x_nna - x_na)^2)
      if(x_na[1]==min(x[,1])||x_na[1]==max(x[,1])||x_na[2]==min(x[,2])||x_na[2]==max(x[,2])){
        inter_sgm2_2[which(is.na(inter_sgm2_2))] <- sgm2_2_nna[dis==min(dis)][1]}}
      
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
############## Select bandwidth based on k-fold cross validation (assuming diagonal bandwidth matrix) for 10 simulations (n=200)
set.seed(389)
n <- 400
x <- cbind(runif(n),runif(n))  ## generate x from uniform distributions 
# x <- rmvnorm(n,rep(0,2),autocorr.mat(2,0.5))  ## generate x from normal distribution 
p <- cbind(pi_1(x),pi_2(x)) 
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))

#### define a function to conduct k-fold cross validation 
#h_11 <- h_22 <-seq(0.01,0.03,0.02)
h_11 <- h_22 <- c(0.004,0.009,0.014,0.019)
h_diag <- expand.grid(h_11,h_22)
k_fold <- 5

bs_CV_diag <- function(j){
  # load libraries
  library(mixtools)
  library(missForest) ## library for imputation
  library(interp)
  y <- rnpnormmix(p,mu,sgm)
  ## generate index for test and training sets
  ind_cv <- sample(x = rep(1:k_fold, each = n / k_fold), size = n)
  cv_error <- matrix(rep(0,nrow(h_diag)*k_fold),ncol=k_fold)
  
  for(i_para in 1:nrow(h_diag)){
    
    for(i_test in 1:k_fold){
      
      y_test <- y[ind_cv==i_test] ; y_train <- y[ind_cv!=i_test]
      x_test <- x[ind_cv==i_test,] ; x_train <- x[ind_cv!=i_test,]
      n_test <- length(y_test) ; n_train <- length(y_train)
      ## set initial values based on true functions
      initl_m_train <- cbind(m_1(x_train)+rnorm(1,sd=0.5),m_2(x_train)+rnorm(1,sd=0.5))  ### set initial values for the mean function 
      initl_sgm_train <- cbind(rep(mean(sigma_1(x_train)),n_train),rep(mean(sigma_2(x_train)),n_train))   ## set the initial values for the standard deviation function
      initl_pi_train <- cbind(rep(mean(pi_1(x_train)),n_train),rep(mean(pi_2(x_train)),n_train))  ## set the initial values for the mixing proportion function
      model_train <- npnormmixEM_grid(x_train,y_train,initl_m=initl_m_train,initl_sgm=initl_sgm_train,initl_pi=initl_pi_train,bw=diag(h_diag[i_para,])) 
      ## interpolation from the training set 
      pi1_test <- interpp(x=x_train[,1],y=x_train[,2],z=model_train$pi1,xo=x_test[,1],yo=x_test[,2],linear=FALSE,extrap=TRUE)$z
      if (anyNA(pi1_test)){
        y_test <- y_test[!is.na(pi1_test)] ; x_test <- x_test[!is.na(pi1_test),]
        pi1_test <- pi1_test[!is.na(pi1_test)]
      }
      if (max(pi1_test)>1){ pi1_test[pi1_test>1] <- 1}
      if (min(pi1_test)<0){ pi1_test[pi1_test<0] <- 0}
      pi2_test <- 1 - pi1_test
      m1_test <- interpp(x=x_train[,1],y=x_train[,2],z=model_train$m1,xo=x_test[,1],yo=x_test[,2],linear=FALSE,extrap=TRUE)$z
      if (anyNA(m1_test)){
        y_test <- y_test[!is.na(m1_test)] ; x_test <- x_test[!is.na(m1_test),]
        m1_test <- m1_test[!is.na(m1_test)]
      }
      m2_test <- interpp(x=x_train[,1],y=x_train[,2],z=model_train$m2,xo=x_test[,1],yo=x_test[,2],linear=FALSE,extrap=TRUE)$z
      if (anyNA(m2_test)){
        y_test <- y_test[!is.na(m2_test)] ; x_test <- x_test[!is.na(m2_test),]
        m2_test <- m2_test[!is.na(m2_test)]
      }
      sgm2_1_test <- interpp(x=x_train[,1],y=x_train[,2],z=model_train$sigma2_1,xo=x_test[,1],yo=x_test[,2],linear=FALSE,extrap=TRUE)$z
      if (anyNA(sgm2_1_test)){
        y_test <- y_test[!is.na(sgm2_1_test)] ; x_test <- x_test[!is.na(sgm2_1_test),]
        sgm2_1_test <- sgm2_1_test[!is.na(sgm2_1_test)]
      }
      sgm2_2_test <- interpp(x=x_train[,1],y=x_train[,2],z=model_train$sigma2_2,xo=x_test[,1],yo=x_test[,2],linear=FALSE,extrap=TRUE)$z
      if (anyNA(sgm2_2_test)){
        y_test <- y_test[!is.na(sgm2_2_test)] ; x_test <- x_test[!is.na(sgm2_2_test),]
        sgm2_2_test <- sgm2_2_test[!is.na(sgm2_2_test)]
      }
      
      r1_test <- pi1_test*dnorm(y_test,mean=m1_test,sd=sqrt(sgm2_1_test))/(pi1_test*dnorm(y_test,mean=m1_test,sd=sqrt(sgm2_1_test))+pi2_test*dnorm(y_test,mean=m2_test,sd=sqrt(sgm2_2_test)))
      r2_test <- 1 - r1_test
      cv_error[i_para,i_test] <- sum((y_test - (r1_test*m1_test+r2_test*m2_test))^2)/length(y_test)
      
    }
    print(i_para)
  }
  mincv <- min(rowMeans(cv_error,na.rm=TRUE),na.rm = TRUE)
  #print(rowMeans(cv_error))
  indexcv <- rowMeans(cv_error,na.rm=TRUE)== mincv
  indexcv[is.na(indexcv)] <- FALSE
  
  return(h_diag[indexcv,])
}

## generate 10 simulations 
the_cluster_bs_CV <- makeCluster(10)
clusterSetRNGStream(the_cluster_bs_CV,27)
clusterExport(the_cluster_bs_CV,c("npnormmixEM_grid","NormKernel","rnpnormmix","x","n","h_diag","k_fold","pi_1",
                                  "pi_2","m_1","m_2","sigma_1","sigma_2","p","mu","sgm"))
h_CV <- clusterCall(cl = the_cluster_bs_CV, bs_CV_diag,1:10)
stopCluster(the_cluster_bs_CV)

## find the average of the 10 selected bandwidth
h_CV <- diag(apply(matrix(unlist(h_CV),nrow=2),1,mean))

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
  y <- rnpnormmix(p,mu,sgm)
  output <- npnormmixEM_grid(x,y,N=100,initl_m=initl_m,initl_sgm=initl_sgm,initl_pi=initl_pi,bw=h_CV) 
  RASE_pi <- sqrt((sum((output$pi1-pi_1(x))^2)+sum((output$pi2-pi_2(x))^2))/n)
  RASE_m  <- sqrt((sum((output$m1-m_1(x))^2)+sum((output$m2-m_2(x))^2))/n)
  RASE_sigma2 <- sqrt((sum((output$sigma2_1-sigma_1(x)^2)^2)+sum((output$sigma2_2-sigma_2(x)^2)^2))/n)
  return(list(RASE_pi=RASE_pi,RASE_m=RASE_m,RASE_sigma2=RASE_sigma2))
}

n_sim <- 500
the_cluster <- makeCluster(16)
clusterSetRNGStream(the_cluster,27)
clusterExport(the_cluster,c("npnormmixEM_grid","NormKernel","rnpnormmix","x","n","h_CV","initl_m","initl_sgm","initl_pi","pi_1",
                            "pi_2","m_1","m_2","sigma_1","sigma_2","p","mu","sgm"))
test_sim <- parLapply(the_cluster,1:n_sim,sim_npmixreg)
stopCluster(the_cluster)

RASEs <- matrix(unlist(test_sim),nrow=3,byrow = FALSE)
RASE_mean <- apply(RASEs,1,mean)
RASE_sd <- apply(RASEs,1,sd)

write.csv(c(h_CV,RASE_mean,RASE_sd),"ex1outputCVn400.csv")


