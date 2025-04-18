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
library(missForest)
library(cowplot)
library(gridGraphics)
library(plot3D)

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

#############################################################################################################################
##############################################################################################################################
##############  Generate plots of true and estimated mean functions 
#############################################################################################################################
#############################################################################################################################

x1_grid <- seq(0,1,l=100)
x2_grid <- seq(0,1,l=100)
M <- mesh(x1_grid,x2_grid)
x_grid <- cbind(c(M$x),c(M$y))
mat_true_m1 <- matrix(m_1(x_grid),100)
mat_true_m2 <- matrix(m_2(x_grid),100)


## set bandwidth matrix based on AMISE 
h_AMISE_1 <- matrix(c(0.0284887185661849,-0.0239429837615277,-0.0239429837615277,0.0265457973454218),nrow=2) ## n=200
h_AMISE_2 <- matrix(c(0.0139342177428565,-0.00856112031775077,-0.00856112031775077,0.01115489550725),nrow=2) ## n=400
h_AMISE_3 <- matrix(c(0.00799299116210106,-0.00520277766436521,-0.00520277766436521,0.00765083225087343),nrow=2) ## n=800

###################################################   Case: n = 200   #######################################################
set.seed(389)
n <- 200
x <- cbind(runif(n),runif(n))  ## generate x from uniform distributions 
# x <- rmvnorm(n,rep(0,2),autocorr.mat(2,0.5))  ## generate x from normal distribution 
p <- cbind(pi_1(x),pi_2(x)) 
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))

set.seed(579)
initl_m <- cbind(m_1(x)+rnorm(1,sd=0.5),m_2(x)+rnorm(1,sd=0.5))  ### set initial values for the mean function 
initl_sgm <- cbind(rep(mean(sigma_1(x)),n),rep(mean(sigma_2(x)),n))   ## set the initial values for the standard deviation function
initl_pi <- cbind(rep(mean(pi_1(x)),n),rep(mean(pi_2(x)),n))

y <- rnpnormmix(p,mu,sgm)
output <- npnormmixEM_grid(x,y,N=100,initl_m=initl_m,initl_sgm=initl_sgm,initl_pi=initl_pi,bw=h_AMISE_1)

inter_m1 <- interpp(x=x[,1],y=x[,2],z=output$m1,xo=c(M$x),yo=c(M$y),linear=FALSE,extrap=TRUE)$z
inter_m2 <- interpp(x=x[,1],y=x[,2],z=output$m2,xo=c(M$x),yo=c(M$y),linear=FALSE,extrap=TRUE)$z
mat_est1_m1 <- matrix(inter_m1,nrow=100)
mat_est1_m2 <- matrix(inter_m2,nrow=100)

###################################################   Case: n = 400   #######################################################
set.seed(389)
n <- 400
x <- cbind(runif(n),runif(n))  ## generate x from uniform distributions 
# x <- rmvnorm(n,rep(0,2),autocorr.mat(2,0.5))  ## generate x from normal distribution 
p <- cbind(pi_1(x),pi_2(x)) 
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))

set.seed(579)
initl_m <- cbind(m_1(x)+rnorm(1,sd=0.5),m_2(x)+rnorm(1,sd=0.5))  ### set initial values for the mean function 
initl_sgm <- cbind(rep(mean(sigma_1(x)),n),rep(mean(sigma_2(x)),n))   ## set the initial values for the standard deviation function
initl_pi <- cbind(rep(mean(pi_1(x)),n),rep(mean(pi_2(x)),n))

y <- rnpnormmix(p,mu,sgm)
output <- npnormmixEM_grid(x,y,N=100,initl_m=initl_m,initl_sgm=initl_sgm,initl_pi=initl_pi,bw=h_AMISE_2)

inter_m1 <- interpp(x=x[,1],y=x[,2],z=output$m1,xo=c(M$x),yo=c(M$y),linear=FALSE,extrap=TRUE)$z
inter_m2 <- interpp(x=x[,1],y=x[,2],z=output$m2,xo=c(M$x),yo=c(M$y),linear=FALSE,extrap=TRUE)$z
mat_est2_m1 <- matrix(inter_m1,nrow=100)
mat_est2_m2 <- matrix(inter_m2,nrow=100)

###################################################   Case: n = 800   #######################################################
#set.seed(389)
set.seed(497)
n <- 800
x <- cbind(runif(n),runif(n))  ## generate x from uniform distributions 
# x <- rmvnorm(n,rep(0,2),autocorr.mat(2,0.5))  ## generate x from normal distribution 
p <- cbind(pi_1(x),pi_2(x)) 
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))

set.seed(579)
initl_m <- cbind(m_1(x)+rnorm(1,sd=0.5),m_2(x)+rnorm(1,sd=0.5))  ### set initial values for the mean function 
initl_sgm <- cbind(rep(mean(sigma_1(x)),n),rep(mean(sigma_2(x)),n))   ## set the initial values for the standard deviation function
initl_pi <- cbind(rep(mean(pi_1(x)),n),rep(mean(pi_2(x)),n))

y <- rnpnormmix(p,mu,sgm)
output <- npnormmixEM_grid(x,y,N=100,initl_m=initl_m,initl_sgm=initl_sgm,initl_pi=initl_pi,bw=h_AMISE_3)

inter_m1 <- interpp(x=x[,1],y=x[,2],z=output$m1,xo=c(M$x),yo=c(M$y),linear=FALSE,extrap=TRUE)$z
inter_m2 <- interpp(x=x[,1],y=x[,2],z=output$m2,xo=c(M$x),yo=c(M$y),linear=FALSE,extrap=TRUE)$z
mat_est3_m1 <- matrix(inter_m1,nrow=100)
mat_est3_m2 <- matrix(inter_m2,nrow=100)

#### generate 3D plots
par(mfrow = c(4, 2))
surf3D(M$x, M$y,mat_true_m1,colvar = mat_true_m1,shade = 0.3, theta = 60, phi = 0, bty = "b2",main="True m1",xlab="x1",ylab="x2",zlab="m1(x)",zlim=c(4,8),clim=c(4,8))
surf3D(M$x, M$y,mat_est1_m1,colvar = mat_est1_m1,shade = 0.3, theta = 60, phi = 0, bty = "b2",main="Estimated m1 (n=200)",xlab="x1",ylab="x2",zlab="m1(x)",zlim=c(4,8),clim=c(4,8))
surf3D(M$x, M$y,mat_est2_m1,colvar = mat_est2_m1,shade = 0.3, theta = 60, phi = 0, bty = "b2",main="Estimated m1 (n=400)",xlab="x1",ylab="x2",zlab="m1(x)",zlim=c(4,8),clim=c(4,8))
surf3D(M$x, M$y,mat_est3_m1,colvar = mat_est3_m1,shade = 0.3, theta = 60, phi = 0, bty = "b2",main="Estimated m1 (n=800)",xlab="x1",ylab="x2",zlab="m1(x)",zlim=c(4,8),clim=c(4,8))
#par(mfrow = c(1, 1))

#par(mfrow = c(2, 2))
surf3D(M$x, M$y,mat_true_m2,colvar = mat_true_m2,shade = 0.3, theta = 60, phi = 0, bty = "b2",main="True m2",xlab="x1",ylab="x2",zlab="m2(x)",zlim=c(-1.2,1.6),clim=c(-1.2,1.6))
surf3D(M$x, M$y,mat_est1_m2,colvar = mat_est1_m2,shade = 0.3, theta = 60, phi = 0, bty = "b2",main="Estimated m2 (n=200)",xlab="x1",ylab="x2",zlab="m2(x)",zlim=c(-1.2,1.6),clim=c(-1.2,1.6))
surf3D(M$x, M$y,mat_est2_m2,colvar = mat_est2_m2,shade = 0.3, theta = 60, phi = 0, bty = "b2",main="Estimated m2 (n=400)",xlab="x1",ylab="x2",zlab="m2(x)",zlim=c(-1.2,1.6),clim=c(-1.2,1.6))
surf3D(M$x, M$y,mat_est3_m2,colvar = mat_est3_m2,shade = 0.3, theta = 60, phi = 0, bty = "b2",main="Estimated m2 (n=800)",xlab="x1",ylab="x2",zlab="m2(x)",zlim=c(-1.2,1.6),clim=c(-1.2,1.6))
par(mfrow = c(1, 1))





