#################################################################################################################################################
##### define a function to implement the modified EM agorithm for two component nonparametric mixture of regressions based on 100 grid points
##### N is the number of grid points, initl_m, initl_sgm,initl_pi are n by c (c=2)
npnormmixEM_grid <- function(x,y,N=100,initl_m,initl_sgm,initl_pi,maxiter=2000,epsilon=0.00001,bw){
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
        inter_m1[which(is.na(inter_m1))[i_na]] <- m1_nna[dis==min(dis)][1]}}

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
        inter_m2[which(is.na(inter_m2))[i_na]] <- m2_nna[dis==min(dis)][1]}}
      
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
        inter_sgm2_1[which(is.na(inter_sgm2_1))[i_na]] <- sgm2_1_nna[dis==min(dis)][1]}}
      
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
        inter_sgm2_2[which(is.na(inter_sgm2_2))[i_na]] <- sgm2_2_nna[dis==min(dis)][1]}}
      
      inter_sgm2_2_narm <- missForest(cbind(x,inter_sgm2_2))$ximp
      inter_sgm2_2 <- inter_sgm2_2_narm[,3]
    }
    
    pi_EM[,1] <- inter_pi1
    pi_EM[,2] <- inter_pi2
    m_EM[,1] <- inter_m1
    m_EM[,2] <- inter_m2
    sgm_EM[,1] <- sqrt(inter_sgm2_1)
    sgm_EM[,2] <- sqrt(inter_sgm2_2)
    print(l)
    print(theta_diff)
  }
  if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(pi1=inter_pi1,pi2=inter_pi2,m1=inter_m1,m2=inter_m2,sigma2_1=inter_sgm2_1,sigma2_2=inter_sgm2_2,iter=l))
}

####################################################################################################
H <- 0.05*diag(2)  ## set the bandwidth matrix 
start_time <- Sys.time()
output1 <- npnormmixEM_grid(x,y,initl_m=initl_m,initl_sgm=initl_sgm,initl_pi=initl_pi,bw=H_opt) 
end_time <- Sys.time()

end_time - start_time
output1 <- npnormmixEM_grid(x,y,initl_m=initl_m,initl_sgm=initl_sgm,initl_pi=initl_pi,bw=H) 


dist_2 <- (x[!is.na(inter_pi1),] - x[is.na(inter_pi1),])^2
z_nona <- inter_pi1[!is.na(inter_pi1)]


library(missForest)
help(missForest)
inter_pi1_narm <- missForest(cbind(x,inter_pi1))$ximp

if (anyNA(inter_pi1)){
  inter_pi1_narm <- missForest(cbind(x,inter_pi1))$ximp
  inter_pi1 <- inter_pi1_narm[,3]
}



