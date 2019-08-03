setwd("C:/Users/Arche/graduate/study/project/kalman")
library(ggplot2)
library(reshape)
#####################################################################################
#                             Doing nothing 
#
#####################################################################################

#generate return matrix
Generate_return <- function(duration=30,mu=0.04,sig=0.1,n_sim=1000){
  set.seed(1000)
  Z <- matrix(data = rnorm(n_sim * duration),nrow = duration,ncol = n_sim)
  rset <- exp(mu - 0.5 * sig^2) * exp(sig)^(Z)
  return(rset)
}


#generate path matrix with input return matrix
Generate_path <- function(return_mat=return1,initial_price=500){
  inner_func <- function(return_vec,S0){
    l <- length(return_vec)
    St <- rep(S0,l+1)
    for (i in 1:l) {
      St[i+1] <- prod(S0,return_vec[1:i])
    }
    return(St)
    
  }
  
  path_mat <- apply(return_mat,2,inner_func,S0=initial_price)
  return(path_mat)
}

#generate plot with a x rows matrix.
Generate_plot <- function(mat,start_y,end_y,colo="#619CFF"){
  mat <- cbind(seq(start_y,end_y),mat)
  reshaped_df <- melt(data.frame(mat),id.vars=1)
  p <- ggplot(reshaped_df,aes(x=X1,y=value,group=variable))+
    geom_line(size=0.2, alpha=0.1,col=colo)
  return(p)
}

# funciton to calculate d0 with retern
myfun <- function(r_vec,initial_asset,I,EndTime=30){
  I[1] <- 1
  x <- prod(initial_asset,r_vec[1:EndTime])
  a <- rep(0,EndTime)
  for (k in 1:EndTime) {
    a[k] <- prod(r_vec[k:EndTime],I[1:k])
  }
  y <- sum(a)
  d0 <- x/y
  return(d0)
}
#####################################################################################
#                             Kalman
#
#####################################################################################
KF <- function(mu0,sigma0,P,Q,dt=1,real_return,R_coeff=c(0,0)){
  
  F0 <- diag(2)
  x <- matrix(c(mu0,sigma0^2),2,1)
  P0 <- P
  n <- length(real_return)
  it <- log(real_return)
  result <- matrix(0,n,2)
  for (i in 1:n) {
    mean_it <- mean(it[1:i])
    var_it <- var(it[1:i])
    var_it[is.na(var_it)] <- 1000
    z <- matrix(c(mean_it,var_it),2,1)
    R <- matrix(c(var_it/i*R_coeff[1],0,0,2*var_it^2/i),2,2)
    
    H <- matrix(c(dt,0,-0.5*dt,dt),2,2)
    x_pred <- F0%*%x
    P_pred <- F0%*%P0%*%t(F0)+Q
    y <- z-H%*%x_pred
    S <- H%*%P_pred%*%t(H)+R
    K <- P_pred%*%t(H)%*%solve(S)
    x <- x_pred+K%*%y
    P0 <- (diag(2)-H%*%K)%*%P_pred
    result[i,] <- x
  }
  return(result)
  
}


NKF <- function(real_return=actural_return){
  result <- matrix(0,length(real_return),2)
  for (i in 1:length(real_return)) {
    result[i,] <-  c((mean(log(real_return[1:i])))+0.5*var(log(real_return[1:i])),var(log(real_return[1:i])))
    
  }
  return(result)

}

#####################################################################################
#                             Kalman and mean
#
#####################################################################################
return1 <- Generate_return(duration = 30)
actural_return <- return1[,50]
KF_result <- KF(mu0 = 0.08,sigma0 = 0.2,P=matrix(c(0.5,0,0,0.01),2,2),Q=matrix(c(0.01,0,0,0.0001),2,2),dt=1,real_return=actural_return,
             R_coeff = c(300,3000))
NKF_result <- NKF(actural_return)
#####################################################################################
#                             plotting
#
#####################################################################################
plot(NKF_result[,1],type = "l",col="red")
lines(KF_result[,1],col="Blue")
abline(h=0.0671)

plot(NKF_result[,2],type = "l",col="red")
lines(KF_result[,2],col="Blue")
abline(h=0.01)
#####################################################################################
#                             calculate d0
#
#####################################################################################
#
calculate_d0_NR <-  function(n_year,quantile=0.05,S0,real_return,d0_vec){
  n <- quantile*length(d0_vec)
  d0_optimal <- sort(d0_vec)[n]
  d0_result <- rep(d0_optimal,n_year)
  remaining_result <- rep(S0,n_year)
  for (i in 1:n_year) {
    S0 <- (S0-d0_optimal)*real_return[i] 
    remaining_result[i] <- S0
  }
  return(list("Withdrwal"=d0_result,"Remaining_value"=remaining_result))
}
  

#----------------------------------------------------------------------------------------------------
calculate_d0_Naive <- function(n_year,quantile=0.05,S0,real_return,mu,sig){
  d0_result <- rep(0,n_year)
  remaining_result <- rep(0,n_year)
  mu0 <- mu
  sig0 <- sig
  for (i in 1:n_year) {
    remaining_year <- n_year+1-i
    return <- Generate_return(duration = remaining_year,mu = mu0,sig = sig0)
    d0 <- apply(X=return,2,FUN = myfun,initial_asset=S0,I=rep(1,remaining_year),EndTime=remaining_year)
    n <- quantile*length(d0)
    d0_optimal <- sort(d0)[n]
    S0 <- (S0-d0_optimal)*real_return[i]
    d0_result[i] <- d0_optimal
    remaining_result[i] <- S0
  }
  return(list("Withdrwal"=d0_result,"Remaining_value"=remaining_result))
}



#----------------------------------------------------------------------------------------------------
calculate_d0_KF <- function(n_year,quantile=0.05,S0,real_return,KF_output){
  d0_result <- rep(0,n_year)
  remaining_result <- rep(0,n_year)
  for (i in 1:n_year) {
    remaining_year <- n_year+1-i
    mu0 <- KF_output[i,1]
    sig0 <- KF_output[i,2]
    return <- Generate_return(duration = remaining_year,mu = mu0,sig = sig0)
    d0 <- apply(X=return,2,FUN = myfun,initial_asset=S0,I=rep(1,remaining_year),EndTime=remaining_year)
    n <- quantile*length(d0)
    d0_optimal <- sort(d0)[n]
    S0 <- (S0-d0_optimal)*real_return[i]
    d0_result[i] <- d0_optimal
    remaining_result[i] <- S0
  }
  return(list("Withdrwal"=d0_result,"Remaining_value"=remaining_result))
}






plot(optimal_nkf$Withdrwal,type = "l",col="red")
lines(optimal_kf$Withdrwal,col="blue")


plot(optimal_nkf$Remaining_value,type = "l",col="red")
lines(optimal_kf$Remaining_value,col="blue")

#####################################################################################
#                             Reporting
#
#####################################################################################
#Method 1
return1 <- Generate_return(duration = 30,mu = 0.04,sig = 0.1,n_sim=1000)
return2 <- Generate_return(duration = 30,mu = 0.08,sig = 0.1,n_sim=1000)
return3 <- Generate_return(duration = 30,mu = 0.015,sig = 0.1,n_sim=1000)

path1 <- Generate_path(return1)
path2 <- Generate_path(return2)
path3 <- Generate_path(return3)


Generate_plot(path1,2000,2030,colo = "#619CFF")+labs(title = "mu=0.04")+theme(plot.title=element_text(hjust = 0.5))
Generate_plot(path2,2000,2030,colo = "#00BA38")+labs(title = "mu=0.08")+theme(plot.title=element_text(hjust = 0.5))
Generate_plot(path3,2000,2030,colo = "#F8766D")+labs(title = "mu=0.0015")+theme(plot.title=element_text(hjust = 0.5))



d0_1 <- apply(X=return1,2,FUN = myfun,initial_asset=500,I=rep(1,30),EndTime=30)
d0_2 <- apply(X=return2,2,FUN = myfun,initial_asset=500,I=rep(1,30),EndTime=30)
d0_3 <- apply(X=return3,2,FUN = myfun,initial_asset=500,I=rep(1,30),EndTime=30)
optimal_NR <- calculate_d0_NR(n_year = 30,quantile=0.05,S0 = 500,real_return =actural_return,d0_vec = d0)

##########  Changing ONLY d0_vec to get ruin probability
mylist <- apply(return1,2,FUN = calculate_d0_NR,n_year = 30,quantile=0.05,S0 = 500,d0_vec = d0_3)
a <- rep(0,1000)
for (i in 1:1000) {
  a[i] <- mylist[[i]]$Remaining_value[30]<=0
}
sum(a)/1000
mylist[[1]]$Withdrwal
rm(mylist)
rm(a)
#########
qplot(d0_1,fill = I("#619CFF"))+labs(title = "mu=0.04")+theme(plot.title=element_text(hjust = 0.5))
qplot(d0_2,fill = I("#00BA38"))+labs(title = "mu=0.08")+theme(plot.title=element_text(hjust = 0.5))
qplot(d0_3,fill = I("#F8766D"))+labs(title = "mu=0.0015")+theme(plot.title=element_text(hjust = 0.5))


#Method 2
actural_return <- return2[,which(order(apply(return2, 2,sum))==500)]
KF_result <- KF(mu0 = 0.05,sigma0 = 0.2,P=matrix(c(0.5,0,0,0.01),2,2),Q=matrix(c(0.01,0,0,0.0001),2,2),
                dt=1,real_return=actural_return,R_coeff = c(300,3000))
optimal_kf <- calculate_d0_KF(30,S0=500,real_return = actural_return,KF_output = KF_result)
mean(optimal_kf$Withdrwal)
qplot(x=1:30,y=KF_result[,1],geom = "line",size=I(1.5))+
  geom_hline(yintercept = 0.015,col="red",linetype="dashed",size=1.5)+
  labs(title = "updated mu")

qplot(x=1:30,y=KF_result[,2],geom = "line",size=I(1.5))+
  geom_hline(yintercept = 0.01,col="red",linetype="dashed",size=1.5)+
  labs(title = "updated sigma^2")

#Method 3
actural_return <- return2[,which(order(apply(return2, 2,sum))==500)]
optimal_naive <- calculate_d0_Naive(30,S0=500,real_return = actural_return,mu = 0.04,sig = 0.1)
mean(optimal_naive$Withdrwal)
qplot(x=1:30,y=optimal_kf$Remaining_value,geom = "line",size=I(1.5))+
  geom_line(y=optimal_naive$Remaining_value,col="red",linetype="dashed",size=1.5)+
  labs(title = "actural mu = 0.08")
