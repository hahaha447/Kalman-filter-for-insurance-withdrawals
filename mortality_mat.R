setwd("C:/Users/Arche/graduate/study/project/kalman/material")
library(readxl)
qxd_f <- read_xlsx("qxd_f.xlsx")$"0"
qxd_m <- read_xlsx("qxd_m.xlsx")$"0"

qxh_mns <- read_xlsx("qxh_mns.xlsx")$"0"
qxh_fns <- read_xlsx("qxh_fns.xlsx")$"0"
qxh_ms <- read_xlsx("qxh_mns.xlsx")$"0"
qxh_fs <- read_xlsx("qxh_fns.xlsx")$"0"

ix_f <- read_xlsx("ix_f.xlsx")$"0"
ix_m <- read_xlsx("ix_m.xlsx")$"0"


Generate_mortvec <- function(start_age,gender,n_sim,smoke){
temp1 <- rep(0,1000)

if (gender=="f") {
  qxd <- qxd_f
  qxh <- qxh_fns
  ix <- ix_f
  if(smoke==T){
    qxh <- qxh_fs
  }
}else{
  qxd <- qxd_m
  qxh <- qxh_mns
  ix <- ix_m
  if(smoke==T){
    qxh <- qxh_ms
  }
}

result <- matrix(0,108-start_age,n_sim)

for (i in start_age:106) {
  j <- i-start_age+1
  healthy <- which(result[j,]==0)
  unhealthy <- which(result[j,]==1)
  dead <- which(result[j,]==2)
  set.seed(1000)
  result[j+1,healthy] <- sample(c(0,1,2),prob = c((1-ix[i+1]-qxh[i+1]),ix[i+1],qxh[i+1]),length(healthy),replace=T)
  result[j+1,unhealthy] <- sample(c(1,2),prob = c((1-qxd[i+1]),qxd[i+1]),length(unhealthy),replace=T)
  result[j+1,dead] <- 2
}
result <- rbind(result,rep(2,n_sim))
return(result)
}

a <- Generate_mortvec(30,"f",50,T)
a

paste("qxd_","m",sep = "")

