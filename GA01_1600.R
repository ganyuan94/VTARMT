#install.packages("doMC")
#install.packages("foreach")

library(MASS)
library(doMC)
library(foreach)
source("fun2.R")
set.seed(12345)

###### Genetic Algorithm to Find the Optimal Sol. ######
n <- 1600
p0 <- 4 # max AR order
print(Sys.time())
n.Island <- 20 # no of islands
n.chor <- 50 # no of chromosomes in each island
step0 <- 4*p0 # least number of sample in each region
pc <- 0.9 #Prob. of crossover
pm <- 1-pc #Prob. of mutation
d <- 1
rept <- 400

resultth <- matrix(rep(NA,2*rept),nrow=2)
resultrh <- matrix(rep(NA,2*rept),nrow=2)
resultre <- matrix(rep(NA,rept),nrow=1) # no. of regions
resultmp <- matrix(rep(NA,rept),nrow=1) # Merge Pattern
resultmd <- matrix(rep(NA,rept),nrow=1)
resultphi <- NULL
resultsig <- NULL
resultao <- matrix(rep(NA,4*rept),nrow=4) # ar order in each region
threshold.th.est=list(NULL,NULL)
threshold.rh.est=list(NULL,NULL)
ar.est=list(NULL,NULL)
xx <- array(0,c(4,4,rept))
registerDoMC(100)

foreach(i = 1:rept) %dopar% {
  #for (i in 1:rept) {
  resultphi <- NULL
  resultsig <- NULL
  print(c("Repetition",i))
  source("tsgenerator_01.R")
  ans <- GA(n.Island,n.chor,pc,p0)
  n.re[i] <- length(ans$ar.E)
  resultmp[i] <- ans$mp.E # merge patterns for each reptitions, only for two threshold lines case, 11 different patterns in total
  resultmd[i] <- ans$md # store the mdl value
  resultphi <- t(ans$phi.E) #AR parameters
  resultsig <- ans$sig.E #AR parameters
  resultre <- ans$re # a vector of length n that stores the regions all observations falls in, with the estimated theta(s) and rho(s)
  truere <- sapply(1:(nrow(yp[[d]])),region,y=yp[[d]],theta,rho) # a vector of length n that stores the regions all observations falls in, with the true theta(s) and rho(s)
  
  n.re.true <- max(truere)
  
  #### Table for regime misclassification ####
  for (ii in 1:n.re.true){
    for (jj in 1:n.re[i]){
      xx[ii,jj,i] <- sum((truere==ii) & (resultre==jj))
    }
  }
  write.table(xx[,,i],file="xx01_1600.csv",append=T,quote=F,row.names=F, col.names=F, sep = ",")
  #### Threshold line result ####
  for (j in 1:length(ans$thres.th.E)){
    resultth[j,i] <- ans$thres.th.E[j]
    resultrh[j,i] <- ans$thres.rh.E[j]
    if (! is.na(resultrh[j,i])) { # adjust to make rho>=0, and 0 <= theta < 2*pi
      if (resultrh[j,i]<0) {resultrh[j,i] <- -resultrh[j,i]; resultth[j,i] <- resultth[j,i] + pi}
      while ((resultth[j,i]>=2*pi)|(resultth[j,i]<0)) {resultth[j,i] <- resultth[j,i] - 2*pi} 
    }
  }
  #### AR order in each regime ####
  for (j in 1:n.re[i]){
    resultao[j,i] <- ans$ar.E[j]
  }
  
  Est <- c(i,resultth[,i],resultrh[,i],resultmp[i],resultao[,i],resultmd[i])
  Est_ <- cbind(i,resultphi,resultsig)
  write.table(t(as.matrix(Est)),file="est01_1600.csv",append=T,quote=F,row.names=F, col.names=F, sep = ",") # Threshold result
  write.table(Est_,file="est01_1600s.csv",append = T, quote = F,row.names=F, col.names=F, sep = ",")  #Threshold result + AR Parameter
}
print(Sys.time())

