#install.packages("doMC")
#install.packages("foreach")

library(MASS)
#library(doMC)
#library(foreach)

Reord <- function(th,rh){
  if (th[1] != th[2]) {
    sth <- sort(th)
    srh <- rh[order(th)]
  } else {
    sth <- th
    srh <- sort(rh)
  }
  return(list(th = sth, rh = srh))
}

#######################################
#    Initialize space for parameters  #
#######################################

Init <- function(n.chor,n.Island){
  temp <- NULL
  
  th <<- NULL # slopes
  rh <<- NULL # intercepts
  mp <<- NULL # merge pattern
  ao <<- NULL # AR order
  re <<- NULL # region numbering
  md <<- NULL # discription length
  
  for (j in 1:n.chor) {temp[[j]] <- rep(0, 2)}
  for (i in 1:n.Island) {
    th[[i]] <<- temp # [[n.Island]] [[n.chor]] [k] (k = 1 or 2)
    rh[[i]] <<- temp # [[n.Island]] [[n.chor]] [k] (k = 1 or 2)
    ao[[i]] <<- temp # [[n.Island]] [[n.chor]] [k] (k = 1, 2, ..., n-p_0)
    re[[i]] <<- temp # [[n.Island]] [[n.chor]] [k] (k = 1, 2, ..., n-p_0)
    mp[[i]] <<- rep(0, n.chor) # [[n.Island]] [k] (k = 1, 2, ..., n.chor)
    md[[i]] <<- rep(0, n.chor) # [[n.Island]] [k] (k = 1, 2, ..., n.chor)
  }
  smd <<- md # To store sorted discription length
  
  elite.th <<- NULL 
  elite.rh <<- NULL
  elite.md <<- NULL
  elite.mp <<- NULL
  elite.ao <<- NULL
  elite.re <<- NULL
  
  for (i in 1:n.Island) {
    elite.th[[i]] <<- rep(0, 2)
    elite.rh[[i]] <<- rep(0, 2)
    elite.ao[[i]] <<- rep(0, 2)
    elite.re[[i]] <<- rep(0, 2)
    elite.md[i] <<- NA
    elite.mp[i] <<- NA
  }
}

##################################
#    Generate New Chromosome     #
##################################

TAR.chromosome <- function(){
  a <- sample(0:2,1) # No. of lines
  #a <- 2
  if (a!=0){
    p <- 0
    while (p == 0) {
      th <- runif(a,0,2*pi)
      rh <- runif(length(th),0,max(z))
      n.r <- 2*a
      re <- sapply(1:(nrow(z)),region,y=z,th,rh)
      if (a==2) {
        temp <- Reord(th,rh)
        th <- temp$th
        rh <- temp$rh
        mp <- sample(c(0,5,6),1)
        temp <- mer(re,mp,n.r)
        n.r <- temp$n.r
        re <- temp$re
      } else  {mp <- NA}
      ao <- sample(p0,n.r,TRUE)
      p <- 1
      for (i in 1:n.r) {
        p <- p * (length(re[re==i])>=step0)
      }
    }
  } else {
    th <- NA # infinity th & rh means no threshold
    rh <- NA
    ao <- sample(p0,1)
    re <- rep.int(1,n-p0)
    mp <- NA
  }
  return(list(th=th,rh=rh,re=re,ao=ao,mp=mp))
}

###################################
#     Store chromosome config     # 
###################################

chromosome <- function(l){
  
  i<-l%%n.Island+1
  j<-l%/%n.Island+1
  
  temp.tar=TAR.chromosome()
  
  th[[i]][[j]]<<-temp.tar$th
  rh[[i]][[j]]<<-temp.tar$rh
  ao[[i]][[j]]<<-temp.tar$ao
  re[[i]][[j]]<<-temp.tar$re    # re[[i]][[j]][t] stores the region of yy[t-d]
  mp[[i]][j] <<-temp.tar$mp
}

###############################################
#    Return the number of region of an obs.   # 
###############################################

region <- function(y,t,theta,rho) { 
  if (is.na(theta[1])) {return(1)} else {
    ab <- (cos(theta) * y[t,1] + sin(theta) * y[t,2]  >=  rho)
    if ((length(theta)==2)){
      if (abs(theta[1] - theta[2]) == pi) return(sum(ab*(1:2))+1) else {
        if (theta[1]==theta[2]) return(sum(ab)+1) else {
          if (ab[1]) return(sum(ab)) else return (length(ab)*2-sum(ab)) }}
    } else {
      if (ab[1]) return(sum(ab)) else return (length(ab)*2-sum(ab))
    }
  }
}

#########################
#    Merge subregions   # 
#########################

mer <- function(re,mp,n.r){
  
  switch(mp,
         {re[re>=2] <- re[re>=2]-1;n.r<-3}, #CASE1: merge 1,2
         {re[re>=3] <- re[re>=3]-1;n.r<-3}, #CASE2: merge 2,3
         {re[re==4] <- 3;n.r<-3}, #CASE3: merge 3,4
         {re[re==4] <- 1;n.r<-3}, #CASE4: merge 4,1
         {re[re==1] <- 3; re <- re-1;n.r<-3}, #CASE5: merge 1,3
         {re[re==4] <- 2;n.r<-3}, #CASE6: merge 2,4
         {re[re>=3] <- re[re>=3]-2; n.r<-2}, #CASE7: merge 1,3 and 2,4
         {re[re<=3] <- 1; re[re==4] <- 2; n.r<-2}, #CASE8: merge 1,2,3
         {re[re>=2] <- 2; n.r <-2}, #CASE9: merge 2,3,4
         {re[re>=3] <- 1;n.r<-2}, #CASE10: merge 3,4,1
         {re[re==2] <- 1; re[re==4] <- 1; re[re==3] <- 2;  n.r<-2} #CASE11: merge 4,1,2
  )
  return(list(re = re,n.r = n.r))
}

AIC <- function(l) {
  i <- l %% n.Island + 1
  j <- l %/% n.Island + 1
  
  cth <- th[[i]][[j]]
  cre <- re[[i]][[j]]
  cao <- ao[[i]][[j]]
  
  AIC <- 0
  for (rr in 1:max(cre)){
    AIC <- AIC + 2*2*(2*cao[rr]+1)
    index <- (1:length(cre))*(cre==rr)
    MatY <- yy[index,]
    n.j <- nrow(MatY)
    MatX <- rep.int(1,n.j)
    for (ii in 1:cao[rr]) {
      MatX <- cbind(MatX,yp[[ii]][index,])
    }
    phi_head <- solve(t(MatX) %*% as.matrix(MatX), t(MatX) %*% as.matrix(MatY))
    A_head <-  MatY - as.matrix(MatX) %*% as.matrix(phi_head)
    sigma_a <- t(A_head) %*% as.matrix(A_head) / n.j
    AIC <- AIC + n.j * (log(det(sigma_a))+log(4*pi*pi) + 2)
  }
  md[[i]][j] <<- AIC
}
################################
#    MDL for l-th chromosome   # 
################################

MDL <- function (l) {
  
  i <- l %% n.Island + 1
  j <- l %/% n.Island + 1
  
  cth <- th[[i]][[j]]
  cre <- re[[i]][[j]]
  cao <- ao[[i]][[j]]
  
  if (! is.na(cth[1])) {n.l=length(cth)} else {n.l <- 0}
  n.r <- max(cre)
  n.ll <- n.l
  if (n.ll==0) n.ll <- 1
  temp.mdl <- n.r * log2(p0) + n.l * (4+2*30)  # line number + two pts/line
  if (length(cth) == 2) {temp.mdl <- temp.mdl + log2(11)}
  for (rr in 1:n.r){
    index <- (1:length(cre))*(cre==rr)
    MatY <- yy[index,]
    n.j <- nrow(MatY)
    MatX <- rep.int(1,n.j)
    for (ii in 1:cao[rr]) {
      MatX <- cbind(MatX,yp[[ii]][index,])
    }
    phi_head <- solve(t(MatX) %*% as.matrix(MatX), t(MatX) %*% as.matrix(MatY))
    A_head <-  MatY - as.matrix(MatX) %*% as.matrix(phi_head)
    sigma_a <- t(A_head) %*% as.matrix(A_head) / n.j
    temp.mdl <- temp.mdl + (2*cao[rr]+2.5)*log(n.j) + 0.5*n.j*(log(4*pi*pi*det(sigma_a))+2)
  }
  md[[i]][j] <<- temp.mdl
}


#############################
#     Sort md into smd      # 
#############################

SIG<-function(i) {
  smd[[i]]<<-sort(md[[i]])
}

###########################################
#     CM Parents Selection/Generation     #
###########################################

CM <- function(l){
  
  i<-l %% n.Island + 1
  j<-l %/% n.Island + 1
  
  #Define the Choosing of Chromosome out function
  rnumber <- 1/(1:n.chor) #prob for parent choosing
  pnumber <- sample(n.chor,1,prob=rnumber)
  index <- (md[[i]]==smd[[i]][pnumber])*(1:n.chor)
  
  C1 <- NULL # The first parent
  C1$th <- th[[i]][index][[1]]
  C1$rh <- rh[[i]][index][[1]]
  C1$ao <- ao[[i]][index][[1]]
  C1$re <- re[[i]][index][[1]]
  C1$mp <- mp[[i]][index][1]
  
  p <- 0
  while (p == 0){
    check1 <- runif(1, 0, 1)
    if (check1 < pc) {
      pnumber <-
        sample(n.chor, 1, prob = rnumber)			# rnumber is 1/n.chor for decreasing prob
      index <- (md[[i]] == smd[[i]][pnumber]) * (1:n.chor)
      C2 <- NULL # The second parent
      C2$th <- th[[i]][index][[1]]
      C2$rh <- rh[[i]][index][[1]]
      C2$ao <- ao[[i]][index][[1]]
      C2$re <- re[[i]][index][[1]]
      C2$mp <- mp[[i]][index][1]
    } else {
      C2 <- TAR.chromosome()
    }
    p <- 1
    temp.TAR <- CO.TAR(C1,C2)
    n.rr <- temp.TAR$n.r
    for (ii in 1:n.rr) {p <- p * (length(temp.TAR$re[temp.TAR$re==ii])>=step0)}
  }
  th[[i]][[j]]<<-temp.TAR$th
  rh[[i]][[j]]<<-temp.TAR$rh
  ao[[i]][[j]]<<-temp.TAR$ao
  re[[i]][[j]]<<-temp.TAR$re
  mp[[i]][j]<<-temp.TAR$mp
}

###############################
#     Crossover Process       #
###############################

CO.TAR <- function(C1,C2){ # AR order part can be finetuned
  
  th1 <- C1$th
  rh1 <- C1$rh
  th2 <- C2$th
  rh2 <- C2$rh
  ao1 <- C1$ao[C1$re] # stores the ar.order of region where each obs. falls in 
  ao2 <- C2$ao[C2$re]
  mp1 <- C1$mp
  mp2 <- C2$mp
  
  if (is.na(th1[1])) {nl.1 <- 0} else {nl.1 <- length(th1)}
  if (is.na(th2[1])) {nl.2 <- 0} else {nl.2 <- length(th2)}
  
  # get new threshold lines
  NL <- sample(c(nl.1,nl.2),1) # number of line 
  index <- sample(1:(nl.1+nl.2),NL) # randomly choose NL line
  if (NL != 0){
    if (is.na(th1[1])) {th1 <- NULL; rh1 <- NULL}
    TH <- c(th1,th2)[index]
    RH <- c(rh1,rh2)[index]
  } else {
    TH <- NA
    RH <- NA
  }
  # get new merge pattern MP
  if (NL==2) {
    temp <- Reord(TH,RH)
    TH <- temp$th
    RH <- temp$rh
    if ((!is.na(mp1))) {
      if (!is.na(mp2)) MP <- sample(c(mp1, mp2), 1) else MP <- mp1
    }
    else MP <- mp2
  } else {MP <- NA}
  
  RE <- sapply(1:(nrow(z)),region,y=z,TH,RH)
  NR <- max(RE)
  temp <- mer(RE,MP,NR)
  RE <- temp$re
  NR <- temp$n.r
  # get new AR order 
  count <<- matrix (rep(0, p0 * NR), nrow = NR)
  for (i in 1:NR) {
    sapply(c(ao1[RE==i],ao2[RE==i]), function(t) {count[i,t] <<- count[i,t] + 1})  
  }
  if (sum(rowSums(count) == 0) == 0){
    AO <- sapply(1:NR, function(t) {sample(1:p0, size = 1, prob = count[t,])} )
  } else {AO <- NULL}
  return(list(th = TH, rh = RH, mp = MP, ao = AO, re= RE, n.r = NR))
}

###############################################
#    Save elite chromosome for all islands    #
###############################################

elitesaving <- function(){
  
  for (i in 1:n.Island) {
    elite.md[i] <<- smd[[i]][1]
    elite.mp[i] <<- mp[[i]][ md[[i]] == elite.md[i] ][1]
    elite.th[[i]] <<- th[[i]][ md[[i]] == elite.md[i] ][[1]]
    elite.rh[[i]] <<- rh[[i]][ md[[i]] == elite.md[i] ][[1]]
    elite.ao[[i]] <<- ao[[i]][ md[[i]] == elite.md[i] ][[1]]
    elite.re[[i]] <<- re[[i]][ md[[i]] == elite.md[i] ][[1]]
  }
}

#########################
#     Elitest Step      #
#########################

ELI<-function(l) {
  
  th[[l]][md[[l]]==smd[[l]][n.chor]][[1]] <<- elite.th[[l]]
  rh[[l]][md[[l]]==smd[[l]][n.chor]][[1]] <<- elite.rh[[l]]
  ao[[l]][md[[l]]==smd[[l]][n.chor]][[1]] <<- elite.ao[[l]]
  re[[l]][md[[l]]==smd[[l]][n.chor]][[1]] <<- elite.re[[l]]
  mp[[l]][md[[l]]==smd[[l]][n.chor]][1] <<- elite.mp[[l]]
  md[[l]][md[[l]]==smd[[l]][n.chor]][1] <<- elite.md[[l]]
  
  m.range=smd[[l]][round(n.chor/3)]-smd[[l]][1]
  for (i in 1:floor(n.chor/3)) {
    temp=(md[[l]]==smd[[l]][n.chor-i])*(1:length(md[[l]]))
    temp0=temp[temp>0]
    temp=temp0[sample(length(temp0),1)]
    ao[[l]][[temp]]<<-elite.ao[[l]]
    mp[[l]][[temp]]<<-elite.mp[[l]]
    re[[l]][[temp]]<<-elite.re[[l]]
    if (is.na(elite.th[[l]][1])) {
      th[[l]][[temp]]<<-elite.th[[l]]
      rh[[l]][[temp]]<<-elite.rh[[l]]
    } else {
      th[[l]][[temp]]<<-elite.th[[l]] + runif(length(elite.th[[l]]),-0.1,0.1)
      rh[[l]][[temp]]<<-elite.rh[[l]] + runif(length(elite.rh[[l]]),-0.1,0.1)
    }
    md[[l]][[temp]] <<- elite.md[[l]] + runif(1) * m.range
  }
}

#############################
#     Migration Function    #
#############################

#1st island good-> 2nd island bad ,..., n.Islandth island good -> 1st island bad

MIG<-function(w) {
  if (w>1) {
    th[[w]][md[[w]]==smd[[w]][n.chor]][[1]]<<-th[[w-1]][md[[w-1]]==smd[[w-1]][1]][[1]]
    th[[w]][md[[w]]==smd[[w]][n.chor-1]][[1]]<<-th[[w-1]][md[[w-1]]==smd[[w-1]][2]][[1]]
    rh[[w]][md[[w]]==smd[[w]][n.chor]][[1]]<<-rh[[w-1]][md[[w-1]]==smd[[w-1]][1]][[1]]
    rh[[w]][md[[w]]==smd[[w]][n.chor-1]][[1]]<<-rh[[w-1]][md[[w-1]]==smd[[w-1]][2]][[1]]
    ao[[w]][md[[w]]==smd[[w]][n.chor]][[1]]<<-ao[[w-1]][md[[w-1]]==smd[[w-1]][1]][[1]]
    ao[[w]][md[[w]]==smd[[w]][n.chor-1]][[1]]<<-ao[[w-1]][md[[w-1]]==smd[[w-1]][2]][[1]]
    re[[w]][md[[w]]==smd[[w]][n.chor]][[1]]<<-re[[w-1]][md[[w-1]]==smd[[w-1]][1]][[1]]
    re[[w]][md[[w]]==smd[[w]][n.chor-1]][[1]]<<-re[[w-1]][md[[w-1]]==smd[[w-1]][2]][[1]]
    mp[[w]][md[[w]]==smd[[w]][n.chor]][1]<<-mp[[w-1]][md[[w-1]]==smd[[w-1]][1]][[1]]
    mp[[w]][md[[w]]==smd[[w]][n.chor-1]][1]<<-mp[[w-1]][md[[w-1]]==smd[[w-1]][2]][[1]]
  } else {
    th[[w]][md[[w]]==smd[[w]][n.chor]][[1]]<<-th[[n.Island]][md[[n.Island]]==smd[[n.Island]][1]][[1]]
    th[[w]][md[[w]]==smd[[w]][n.chor-1]][[1]]<<-th[[n.Island]][md[[n.Island]]==smd[[n.Island]][2]][[1]]
    rh[[w]][md[[w]]==smd[[w]][n.chor]][[1]]<<-rh[[n.Island]][md[[n.Island]]==smd[[n.Island]][1]][[1]]
    rh[[w]][md[[w]]==smd[[w]][n.chor-1]][[1]]<<-rh[[n.Island]][md[[n.Island]]==smd[[n.Island]][2]][[1]]
    ao[[w]][md[[w]]==smd[[w]][n.chor]][[1]]<<-ao[[n.Island]][md[[n.Island]]==smd[[n.Island]][1]][[1]]
    ao[[w]][md[[w]]==smd[[w]][n.chor-1]][[1]]<<-ao[[n.Island]][md[[n.Island]]==smd[[n.Island]][2]][[1]]
    re[[w]][md[[w]]==smd[[w]][n.chor]][[1]]<<-re[[n.Island]][md[[n.Island]]==smd[[n.Island]][1]][[1]]
    re[[w]][md[[w]]==smd[[w]][n.chor-1]][[1]]<<-re[[n.Island]][md[[n.Island]]==smd[[n.Island]][2]][[1]]
    mp[[w]][md[[w]]==smd[[w]][n.chor]][1]<<-mp[[n.Island]][md[[n.Island]]==smd[[n.Island]][1]][[1]]
    mp[[w]][md[[w]]==smd[[w]][n.chor-1]][1]<<-mp[[n.Island]][md[[n.Island]]==smd[[n.Island]][2]][[1]]
  }
}

####################
#    Main Flow     #
####################

GA <- function(n.Island,n.chor,pc,p0) {
  
  # Initialzie storage space
  Init(n.chor,n.Island)
  
  # Generate the population
  print("Generate Chromosome")
  sapply(0:(n.Island*n.chor-1),chromosome)
  print(c("Find MDL",date()))
  sapply(0:(n.Island*n.chor-1),MDL)
  print(c("Finish MDL",date()))
  sapply(1:n.Island,SIG)
  elitesaving() 
  
  # Main process
  for (b in 1:40) {
    print(c("Generation",b,"start",date()))
    
    sapply(0:(n.Island*n.chor-1),MDL)
    # print(c("Generation",b,"MDL.done",date()))
    sapply(1:n.Island,SIG)
    # Doing Elite Step for Monotonic increasing of Chromosome
    sapply(1:n.Island,ELI)
    # Do SIG function again as the sigma value changes due to the elite step
    sapply(1:n.Island,SIG)
    # Saving the Elite Chromosome
    elitesaving()
    print(c("best.th",elite.th[elite.md==min(elite.md)][[1]]))
    print(c("best.rh",elite.rh[elite.md==min(elite.md)][[1]]))
    print(c("best.mp",elite.mp[elite.md==min(elite.md)][1]))
    print(c("best.md",min(elite.md)[1]))
    print(c("best.ao",elite.ao[elite.md==min(elite.md)][[1]]))
    sapply(0:(n.Island*n.chor-1),CM)
    if (b%%4==0){
      print("Do migration every four generations")
      sapply(0:(n.Island*n.chor-1),MDL)
      sapply(1:n.Island,SIG)
      sapply(1:n.Island,MIG)
    }
    
    # Close the repeating b loops
  }
  
  # Do one more calculation for finding threshold
  sapply(0:(n.Island*n.chor-1),MDL)
  sapply(1:n.Island,SIG)
  elitesaving()
  
  thres.th.est<-elite.th[elite.md==min(elite.md)][[1]]
  thres.rh.est<-elite.rh[elite.md==min(elite.md)][[1]]
  thres.ao<-elite.ao[elite.md==min(elite.md)][[1]]
  thres.mp<-elite.mp[elite.md==min(elite.md)][[1]]
  thres.md<-min(elite.md)[1]
  thres.re<-sapply(1:(nrow(z)),region,y=z,thres.th.est,thres.rh.est)
  # One more fine-tuning step
  #  for (i in seq(-0.02,0.02,0.0002)) {
  #	  for (j in seq(-0.02,0.02,0.0002)){
  #		  temp.re <- sapply(1:(nrow(yp[[d]])),region,y=yp[[d]],elite.th[elite.md==min(elite.md)][[1]]+i,elite.rh[elite.md==min(elite.md)][[1]]+j)
  #		  if (temp.re!=thres.re) {
  #  		  temp.md <- MDLCal(temp.re,thres.ao)
  #  		  if (temp.md<thres.md) {
  #          thres.md <- temp.md
  #          thres.re <- temp.re
  #          thres.th.est <- elite.th[elite.md==min(elite.md)][[1]]+i
  #          thres.rh.est <- elite.rh[elite.md==min(elite.md)][[1]]+j
  #  		  }#
  #		  }
  #	  }
  #  }
  NR <- max(thres.re)
  temp <- mer(thres.re,thres.mp,NR)
  thres.re <- temp$re
  temp <- ModelEst(thres.re,thres.ao)
  
  
  print(c("result, th=",thres.th.est,"rh =",thres.rh.est))
  print(c("mp",thres.mp,"MDL",thres.md))
  print(c("ar.order",thres.ao))
  print(c("phi",temp$phi_head))
  return(list(re = thres.re, md = thres.md,thres.th.E=thres.th.est, thres.rh.E = thres.rh.est,ar.E=thres.ao, mp.E = thres.mp, phi.E = temp$phi_head, sig.E = temp$sigma_head))
}


  ###### Genetic Algorithm to Find the Optimal Sol. ######
  p0 <- 7 # max AR order
  print(Sys.time())
  n.Island <- 30 # no of islands
  n.chor <- 30 # no of chromosomes in each island
  step0 <- 3*p0 # least number of sample in each region
  pc <- 0.9 #Prob. of crossover
  pm <- 1-pc #Prob. of mutation
  d <- 1
  
  resultth <- NULL
  resultrh <- NULL
  resultre <- NULL
  resultmp <- NULL
  resultmd <- NULL
  resultphi <- NULL
  resultsig <- NULL
  resultao <- NULL
  threshold.th.est=list(NULL,NULL)
  threshold.rh.est=list(NULL,NULL)
  ar.est=list(NULL,NULL)
  xx <- matrix(0,4,4)
  
  #READ-IN DATA
  Data <- read.csv("IR.csv",header=T)
  Data <- Data[1:564,]
  n <- length(Data[,1])
  par(mfrow=c(1,1))
  plot(Data[,1:2],type="p")
  abline(h=0,col="blue")
  abline(v=0,col="blue")
  par(mfrow=c(2,1))
  plot(Data[,1],type="l",xlab = "t", ylab = "y_{1,t}")
  plot(Data[,2],type="l",xlab = "t", ylab = "y_{2,t}")
  z <- Data[-c(1:(p0-d),(n+1-d):n),3:4]
  y <- Data[,1:2]
  #LAG/LEAD SERIES
  yy <- y[-(1:p0),]
  yp <- NULL  #yp is initiated to store lead series
  for (i in 1:(p0-1)){
    yp[[i]] <- y[-c(1:(p0-i),(n+1-i):n),]
  }
  yp[[p0]] <- y[-((n-p0+1):n),]

  ans <- GA(n.Island,n.chor,pc,p0)
  n.re <- length(ans$ar.E)
  resultmp <- ans$mp.E # merge patterns for each reptitions, only for two threshold lines case, 11 different patterns in total
  resultmd <- ans$md # store the mdl value
  resultphi <- t(ans$phi.E) #AR parameters
  resultsig <- ans$sig.E #AR parameters
  resultre <- ans$re # a vector of length n that stores the regions all observations falls in, with the estimated theta(s) and rho(s)
  #truere <- sapply(1:(nrow(z)),region,y=z,theta,rho) # a vector of length n that stores the regions all observations falls in, with the true theta(s) and rho(s)
  
  #n.re.true <- max(truere)
  
  #### Table for regime misclassification ####
  #for (ii in 1:n.re.true){
  #  for (jj in 1:n.re[i]){
  #    xx[ii,jj] <- sum((truere==ii) & (resultre==jj))
  #  }
  #}
 
  #### Threshold line result ####
  for (j in 1:length(ans$thres.th.E)){
    resultth[j] <- ans$thres.th.E[j]
    resultrh[j] <- ans$thres.rh.E[j]
    if (! is.na(resultrh[j])) { # adjust to make rho>=0, and 0 <= theta < 2*pi
      if (resultrh[j]<0) {resultrh[j] <- -resultrh[j]; resultth[j] <- resultth[j] + pi}
      while ((resultth[j]>=2*pi)|(resultth[j]<0)) {resultth[j] <- resultth[j] - 2*pi} 
    }
  }
  #### AR order in each regime ####


print(Sys.time())
