
###### Parameter setting ######

dimen <- 2
theta <- c(atan(0.5)) #Slope of real threshold lines
rho <- c(2) # y-intercepts of lines
n.re <- 2*length(theta)
if ((n.re == 4) & (theta[1] == theta[2])) {n.re <- 3}
mp <- 0
d <- 1 # Lag length
phi <- list(
  phi1 = list(
    matrix(c(0,0),nrow=2),
    matrix(c(0.1,-0.6,0.5,-0.1),byrow=T,nrow=2)
  ),
  phi2 = list(
    matrix(c(1,0),nrow=2),
    matrix(c(0.4,-0.4,-0.3,0.5),byrow=T,nrow=2),
    matrix(c(-0.6,0.2,0.4,-0.6),byrow=T,nrow=2)
  )
)
sigma <- list(
  matrix(c(4,2,2,4),nrow=2),
  matrix(c(4,2,2,4),nrow=2)
)
###### Check ######

if (n.re!=length(phi)){print("Error! Length of threshold parameter and AR orders does not match!!")}
if (length(phi)!=length(sigma)){print("Error! Length of AR orders and Sigma does not match!!")}

###### Generate a Realization of the series ######
Realization(n,theta,rho,d,phi,sigma)
ao <- as.vector(sapply(phi[1:length(phi)],length)-1)