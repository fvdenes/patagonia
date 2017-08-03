###############################################################################
#### Distance sampling simulation: non-perpendicular line transect dataset ####
### This code is based on code provided in Kéry and Royle (2015), Chapter 8 ###
###############################################################################


#######################################################################################################################################################
#### A brief description of the rationale for obtaining random values for the non-perpendicular -sighting- distances of objects in a line transect ####
#######################################################################################################################################################
library(plotrix)
### A line transect of 500m is to be surveyed. There are 10 objects that can potentially be detected, and the maximum detection distance is 100m:
set.seed(1234)
N=10
z1 <- runif(N,-100,100)
z2 <- runif(N,0,500)

par(mfrow=c(1,3))

plot(z1,z2,asp=1,pch=1,main="Standard distance measurements",xlab="Transect width",ylab="Transect length") 
abline(v=0,col=2)
abline(v=-100,col=2, lty=2)
abline(v=100,col=2, lty=2)
segments(x0=z1,x1=0,y0=z2,y1=z2)
# This is standard distance sampling, with distance measurements taken at right angles in relation to the transect line


plot(z1,z2,asp=1,pch=1,main="Object detection radius",xlab="Transect width",ylab="Transect length")
abline(v=0,col=2)
draw.circle(z1[1],z2[1],radius=100, lty=2) 
# Assuming that maximum detection distance is 100m, the process of detection can also be viewed in relation to each object, where each has a "maximum detection radius" of 100m.


plot(z1,z2,asp=1,pch=1,main="Transect detection arc",xlab="Transect width",ylab="Transect length")
abline(v=0,col=2, lty=2)
draw.arc(z1,z2,radius=100,lty=3,deg1=180,deg=360)
# Considering that the observer traverses the transect from one end to the other, and that his attention is focused on what lies ahead (and on each side, but not behind!), I define a detection arc for each object. 

i <- which(z2>100)
z3 <- sqrt((100^2)-(z1[i]^2))
segments(x0=z1,x1=0,y0=z2,y1=z2)
segments(x0=z1[i],x1=0,y0=z2[i],y1=z2[i]-z3)
segments(x0=0,x1=0,y0=z2[i],y1=z2[i]-z3,col=3)
# The intersection between the arcs and the transect line mark the starting point of a segment "S" of the transect from which the object can potentially be detected by the observer (shown in green). The segment ends at the point where the object is perpendicular to the transect (farther than that the object is behind the observer).

z <- which(z2<100)
segments(x0=z1[z],x1=0,y0=z2[z],y1=0)
segments(x0=0,x1=0,y0=z2[z],y1=0,col=3)
# For objects located in the begining of the transect (before the first 100), detection is possible from the start of the transect (the Y segment is shorter).

# Each object can be detected from the begining of the respective S segment to the end of it. The probability of detection along the segment increases as the observer gets closer to the end, which is the point where the distance to the object is shortest. To simulate random sighting distance data for each detected object, I divide each segment by 100. On each 100th, I calculate the distance to the object (the hypotenuse of the triangle). Based on a pre-determined half-normal detection function (i.e. with known scale parameter sigma), I calculate the probability of detection for each 100th, and sample one out the the 100 values following their respective probabilities. 

###################################################################################################################################################################################
#### Create dataset, subject the data to sampling, perform likelihood analysis (line and point transect, conditional and full likelihood models) and obtain abundance estimate ####
###################################################################################################################################################################################

N=200 # Population size
sigma=30 # Half-normal shape parameter
B = 100 # Maximum detection distance for objects
g <- function(x, sig) exp(-x^2/(2*sig^2)) # Half-normal detection function

###############################################################################
### Function to generate dataset and obtain samples with the two approaches ###
###############################################################################

sim.ldata <- function(N=200, sigma=30, B=100){
  
  u1 <- runif(N,-B,B) # Distances to transect line for each object
  u2 <- runif(N,0,500) # Perpendicular position along transect line of each object 
  
  g <- function(x, sig) exp(-x^2/(2*sig^2)) # The half-normal function
  
  ### Regular distance sampling protocol - perpendicular distances only ###
  p <- g(u1, sig=sigma) # Detection probabilities for sightings at right angles
  y <- rbinom(N,1,p)    # Some inds. are detected and their distance measured
  x <- abs(u1[y==1])    # Perpendicular distances of detected individuals
  
  
  # For objects located u2>100, detection is possible once the observer is within 100 radius (i.e. maximum detection distance = B = 100m)
  i <- which(u2>B)
  # For objects located u2<100, maximum detection distance is slightly smaller than 100m, because the intersection of the detection arc to the transect line is before the point where the survey starts. 
  z <- which(u2<B)
  u3 <- sqrt((u1[z]^2)+(u2[z]^2))  # The maximum detection distance for objects within the first 100m of the transect, for which detection is possible from the begining of the transect
  
  sight.x <- rep(NA,length(x)) # a vector to hold the "sighting distances"
  
  for(k in 1:length(x)){
    if(k%in%i){
      vector <- seq(x[k],B,length.out = 100) # maximum detection distance is B
      sight.x[k] <- sample(x=vector,size=1,prob=g(vector,sig=sigma))
    }  
    if(k%in%z){
      vector <- seq(x[k],u3[which(z==k)],length.out = 100) # maximum detection distance is slightly smaller than B (because the intersection of the detection arc to the transect line is before the point where the survey starts)
      sight.x[k] <- sample(x=vector,size=1,prob=g(vector,sig=sigma))
    }
  }

  
  return(list(x=x,sight.x=sight.x))
}

#######################################################
### Simulation loop: analysis of simulated datasets ###
#######################################################

#set.seed(2015)

simrep <- 1000
simout <- matrix(NA, nrow=simrep, ncol=6)
colnames(simout) <- c("Nhat.standard.lt","sigma.std.lt","Nhat.sight.lt.cond","sigma.sight.lt","Nhat.sight.pt.cond","sigma.sight.pt")
N=200
B=100
sigma=30

for (sim in 1:simrep){
  
  result <- sim.ldata(N=200, sigma = 30)  
  
  ## Line transect conditional likelihood (Likelihood formulation from Kéry and Royle (2015), Chapter 8,  p. 402)##
  Lcond1 <- function(lsigma,x){  # Define conditional nll
    sigma <- exp(lsigma)
    -1*sum(log(g(x, sig=sigma)/integrate(g, 0,B, sig=sigma)$value/B)) 
  }
  
  g <- function(x, sig) exp(-x^2/(2*sig^2)) # The half-normal function
  
  # Standard protocol (perpendicular) data
  result1 <- optim(log(30), Lcond1, x=result$x , hessian=TRUE, method="Brent", lower=-5, upper=10)
  pbar1 <- integrate(g, 0,B, sig=exp(result1$par))$value/B
  n1 <- length(result$x)
  (Nhat.standard <- n1/pbar1)
  
  simout[sim,1] <- Nhat.standard
  simout[sim,2] <- exp(result1$par)
  
  # Sight distance data
  result2 <- optim(log(30), Lcond1,x=result$sight.x ,  hessian=TRUE, method="Brent", lower=-5, upper=20)
  pbar2 <- integrate(g, 0,B, sig=exp(result2$par))$value/B
  n2 <- length(result$sight.x)
  (Nhat.sight <- n2/pbar2)
  
  simout[sim,3] <- Nhat.sight
  simout[sim,4] <- exp(result2$par)
  
  ## Point transect conditional likelihood (Likelihood formulation from Kéry and Royle (2015), Chapter 8,  p. 414) .##
  # Estimate sigma using the point transect likelihood
  
  Lik.cond.point <- function (parm, data, B){ 
    sigma <- exp(parm)
    p <- exp(-data*data/(2*sigma*sigma))
    f <- 2*data/(B^2)
    pbar <- integrate(function(x, s= sigma) exp(-x^2/(2*s^2))*x, 0, B)$value*2/(B^2)
    negLL <- -1*sum(log(p*f/pbar))
    return(negLL)
  }
  
  # Sight distance data
  result3 <- optim(c(0), Lik.cond.point, data=result$sight.x, B=B, method="Brent", hessian=TRUE, lower=-10, upper=10)
  sigma.hat <- exp(result3$par) # Estimated sigma
  
  # Estimated average detection probability and conditional estimator of N
  pbar3 <- integrate(g, 0,B, sig=sigma.hat)$value/B
  
  Nhat.p.sight <- length(result$sight.x)/pbar3
  
  simout[sim,5] <- Nhat.p.sight
  simout[sim,6] <- sigma.hat
}

simout

####################################################################################
### Plot estimated abundance values from each method and compare with real value ###
####################################################################################

par(mfrow=c(2,3))
boxplot(simout[,2],ylim=c(0,80),main="Standard Line Transect - cond. lik.", xlab="Sigma")
abline(h=sigma,col=2)
boxplot(simout[,4],ylim=c(0,80),main="Sight Dist. Line Transect - cond. lik.", xlab="Sigma")
abline(h=sigma,col=2)
boxplot(simout[,6],ylim=c(0,80),main="Sight Dist. Point Transect - cond. lik.", xlab="Sigma")
abline(h=sigma,col=2)

hist(simout[,1], main="Standard Line Transect - cond. lik.", xlab="Nhat", xlim=c(0,500))
abline(v=N, col=2)
hist(simout[,3], main="Sight Dist. Line Transect - cond. lik.", xlab="Nhat", xlim=c(0,500))
abline(v=N, col=2)
hist(simout[,5], main="Sight Dist. Point Transect - cond. lik.", xlab="Nhat", xlim=c(0,500))
abline(v=N, col=2)

################################
### Discussion and questions ###
################################
# 1: As expected, standard line transect protocol resulted in abundance estimates that were close to the true value.
# 2: When applied to the non-perpendicular data, however, line transect analysis (conditional likelihood) resulted in overestimated of detection probabilities (~1) and underestimated abundance.
# 3: Using the ponit-transect likelihood formulation, sigma estimates are much closer to the true values. However, there still seems to be a small positive bias, which results in a slightly understimated abundance
# 4: I am still trying to undrstand why this bias occurs. Perhaps there is a correction to the likelihood formulation that could result in the correct estimates? Or maybe it the data-generating algorithm is creating biased sighting distance values? 



#########################################################
##### Simulate new data and analyze with 'unmarked' #####
#########################################################
library(unmarked)

## Simulate data. Function 'sim.sight.DS.unmkd' simulates detection and sighting distances for j = nsites, equally divided between 3 habitats. Detection is a half-normal function of distance. One can define covariate effects on detection for habitats 2 and 3 in parameters beta1 and beta2 (on the log scale). The distances are binned in distance categories according to the max distance B and interval.width parameters. The function also formats the data into unmarkedframeDS and unmarkedframeGDS (with 1 primary period) for use with distance sampling model functions in 'unmarked'. The output is a list with 3 unmarkedframeDS objects ("detection distance - line transect", "sighting distance - line transect", and "sighting distance - point transect" datasets) and 3 unmarkedframeGDS objects (same setup).
#obs.1: 'nsites' should be a multiple of three (so samples are equally divided among each habitat type). 


sim.sight.DS.unmkd <- function(N=100,sigma=30,B=100,nsites=30,beta1=0,beta2=0,interval.width=10){
  
  siteperhabitat <- nsites/3
  habitat <- factor(c(rep("hab1",siteperhabitat),rep("hab2",siteperhabitat),rep("hab3",siteperhabitat)),levels=c("hab1","hab2","hab3"))
  
  beta0 <- log(sigma) # beta for stream on log-scale (here equal to baseline sigma on log-scale)
  beta1 <- beta1 # Effect (slope) of pond habitat 
  beta2 <- beta2 # Effect (slope) of lagoon habitat 
  
  sigma0 <- exp(beta0) # sigma parameter for streams
  sigma1 <- exp(beta0 + beta1) # sigma for pond
  sigma2 <- exp(beta0 + beta2) # sigma for lagoon
  
  
  sim.ldata <- function(N.=N,sigma.=sigma,B.=B){
    
    u1 <- runif(N.,-B.,B.) # Distances to transect line for each object
    u2 <- runif(N.,0,500) # Perpendicular position along transect line of each object 
    
    g <- function(x, sig) exp(-x^2/(2*sig^2)) # The half-normal function
    
    ### Regular distance sampling protocol - perpendicular distances only ###
    p <- g(u1, sig=sigma.) # Detection probabilities for sightings at right angles
    y <- rbinom(N.,1,p)    # Some inds. are detected and their distance measured
    x <- abs(u1[y==1])    # Perpendicular distances of detected individuals
    
    
    # For objects located u2>100, detection is possible once the observer is within 100 radius (i.e. maximum detection distance = B = 100m)
    i <- which(u2>B.)
    # For objects located u2<100, maximum detection distance is slightly smaller than 100m, because the intersection of the detection arc to the transect line is before the point where the survey starts. 
    z <- which(u2<B.)
    u3 <- sqrt((u1[z]^2)+(u2[z]^2))  # The maximum detection distance for objects within the first 100m of the transect, for which detection is possible from the begining of the transect
    
    sight.x <- rep(NA,length(x)) # a vector to hold the "sighting distances"
    
    for(k in 1:length(x)){
      if(k%in%i){
        vector <- seq(x[k],B.,length.out = 100) # maximum detection distance is B
        sight.x[k] <- sample(x=vector,size=1,prob=g(vector,sig=sigma.))
      }  
      if(k%in%z){
        vector <- seq(x[k],u3[which(z==k)],length.out = 100) # maximum detection distance is slightly smaller than B (because the intersection of the detection arc to the transect line is before the point where the survey starts)
        sight.x[k] <- sample(x=vector,size=1,prob=g(vector,sig=sigma.))
      }
    }
    
    return(list(x=x,sight.x=sight.x))
  }
  
  det.distances <- matrix(NA,nsites,N)
  sight.distances <- matrix(NA,nsites,N)
  
  
  
  for(i in 1:siteperhabitat){ # a loop to simulate sampling for each site (computes values for hab1, hab2 and hab3 in each iteration, so i in 1:10)
    hab1 <- sim.ldata() # with baseline sigma value
    hab2 <- sim.ldata(sigma=sigma1) # with slope for pond habitat
    hab3 <- sim.ldata(sigma=sigma2) # with slope for lagoon habitat
    
    det.distances[i,1:length( hab1$x)] <-  hab1$x
    sight.distances[i,1:length( hab1$sight.x)] <-  hab1$sight.x
    
    det.distances[i+siteperhabitat,1:length(hab2$x)] <- hab2$x
    sight.distances[i+siteperhabitat,1:length(hab2$sight.x)] <- hab2$sight.x
    
    det.distances[i+(2*siteperhabitat),1:length(hab3$x)] <- hab3$x
    sight.distances[i+(2*siteperhabitat),1:length(hab3$sight.x)] <- hab3$sight.x
  }
  
  nbins<- B %/% interval.width # Number of distance categories
  det.bin <- det.distances %/% interval.width+1 # note integer division function %/%
  sight.bin <- sight.distances %/% interval.width+1 # note integer division function %/%
  
  if (length(sight.bin[which(sight.bin>10)])>0) sight.bin[which(sight.bin>10)]<-10
  
  det.obs <- matrix(NA,nsites,nbins)
  colnames(det.obs) <- 1:nbins
  sight.obs <- matrix(NA,nsites,nbins)
  colnames(sight.obs) <- 1:nbins
  
  for(z in 1:nsites){ # a loop to bin detection and sight distances on each site according to the distance categories
    det.freq <- table(det.bin[z,])
    det.padded <- rep(0,nbins) # Pad the frequencies to include those with 0 detections
    names(det.padded) <- 1:nbins
    det.padded[names(det.freq)]<-det.freq
    det.obs[z,] <- det.padded
    
    sight.freq <- table(sight.bin[z,])
    sight.padded <- rep(0,nbins) # Pad the frequencies to include those with 0 detections
    names(sight.padded) <- 1:nbins
    sight.padded[names(sight.freq)]<-sight.freq
    sight.obs[z,] <- sight.padded
  }
  
  det.obs <- cbind(as.data.frame(det.obs),habitat)
  sight.obs <- cbind(as.data.frame(sight.obs),habitat)
  
  ## Create umkarfedframe objects
  # line transect format for detection (perpendicular) distance data 
  det.lt.DS <- unmarkedFrameDS(y=as.matrix(det.obs[,1:nbins]),
                              siteCovs = as.data.frame(habitat),
                              dist.breaks = seq(0,B,by=nbins),
                              unitsIn="m",
                              survey="line",
                              tlength=rep(500,nsites)) # transect length in meters.  
  
  det.lt.GDS <- unmarkedFrameGDS(y=as.matrix(det.obs[,1:nbins]),
                               siteCovs = as.data.frame(habitat),
                               numPrimary = 1, # only one primary period!
                               dist.breaks = seq(0,B,by=nbins),
                               unitsIn="m",
                               survey="line",
                               tlength=rep(500,nsites))  
  
  # line transect format for sighting (non-perpendicular) distance data 
  sight.lt.DS <- unmarkedFrameDS(y=as.matrix(sight.obs[,1:nbins]),
                              siteCovs = as.data.frame(habitat),
                              dist.breaks = seq(0,B,by=nbins),
                              unitsIn="m",
                              survey="line",
                              tlength=rep(500,nsites)) # transect length in meters.  
  
  sight.lt.GDS<- unmarkedFrameGDS(y=as.matrix(sight.obs[,1:nbins]),
                               siteCovs = as.data.frame(habitat),
                               numPrimary = 1, # only one primary period!
                               dist.breaks = seq(0,B,by=nbins),
                               unitsIn="m",
                               survey="line",
                               tlength=rep(500,nsites))  
  
  # point transect format for sighting (non-perpendicular) distance data 
  sight.pt.DS <- unmarkedFrameDS(y=as.matrix(sight.obs[,1:nbins]),
                                siteCovs = as.data.frame(habitat),
                                dist.breaks = seq(0,B,by=nbins),
                                unitsIn="m",
                                survey="point")
  
  sight.pt.GDS <- unmarkedFrameGDS(y=as.matrix(sight.obs[,1:nbins]),
                                 siteCovs = as.data.frame(habitat),
                                 numPrimary = 1, # only one primary period!
                                 dist.breaks = seq(0,B,by=nbins),
                                 unitsIn="m",
                                 survey="point")
  
  return(list("det.lt.DS"=det.lt.DS,"det.lt.GDS"=det.lt.GDS,"sight.lt.DS"=sight.lt.DS,"sight.lt.GDS"=sight.lt.GDS,"sight.pt.DS"= sight.pt.DS," sight.pt.GDS"= sight.pt.GDS))
} 
#sim.dist <- sim.sight.DS.unmkd() # Simulate data with default values (N=20, sigma=30, 30 sites, no habitat effects on detection)
#fm1 <- distsamp(~1~1, data=sim.dist$det.umf1,output="abund",keyfun = "halfnorm") # fit distsamp line transect model on the detection (perpendicular) distance dataset.
#summary(fm1)


# Simulate 1000 datasets and see if estimates of abundance and sigma (scale parameter of the half-normal detection function) are close to true values for detection and sighting distance datasets, using the 'distsamp' line transect and point transect models. No covariate effects. Nhat = 20.
set.seed(2016)

simrep <- 1000
simout2 <- matrix(NA, nrow=simrep, ncol=6)
colnames(simout2) <- c("Nhat.det.lt","sigma.det.lt","Nhat.sight.lt","sigma.sight.lt","Nhat.sight.pt","sigma.sight")
Nhat = 25 # Average population per site. To avoid randomly generating 0 observations (which causes the simulation to crash), set the Nhat argument to at least 25.
sigma = 30

t1 <- proc.time()
for (sim in 1:simrep){ # Run time ~ 5 mins.

  sim.dist <- sim.sight.DS.unmkd(N=Nhat,nsites=60)
  
  fm1 <- distsamp(~1~1, data=sim.dist$det.lt.DS,output="abund",keyfun = "halfnorm")
  fm2 <- distsamp(~1~1, data=sim.dist$sight.lt.DS,output="abund",keyfun = "halfnorm")
  fm3 <- distsamp(~1~1, data=sim.dist$sight.pt.DS,output="abund",keyfun = "halfnorm")
  
  sigma1 <- exp(coef(fm1)[2])
  sigma2 <- exp(coef(fm2)[2])
  sigma3 <- exp(coef(fm3)[2])
  
  p1 <- integrate(gxhn, 0,B, sigma=sigma1)$value/B # Notice that here maximum detection distance B=100
  p2 <- integrate(gxhn, 0,B, sigma=sigma2)$value/B
  p3 <- integrate(gxhn, 0,B, sigma=sigma3)$value/B # And here, the detection function is integrated over B, not the point transect circle area.
  
  Nhat1 <- sum(sim.dist$det.lt.DS@y)/(p1*60)
  Nhat2 <- sum(sim.dist$sight.lt.DS@y)/(p2*60)
  Nhat3 <- sum(sim.dist$sight.pt.DS@y)/(p3*60)
  
  simout2[sim,1] <- Nhat1
  simout2[sim,2] <- sigma1
  simout2[sim,3] <- Nhat2
  simout2[sim,4] <- sigma2
  simout2[sim,5] <- Nhat3
  simout2[sim,6] <- sigma3
}
t2 <- proc.time()
print(t2 - t1)

simout2

####################################################################################
### Plot estimated abundance values from each method and compare with true value ###
####################################################################################

par(mfrow=c(2,3))

hist(simout2[,1], main="Detection Dist. Line Transect", xlab="Nhat", xlim=c(0,(2*Nhat)))
abline(v=Nhat, col=2)
hist(simout2[,3], main="Sight Dist. Line Transect", xlab="Nhat", xlim=c(0,(2*Nhat)))
abline(v=Nhat, col=2)
hist(simout2[,5], main="Sight Dist. Point Transect", xlab="Nhat", xlim=c(0,(2*Nhat)))
abline(v=Nhat, col=2)

boxplot(simout2[,2],ylim=c(0,(2*sigma)),main="Detection Dist. Line Transect", xlab="Sigma")
abline(h=sigma,col=2)
boxplot(simout2[,4],ylim=c(0,(2*sigma)),main="Sight Dist. Line Transect", xlab="Sigma")
abline(h=sigma,col=2)
boxplot(simout2[,6],ylim=c(0,(2*sigma)),main="Sight Dist. Point Transect", xlab="Sigma")
abline(h=sigma,col=2)


#############################################################################
### Include detection covariates - assess performance of slope estimators ###
#############################################################################
# Simulate 1000 datasets, in which factor Habitats has effects on the scale parameter sigma of the half-normal detection function. Beta1 (habitat 2) = -0.2, beta2 (habitat 3) = 0.4. Habitat 1 is the baseline level corresponding to the log(sigma) intercept estimate. Assess is coefficient estimates are close to true values.


set.seed(2016)

simrep <- 1000
simout3 <- matrix(NA, nrow=simrep, ncol=9)
colnames(simout3) <- rep(c("Det.lt","Sight.lt","Sight.pt"),3)
Nhat = 50 # Average population per site. We set this to a higher value because of negative betas for detection covariates.
sigma = 30

t1 <- proc.time()
for (sim in 1:simrep){ # Run time ~9 mins.
  
  sim.dist.cov <- sim.sight.DS.unmkd(N=Nhat,beta1=-0.2,beta2=0.4,nsites=60)
  
  fm4 <- distsamp(~habitat~1, data=sim.dist.cov$det.lt.DS,output="abund",keyfun = "halfnorm")
  fm5 <- distsamp(~habitat~1, data=sim.dist.cov$sight.lt.DS,output="abund",keyfun = "halfnorm")
  fm6 <- distsamp(~habitat~1, data=sim.dist.cov$sight.pt.DS,output="abund",keyfun = "halfnorm")
  
  simout3[sim,1] <- coef(fm4)[2]
  simout3[sim,2] <- coef(fm5)[2]
  simout3[sim,3] <- coef(fm6)[2]
  simout3[sim,4] <- coef(fm4)[3]
  simout3[sim,5] <- coef(fm5)[3]
  simout3[sim,6] <- coef(fm6)[3]
  simout3[sim,7] <- coef(fm4)[4]
  simout3[sim,8] <- coef(fm5)[4]
  simout3[sim,9] <- coef(fm6)[4]
}
t2 <- proc.time()
print(t2 - t1)

simout3

##############################################################################################
### Plot estimated intercept and coefficients from each method and compare with true value ###
##############################################################################################

par(mfrow=c(1,3))

boxplot(simout3[,1:3],ylim=c(3,4.5),main="log(sigma) - Intercept", xlab="Model")
abline(h=log(sigma),col=2)
boxplot(simout3[,4:6],ylim=c(-0.5,0.5),main="Beta 1 (Habitat 2)", xlab="Model")
abline(h=-0.2,col=2)
boxplot(simout3[,7:9],ylim=c(-2,8),main="Beta 2 (Habitat 3)", xlab="Model")
abline(h=0.4,col=2)


