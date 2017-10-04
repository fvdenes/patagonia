## Bayesian point-transect HDS with groupsize and DA ##
#Simulating dataset with habitat as a 3-level factor #
# ------------------------------------------------------------------------
simHDSg_factor <- function(type = "line", nsites = 100, lambda.group = 0.75, alpha0 = 0, alpha1 = 0.5, beta0 = 1, beta1 = -0.5,beta2 = 0.5, beta3=0.8, B = 4, discard0 = TRUE){
  #
  # Function simulates hierarchical distance sampling (HDS) data for groups under 
  #   either a line (type = "line") or a point (type = "point") transect protocol 
  #   and using a half-normal detection function (Buckland et al. 2001).
  #   Other function arguments:
  #     nsites: Number of sites (spatial replication)
  #     lambda.group: Poisson mean of group size
  #     alpha0, alpha1: intercept and slope of log-linear model relating sigma of
  #        half-normal detection function to group size
  #     beta0, beta1, etc: intercept and slope of log-linear model relating the Poisson
  #        mean of the number of groups per unit area to habitat
  #     B: strip half width
  #
  # Get covariates
  
  habitat_factor <- as.factor(sample(c(1,2,3),nsites,replace=T))  # Simulated habitat factor with 3 levels  
  habitat_matrix <- data.frame(model.matrix(~habitat_factor))             
  habitat2 <- habitat_matrix[,2] # habitat 1 is present when habitat2=0 and habitat3=0.
  habitat3 <- habitat_matrix[,3]
  elevation <- rnorm(nsites) 
  
  # Simulate abundance model for groups (Poisson GLM for N)
  lambda <- exp(beta0 + beta1*habitat2 + beta2*habitat3 + beta3*elevation)  # Density of groups per "square"
  N <- rpois(nsites, lambda)            # site-specific number of groups
  N.true <- N                           # for point: inside of B
  
  # Simulate observation model
  data <- groupsize <- NULL
  
  for(i in 1:nsites){
    if(N[i]==0){
      data <- rbind(data,c(i,NA,NA,NA,NA,NA)) # save site, y=1, u, v, d
      next
    }
    
    if(type=="line"){
      # Simulation of distances, uniformly, for each individual in the population
      d <- runif(N[i], 0, B)
      gs <- rpois(N[i],lambda.group) +1  # Observable group sizes >= 1
      groupsize<-c(groupsize,gs)
      sigma.vec <- exp(alpha0 + alpha1*(gs-1))  # Subtract 1 for interpretation
      # Detection probability for each group
      p <- exp(-d*d/(2*(sigma.vec^2)))
      # Determine if individuals are captured or not
      y <- rbinom(N[i], 1, p)
      u1 <- u2 <- rep(NA,N[i])
      # Subset to "captured" individuals only
      d <- d[y==1] ; u1 <- u1[y==1] ; u2 <- u2[y==1] ; gs <- gs[y==1] ; y <- y[y==1]
    }
    
    if(type=="point"){
      # Simulation of data on a circle of radius B (algorithm of Wallin)
      angle <- runif(N[i], 0, 360)
      r2 <- runif(N[i], 0, 1)
      r <-  B*sqrt(r2)
      u1 <-  r*cos(angle)  + B
      u2 <-  r*sin(angle)  + B
      
      d <- sqrt((u1 - B)^2 + (u2-B)^2)
      N.true[i] <- sum(d<= B)    # Population size inside of count circle, should be N[i] here.
      gs <- rpois(N[i], lambda.group) + 1
      groupsize <-c(groupsize,gs)
      sigma.vec <- exp(alpha0 + alpha1*(gs-1))
      # For counting individuals on a circle so we truncate p here
      p <- ifelse(d<(B), 1, 0)*exp(-d*d/(2*(sigma.vec^2)))
      y <- rbinom(N[i], 1, p)
      # Subset to "captured" individuals only
      d <- d[y==1] ; u1 <- u1[y==1] ; u2 <- u2[y==1] ; gs <- gs[y==1] ; y <- y[y==1]
    }
    # Now compile things into a matrix and insert NA if no individuals were 
    # captured at site i. Coordinates (u,v) are preserved.
    if(sum(y) > 0)
      data <- rbind(data,cbind(rep(i, sum(y)), y, u1, u2, d, gs))
    else
      data <- rbind(data,c(i,NA,NA,NA,NA,NA)) # make a row of missing data
  }
  # Subset to sites at which individuals were captured. You may or may not
  #  do this depending on how the model is formulated so be careful. 
  if(discard0)
    data <- data[!is.na(data[,2]),]
  
  # Visualisation
  if(type=="line"){       # For line transect
    par(mfrow = c(2,2))
    hist(data[,"d"], col = "lightblue", breaks = 20, main = 
           "Frequency of distances to groups", xlab = "Distance")
    ttt <- table(data[,1])
    n <- rep(0, nsites)
    n[as.numeric(rownames(ttt))] <- ttt
    plot(habitat_factor, n, main = "Observed group counts (n) vs. habitat", frame = F)
    plot(elevation, n, main = "Observed group counts (n) vs. elevation", frame = F)
    plot(table(data[,"gs"]), main = "Observed group sizes", ylab = "Frequency", frame = F)
  }
  
  if(type=="point"){       # For point transect
    par(mfrow = c(3,2))
    plot(data[,"u1"], data[,"u2"], pch = 16, main = 
           "Located groups in point transects", xlim = c(0, 2*B),
         ylim = c(0, 2*B), col = data[,1], asp = 1)
    points(B, B, pch = "+", cex = 3)
    library(plotrix)
    draw.circle(B, B, B)
    hist(data[,"d"], col = "lightblue", breaks = 20, main = 
           "Frequency of distances to groups", xlab = "Distance")
    ttt <- table(data[,1])
    n <- rep(0, nsites)
    n[as.numeric(rownames(ttt))] <- ttt
    plot(habitat_factor, n, main = "Observed group counts (n) vs. habitat", frame = F)
    plot(elevation, n, main = "Observed group counts (n) vs. elevation", frame = F)
    plot(table(data[,"gs"]), main = "Observed group sizes", ylab = "Frequency", frame = F)
  }
  
  # Output
  list(type = type, nsites = nsites, lambda.group = lambda.group, alpha0 = alpha0, alpha1 = alpha1, beta0 = beta0, beta1 = beta1,beta2=beta2, B = B, data=data, habitat_matrix=habitat_matrix, elevation=elevation,N = N, N.true = N.true, groupsize=groupsize)
}

# Analysis in BUGS - Model for point transect data #
# ------------------------------------------------------------------------
set.seed(1234)                 # we all create same data set
temp <- simHDSg_factor(nsites=100,lambda.group = 0.75,alpha0 = 0,type="point")

data <- temp$data              # harvest data
B <- temp$B                    # Get strip half width
habitat2 <- temp$habitat_matrix[,2] # habitat 1 is present when habitat2=0 and habitat3=0.
habitat3 <- temp$habitat_matrix[,3]
elevation <- temp$elevation
nsites <- temp$nsites          # Number of spatial replicates
groupsize <- data[,"gs"] -1    # Input groupsize-1 as data

delta <- 0.1 #width of distance bins for aproximation
d <- data[,5]
dclass <- d%/%delta+1 # assign observed distances to distance intervals
xg <- seq(0,B,delta)
midpt <- xg[-1] - delta/2 # get midpoint of distance intervals
nD <- length(midpt) # how many intervals

M <- 600                        # Size of augmented data set is M
nz <- M-nrow(data)              # Number of "pseudo-groups" added
y <- c(data[,2],rep(0,nz))      # Indicator of capture (== 1 for all obs. groups)
nind <- nrow(data)              # Number of observed groups
site <- c(data[,1], rep(NA,nz)) # Site they belong to is unknown 

groupsize <- c(groupsize, rep(NA,nz)) # .... as is their size
zst <- y                        # Starting values for data augmentation variable

dclass <- c(dclass,rep(NA,nz))

# Bundle data and produce summary
str(bugs.data <- list (y=y, B=B, nind=nind, nsites=nsites, midpt=midpt, delta=delta, dclass=dclass, habitat2=habitat2,habitat3=habitat3,elevation=elevation,site=site, nz=nz, nD=nD, groupsize=groupsize))



# Define model in BUGS langauge
cat("
    model{ 
    
    # Prior distributions for model parameters
    alpha0 ~ dunif(-10,10)               # groupsize intercept
    alpha1 ~ dunif(-10,10)               # groupsize slope
    beta0 ~ dunif(-10,10)                # intercept for log(lambda)
    beta1 ~ dunif(-10,10)                # slope for habitat 2
    beta2 ~ dunif(-10,10)                # slope for habitat 3
    beta3 ~ dunif(-10,10)                # slope for elevation
    lambda.group ~ dgamma(0.1, 0.1)
    
    # psi is a derived parameter
    psi <- sum(lambda[])/(nind+nz)
    
    # Individual level model: observations and process
    
    for(i in 1:(nind+nz)){
    
    z[i] ~ dbern(psi)                    # Data augmentation variables (site)
    dclass[i] ~ dcat(pi.probs[site[i],]) # Population distribution of distance class
    mu[i] <- z[i]*p[site[i],dclass[i]]   # p depends on site and distance class
    groupsize[i] ~ dpois(lambda.group)   # Group size is Poisson
    
    log(sigma[i]) <- alpha0 +  alpha1*groupsize[i]
    y[i] ~ dbern(mu[i])
    site[i] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
    zg[i]<- z[i]*(groupsize[i] + 1)      # Number of individuals in that group
    
    }
    
    for(s in 1:nsites){
    # Construct cell probabilities for nD cells (distance classes)
    for (g in 1:nD){                     # midpt[g] = midpoint of each cell; g is each distance class  
    log(p[s,g]) <- -midpt[g]*midpt[g]/(2*sigma[s]*sigma[s])
    pi[s,g] <- ((2*midpt[g])/(B*B))*delta # prob. per interval
    pi.probs[s,g] <- pi[s,g]/norm[s]
    f[s,g] <- p[s,g] * pi[s,g]
    fc[s,g] <- f[s,g]/pcap[s]          # Conditional probabilities
    }
    pcap[s] <- sum(f[s,])                # Capture prob. is the sum of all rectangular areas
    norm[s] <- sum(pi[s,])
    
    # Model for population size of groups
    N[s] ~ dpois(lambda[s])
    log(lambda[s])<- beta0 + beta1*habitat2[s] + beta2*habitat3[s] + beta3*elevation[s]
    site.probs[s]<- lambda[s]/sum(lambda[])
    }
    
    # Derived quantities
    G <- sum(z[])        # Total number of groups
    Ntotal <- sum(zg[])  # Total population size (all groups combined)
    }
    ",fill=TRUE, file="model1.txt")

# Load some libraries, define MCMC settings, inits function and parameters to save

library("jagsUI")  # 
ni <- 8000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3
inits <- function(){list(alpha0=0, alpha1=0.5, beta0=0.5, beta1=-0.5,beta2=0.8,beta3=0.8, z=zst)}
params <- c("alpha0", "alpha1", "beta0", "beta1","beta2","beta3", "psi", "Ntotal", "G", "lambda.group")

# Call JAGS, check convergence and summarize posterior distributions
out1 <- jags(bugs.data, inits, params, "model1.txt", n.thin=nt,  n.chains=nc, n.burnin=nb,n.iter=ni, parallel=FALSE)

par(mfrow = c(3,2))
traceplot(out1)
print(out1, 3)

# Data-generating values for model parameters: alpha0 = 0, alpha1 = 0.5, beta0 = 1, beta1 = -0.5,beta2 = 0.5, beta3 = 0.8

sum(temp$N.true) # = number of groups (G)
sum(temp$groupsize) # = number of individuals (Ntotal)
