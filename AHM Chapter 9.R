# =========================================================================
#
# 9. Advanced Hierarchical Distance Sampling 
#
# =========================================================================





# 9.1 Introduction
# ------------------------------------------------------------------------



# 9.2 Distance sampling with clusters, groups, or other individual covariates
# ------------------------------------------------------------------------


# 9.2.1 Simulating HDS data with group size
# ------------------------------------------------------------------------
# Function to simulate data under HDS protocol with groups
simHDSg <- function(type = "line", nsites = 100, lambda.group = 0.75, alpha0 = 0, alpha1 = 0.5, beta0 = 1, beta1 = 0.5, B = 4, discard0 = TRUE){
  #
  # Function simulates hierarchical distance sampling (HDS) data for groups under 
  #   either a line (type = "line") or a point (type = "point") transect protocol 
  #   and using a half-normal detection function (Buckland et al. 2001).
  #   Other function arguments:
  #     nsites: Number of sites (spatial replication)
  #     lambda.group: Poisson mean of group size
  #     alpha0, alpha1: intercept and slope of log-linear model relating sigma of
  #        half-normal detection function to group size
  #     beta0, beta1: intercept and slope of log-linear model relating the Poisson
  #        mean of the number of groups per unit area to habitat
  #     B: strip half width
  #
  # Get covariates
  habitat <- rnorm(nsites)              # Simulated covariate  
  
  # Simulate abundance model for groups (Poisson GLM for N)
  lambda <- exp(beta0 + beta1*habitat)  # Density of groups per "square"
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
    par(mfrow = c(1, 3))
    hist(data[,"d"], col = "lightblue", breaks = 20, main = 
           "Frequency of distances to groups", xlab = "Distance")
    ttt <- table(data[,1])
    n <- rep(0, nsites)
    n[as.numeric(rownames(ttt))] <- ttt
    plot(habitat, n, main = "Observed group counts (n) vs. habitat", frame = F)
    plot(table(data[,"gs"]), main = "Observed group sizes", ylab = "Frequency", frame = F)
  }
  
  if(type=="point"){       # For point transect
    par(mfrow = c(2,2))
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
    plot(habitat, n, main = "Observed group counts (n) vs. habitat", frame = F)
    plot(table(data[,"gs"]), main = "Observed group sizes", ylab = "Frequency", frame = F)
  }
  
  # Output
  list(type = type, nsites = nsites, lambda.group = lambda.group, alpha0 = alpha0, alpha1 = alpha1, beta0 = beta0, beta1 = beta1, B = B, data=data, habitat=habitat, N = N, N.true = N.true, groupsize=groupsize)
}


data <- simHDSg(type = "line")     # Defaults for line transect data
data <- simHDSg(type = "point")    # Default for point transect data
data <- simHDSg(lambda.group = 5)  # Much larger groups
data <- simHDSg(lambda.group = 5, alpha1 = 0) # No effect of groups size on p


# 9.2.2 Analysis in BUGS
# ------------------------------------------------------------------------
set.seed(1234)                 # we all create same data set
temp <- simHDSg(type="line")   # Execute function
data <- temp$data              # harvest data
B <- temp$B                    # Get strip half width
habitat <- temp$habitat        # habitat covariate
nsites <- temp$nsites          # Number of spatial replicates
groupsize <- data[,"gs"] -1    # Input groupsize-1 as data

M <- 400                        # Size of augmented data set is M
nz <- M-nrow(data)              # Number of "pseudo-groups" added
y <- c(data[,2],rep(0,nz))      # Indicator of capture (== 1 for all obs. groups)
nind <- nrow(data)              # Number of observed groups
site <- c(data[,1], rep(NA,nz)) # Site they belong to is unknown 
d <- c(data[,5], rep(NA,nz))    # Their distance data are missing ...
groupsize <- c(groupsize, rep(NA,nz)) # .... as is their size
zst <- y                        # Starting values for data augmentation variable

# Bundle data and produce summary
str(bugs.data <- list (y=y, B=B, nind=nind, nsites=nsites, d=d, habitat=habitat, 
                       site=site, nz=nz, groupsize=groupsize))


# Define model in BUGS langauge
cat("
    model{ 
    
    # Prior distributions for model parameters
    alpha0 ~ dunif(-10,10)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-10,10)
    beta1 ~ dunif(-10,10)
    lambda.group ~ dgamma(0.1, 0.1)
    # psi is a derived parameter
    psi <- sum(lambda[])/(nind+nz)
    
    # Individual level model: observations and process
    for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)                   # Data augmentation variables
    d[i] ~ dunif(0, B)                  # Distance is uniformly distributed
    groupsize[i] ~ dpois(lambda.group)  # Group size is Poisson
    
    log(sigma[i]) <- alpha0 +  alpha1*groupsize[i]
    mu[i] <- z[i]*exp(-d[i]*d[i]/(2*sigma[i]*sigma[i])) #p dep on dist class 
    # here using the half normal detection function (Buckland et al. 2001)
    y[i] ~ dbern(mu[i])
    
    site[i] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
    zg[i]<- z[i]*(groupsize[i] + 1)      # Number of individuals in that group
    }
    
    for(s in 1:nsites){
    # Model for population size of groups
    N[s] ~ dpois(lambda[s])
    log(lambda[s])<- beta0 + beta1*habitat[s]
    site.probs[s]<- lambda[s]/sum(lambda[])
    }
    
    # Derived quantities
    G <- sum(z[])        # Total number of groups
    Ntotal <- sum(zg[])  # Total population size (all groups combined)
    }
    ",fill=TRUE, file="model1.txt")

# Load some libraries, define MCMC settings, inits function and parameters to save
library("R2WinBUGS")
library("jagsUI")  # 
ni <- 6000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3
inits <- function(){list(alpha0=0, alpha1=0.5, beta0=0, beta1=0, z=zst)}
params <- c("alpha0", "alpha1", "beta0", "beta1", "psi", "Ntotal", "G", 
            "lambda.group")

# Call JAGS (ART 1.4 min), check convergence and summarize posterior distributions
out1 <- jags(bugs.data, inits, params, "model1.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb,n.iter=ni)
traceplot(out1)   ;   print(out1, 3)


# 9.2.3 Imperfect observation of cluster size
# ------------------------------------------------------------------------



# 9.3 Time-removal and distance sampling combined
# ------------------------------------------------------------------------


# 9.3.1 The four-part hierarchical model
# ------------------------------------------------------------------------


# 9.3.2 Simulating some time-removal/DS data
# ------------------------------------------------------------------------
# Obtain a data set and harvest the results
set.seed(1235)                 # so we all create the same data set 
temp <- simHDStr(type="point") # Simulate point count-removal data set
data <- temp$data              # harvest data
B <- temp$B                    # upper limit of counting (maximum count distance)
nsites <- temp$nsites          # Number of sites
habitat <- temp$habitat        # habitat covariate
K <- temp$K                    # Number of removal periods


# Create the observed encounter frequencies per site (include the zeros! )
data <- data[!is.na(data[,2]),]   # Sites where detections did occur
n <- rep(0,nsites)                # The full site vector
names(n) <- 1:nsites
n[names(table(data[,1]))] <- table(data[,1])  # Put in the counts
site <- data[,1]
nobs <- nrow(data) 

# Create the distance class data
nD <- 10             # Number of distance classes 
delta <- B/nD        # bin size or width
mdpts <- seq(delta/2,B,delta) # midpoint distance of bins up to max distance
dclass <- data[,"d"] # distance class for each observation
dclass <- dclass%/%delta  +1
tint <- data[,"aux"]

# Bundle data and summarize
str( win.data<-list(n=n, site=site, dclass=as.numeric(dclass),nsites=nsites, nobs=nobs, delta=delta, nD=nD,mdpts=mdpts,B=B, K=K, tint=tint, habitat=habitat) )


cat("
    model {
    # Prior distributions for basic parameters
    # Intercepts
    beta.a0 ~ dnorm(0,0.01)    # intercept for availability
    alpha0 ~ dnorm(0, 0.01)    # intercept for sigma
    alpha1 ~ dnorm(0,0.01)     # slope on sigma covariate
    # Coefficients
    # beta.a1 ~ dnorm(0,0.01)  # slope for availability covariate
    beta0 ~ dnorm(0,0.01)      # intercept for lambda
    beta1 ~dnorm(0,0.01)       # slope for lambda covariate
    
    for(s in 1:nsites){
    # Add covariates to scale parameter DISTANCE (perceptibility)
    log(sigma[s]) <- alpha0 +  alpha1*habitat[s] 
    # Add covariates for availability here TIME-REMOVAL (availability)
    p.a[s] <- exp(beta.a0) / (1+exp(beta.a0)) 
    # Optional covariates on availability
    # exp(beta.a0 + beta.a1*date[s])/(1+exp(beta.a0+beta.a1*date[s]))
    # Distance sampling detection probability model
    for(b in 1:nD){
    log(g[b,s]) <- -mdpts[b]*mdpts[b]/(2*sigma[s]*sigma[s])  # Half-normal  
    f[b,s] <- ( 2*mdpts[b]*delta )/(B*B) # Radial density function
    pi.pd[b,s] <- g[b,s]*f[b,s]  #  Product Pr(detect)*Pr(distribution)
    pi.pd.c[b,s] <- pi.pd[b,s]/pdet[s]  # Conditional probabilities
    }
    pdet[s] <- sum(pi.pd[,s])  # Probability of detection at all 
    
    # Time-removal probabilities
    for (k in 1:K){
    pi.pa[k,s] <- p.a[s] * pow(1-p.a[s], (k-1))  
    pi.pa.c[k,s] <- pi.pa[k,s]/phi[s] # Conditional probabilities of availability
    }
    phi[s] <- sum(pi.pa[,s]) # Probability of ever available
    }
    # Conditional observation model for categorical covariates
    for(i in 1:nobs){  
    dclass[i] ~ dcat(pi.pd.c[,site[i]]) 
    tint[i] ~ dcat(pi.pa.c[,site[i]])
    }
    # Abundance model
    for(s in 1:nsites){ 
    # Binomial model for # of captured individuals
    # n[s] ~ dbin(pmarg[s], M[s]) # Formulation b, see text
    # pmarg[s] <- pdet[s]*phi[s] 
    n[s] ~ dbin(pdet[s], N[s])    # Formulation a, see text
    N[s] ~ dbin(phi[s],M[s])      # Number of available individuals
    M[s] ~ dpois(lambda[s])       # Abundance per survey/site/point
    # Add site-level covariates to lambda
    log(lambda[s]) <- beta0 + beta1*habitat[s] 
    }
    # Derived quantities
    Mtot <- sum(M[])  # Total population size
    Ntot <- sum(N[])  # Total available population size
    PDETmean <- mean(pdet[]) # Mean perceptibility across sites
    PHImean <- mean(phi[]) # Mean availability across sites
    }
    ", fill=TRUE, file="tr-ds.txt")

# Create initial values (including for M and N) and list parameters to save
Mst <- Nst <- n + 1
inits <- function(){list(M=Mst, N=Nst, alpha0=1, beta0=runif(1,-1,1), 
                         beta.a1=runif(1,-1,1), beta1=runif(1,-1,1), alpha1=runif(1,-1,1),
                         beta.a0=runif(1,-1,1))}
params <- c("beta.a0", "beta.a1", "alpha0", "alpha1", "beta0", "beta1", "PDETmean", "PHImean", "Mtot", "Ntot")

# MCMC settings
ni <- 50000   ;   nb <- 10000   ;   nt <- 4   ;   nc <- 3

# Run JAGS in parallel (ART 7.3 min), check convergence and summarize posteriors
out2a <- jags(data=win.data, inits=inits, parameters=params, 
              model.file ="tr-ds.txt",n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, 
              parallel = TRUE)
traceplot(out2a)   ;   print(out2a, 3)

sum(temp$M) 

print(out2b,3)



# 9.4 Mark-Recapture/Double observer Distance Sampling
# ------------------------------------------------------------------------


# 9.4.1 Simulating MRDS data
# ------------------------------------------------------------------------
# Simulate a double-observer sampling data set
set.seed(1235)
temp <- simHDStr(type="point", method="double") # simulate double observer point count data set
data <- temp$data         # harvest data
B <- temp$B               # upper limit of counting (maximum count distance)
nsites <-temp$nsites      # number of sites
habitat <-temp$habitat    # habitat covariate


# Processing of the data: pad the count vector with 0s etc.
data <- data[!is.na(data[,2]),]
n <- rep(0,nsites)
names(n) <- 1:nsites
n[names(table(data[,1]))] <- table(data[,1])
site <- data[,1]
dclass <- data[,"d"]      # categorical distance class for each observation
aux <- data[,"aux"]       # the auxiliary variable is capture history

# Create the categorical distance variable, use 10 classes here.
nD <- 10
delta <- B/nD # bin width
mdpts <-seq(delta/2,B,delta) # midpoint of bins up to max distance
nobs <- nrow(data) 
dclass <- dclass%/%delta  +1

# Bundle data and look at overview of data
str( win.data <-list(n=n,site=site, dclass=as.numeric(dclass), nsites=nsites, 
                     nobs=nobs, delta=delta, nD=nD, mdpts=mdpts, B=B, aux=aux, habitat=habitat) )


# 9.4.2 Analysis in BUGS
# ------------------------------------------------------------------------
# Define model in BUGS langauge
cat("
    model {
    
    #Priors for fixed detection parameters
    # 2 observer detection probability parameters
    logitp1 ~ dnorm(0, 0.01)
    logitp2 ~ dnorm(0, 0.01)
    # Intercepts
    alpha0 ~ dnorm(0, 0.01)    # intercept for sigma
    alpha1 ~ dnorm(0, 0.01)    # slope on sigma covariate
    # Coefficients
    beta0 ~ dnorm(0,0.01)      # intercept for lambda
    beta1 ~ dnorm(0,0.01)      # slope for lambda covariate
    
    # Detection scale parameter model
    for(s in 1:nsites){ 
    # Covariates on scale parameter (perceptibility)
    log(sigma[s]) <- alpha0 + alpha1*habitat[s] 
    # Double observer cell probabilities here if there are covariates
    logit(pobs[1,s]) <- logitp1 # + covariates
    logit(pobs[2,s]) <- logitp2 # + covariates
    
    # Distance sampling model and cell probabilities
    for(b in 1:nD){
    log(g[b,s]) <- -mdpts[b]*mdpts[b]/(2*sigma[s]*sigma[s])  # half-normal 
    f[b,s] <- ( 2*mdpts[b]*delta )/(B*B) # Scaled radial density function
    pi.pd[b,s] <- g[b,s]*f[b,s]          #  Product Pr(detect)*Pr(distribution)
    pi.pd.c[b,s] <- pi.pd[b,s]/pdet[s]   # Conditional cell probabilities
    }
    pdet[s] <- sum(pi.pd[,s])              # Marginal probability of detection
    
    # Double observer cell probabilities and conditional probabilities
    doprobs[1,s] <- pobs[1,s]*(1-pobs[2,s]) 
    doprobs.condl[1,s] <- doprobs[1,s]/sum(doprobs[,s])
    doprobs[2,s] <- (1-pobs[1,s])*pobs[2,s]
    doprobs.condl[2,s] <- doprobs[2,s]/sum(doprobs[,s])
    doprobs[3,s] <- pobs[1,s]*pobs[2,s] 
    doprobs.condl[3,s] <- doprobs[3,s]/sum(doprobs[,s])
    pavail[s] <- sum(doprobs[,s])  # probability of availability AT ALL
    }
    
    # Observation model for two categorical covariates
    for(i in 1:nobs){  
    dclass[i] ~ dcat(pi.pd.c[,site[i]]) 
    aux[i] ~ dcat(doprobs.condl[,site[i]])
    }
    
    # Abundance model
    for(s in 1:nsites){ 
    # Binomial model for # of captured individuals
    n[s] ~ dbin(pdet[s], N[s]) 
    N[s] ~ dbin(pavail[s], M[s])   # binomial availability model
    # Abundance model
    M[s] ~ dpois(lambda[s])        # predicted abundance per survey/site/point
    # Add site-level covariates to lambda
    log(lambda[s])<- beta0 + beta1*habitat[s] 
    }
    # Derived parameters
    Mtot <- sum(M[])
    Ntot <- sum(N[])
    logit(p1) <- logitp1
    logit(p2) <- logitp2
    sigma0 <- exp(alpha0)            # Baseline sigma
    }
    ", fill=TRUE,file="do_model.txt")

# Inits function
Nst <- n + 1         # inits for N
inits <- function(){list(M=Nst+1, N=Nst, alpha0=runif(1,1,2),
                         beta0=runif(1,-1,1), beta1=runif(1,-1,1), alpha1=runif(1,-1,1),
                         logitp1=0, logitp2=0)} 

# Parameters to monitor
params <- c("alpha0", "alpha1", "beta0", "beta1", "Ntot", "Mtot", "logitp1", 
            "logitp2", "p1", "p2", "sigma0")

# MCMC settings
ni <- 50000   ;   nb <- 10000   ;   nt <- 4   ;   nc <- 3

# Run JAGS in parallel (ART 6.8 min), check convergence and summarize the results
out3 <- jags(data=win.data, inits=inits, parameters.to.save=params, 
             model.file="do_model.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
             parallel = TRUE)
traceplot(out3)   ;   print(out3, 3)


# Put true values into a vector
truth <- temp$parms
psi <- 1-(1-truth["p.double1"])*(1-truth["p.double2"])# Compute availability
truth <- c(truth[c("p.double1", "p.double2")], exp(truth["alpha0"]), 
           truth["beta0"], "Mtot" = sum(temp$M),
           "Ntot" = sum(temp$M)*as.numeric(psi), 
           truth[c("alpha0","alpha1","beta1")])

# Get posterior means and 2.5% and 97.5% percentiles (95% CRI)
post <- out3$summary[c("p1", "p2", "sigma0", "beta0", "Mtot", "Ntot", 
                       "alpha0", "alpha1", "beta1"), c(1,3,7)]

# Table compares truth with posterior mean and 95% CRI from JAGS
cbind(truth, posterior = round(post, 3))


# 9.4.3 Remarks
# ------------------------------------------------------------------------



# 9.5 Open HDS models: temporary emigration
# ------------------------------------------------------------------------


# 9.5.1 Data and model structure
# ------------------------------------------------------------------------


# 9.5.2 Cautionary note on temporary emigration processes
# ------------------------------------------------------------------------


# 9.5.3 Modeling temporary emigration with distance sampling in unmarked using the function gdistsamp
# ------------------------------------------------------------------------
# Load the wagtail data, investigate NA patterns in detection data
load("wagtail.RData")
Y <- wagtail$Y
table(n.missing <- apply(y, 1, function(x) sum(is.na(x)))) # Frequency distribution of number of NAs per site
n.missing
keep <- which(n.missing == 0)     # Sites with complete distance data
Y <- Y[keep,]                     # restrict analysis to those

# Harvest other data for sites with complete distance data
potato <- wagtail$potato[keep]   ;   grass <- wagtail$grass[keep]
lscale <- wagtail$lscale[keep]   ;   hour <- wagtail$hour[keep,]
date <- wagtail$date[keep,]   ;   rep <- wagtail$rep[keep,]
breaks <- wagtail$breaks

# Look at the distance data
str(Y)
tmp <- apply(Y, 2, sum, na.rm = T) 
matplot(1:6, t(matrix(tmp, nrow = 4, byrow= T)), type = "b", ylim = c(0, 90), xlab = "Distance class", ylab = "Number of wagtails", frame = F, lwd = 3, lty = 1)

# Standardize all continuous covariates
mn.potato <- mean(potato)   ;   sd.potato <- sd(potato)
mn.grass <- mean(grass)   ;   sd.grass <- sd(grass)
mn.lscale <- mean(lscale)   ;   sd.lscale <- sd(lscale)
mn.date <- mean(date)    ;   sd.date <- sd(c(date))
mn.hour <- mean(hour)   ;   sd.hour <- sd(c(hour))
POTATO <- (potato - mn.potato) / sd.potato
GRASS <- (grass - mn.grass) / sd.grass
LSCALE <- (lscale - mn.lscale) / sd.lscale
DATE <- (date - mn.date) / sd.date
HOUR <- (hour - mn.hour) / sd.hour

# Package into unmarked GDS data frame and inspect the data
umf <- unmarkedFrameGDS(y = Y[,1:24], survey="point", unitsIn="m",
                        dist.breaks=breaks, numPrimary = 4,
                        siteCovs = data.frame(POTATO, GRASS, LSCALE),
                        yearlySiteCovs=list(rep = rep, DATE = DATE, HOUR = HOUR))
str(umf)
summary(umf)


# Model fitting: Null models fm0
# exponential detection function
summary(fm0.exp <- gdistsamp(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
                             keyfun = "exp", output = "density", unitsOut = "ha", 
                             mixture = "P", K = 100, se = TRUE, data = umf) )

# hazard detection function
summary(fm0.haz <- gdistsamp(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
                             keyfun = "haz", output = "density", unitsOut = "ha", 
                             mixture = "P", K = 100, se = TRUE, data = umf ) )


# half-normal detection function 
summary(fm0.hn <- gdistsamp(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
                            keyfun = "halfnorm", output = "density", unitsOut = "ha",
                            mixture = "P", K = 100, se = TRUE, data = umf,control=list(trace=TRUE, REPORT=1)) )

# Compare AIC scores for 3 detection functions
rbind('AIC exp' = fm0.exp@AIC, 'AIC haz' = fm0.haz@AIC, 'AIC hn' = fm0.hn@AIC)

backTransform(fm0.haz, type="lambda")

backTransform(fm0.haz, type="det")

plot(1:300, gxhaz(1:300, shape = exp(5.13), scale=1.59), frame = F, type = "l", xlab = "Distance WagtailñObserver (metres)", ylab = "Detection probability", lwd=3)


# Model with time-dependent phi
fm1 <- gdistsamp(lambdaformula = ~1, phiformula = ~rep-1, 
                 pformula = ~1, keyfun = "haz", output = "density", unitsOut = "ha",
                 mixture = "P", K = 100, se = TRUE, data = umf)

# Compare AIC for models with phi constant and phi time-dependent
rbind('AIC phi constant' = fm0.haz@AIC, 'AIC phi time-dep' = fm1@AIC)


# Add covariates on lambda:  2-phase fitting to assist convergence using
#   first K = 20
summary( fm2.init <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
                               phiformula = ~rep-1, pformula = ~1, keyfun = "haz", output = "density",
                               unitsOut = "ha",   control=list(trace=TRUE, REPORT=1),
                               mixture = "P", K = 20, se = TRUE, data = umf))

starts <- coef(fm2.init)
summary( fm2  <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
                           phiformula = ~rep-1, pformula = ~1, keyfun = "haz", output = "density",
                           unitsOut = "ha",  starts=starts, control=list(trace=TRUE, 
                                                                         REPORT=1), mixture = "P", K = 100, se = TRUE, data = umf))


# Add covariates on lambda: optimisation with default, no inits
summary( fm2a <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
                           phiformula = ~rep-1, pformula = ~1, keyfun = "haz", output = "density",
                           unitsOut = "ha", control=list(trace=TRUE, REPORT=1),
                           mixture = "P", K = 100, se = TRUE, data = umf))


# Estimates of availability probabilities
plogis(fm2@estimates@estimates$phi@estimates)

# Models with time-dependent phi, AIC-best key functions and 
#    three covariates on phi. Use previous estimates as starting values.
(tmp <- coef(fm2))
starts <- c(tmp[1:4], tmp[5:8], 0,0,0, tmp[9], tmp[10])

summary(fm3 <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
                         phiformula = ~(rep-1)+ POTATO+GRASS+LSCALE, pformula = ~1,
                         keyfun = "haz", output = "density", unitsOut = "ha",
                         mixture = "P", K = 100, control=list(trace=TRUE, REPORT=1),
                         se = TRUE, data = umf, starts = starts))

# Models with time-dependent phi, AIC-best key function, 3 covariates on phi 
#     and, in addition, date and hour on detection
# linear effects on detection
tmp <- fm3@estimates@estimates
starts <- c(tmp$lambda@estimates, tmp$phi@estimates, tmp$det@estimates, 0, 0, tmp$scale@estimates)
summary(fm4A <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
                          phiformula = ~(rep-1)+ POTATO+GRASS+LSCALE, pformula = ~ DATE + HOUR,
                          keyfun = "haz", output = "density", unitsOut = "ha",
                          mixture = "P", K = 100, control=list(trace=TRUE, REPORT=1),
                          se = TRUE, data = umf, starts = starts) )

# quadratic effects on detection
tmp <- fm4A@estimates@estimates
p.start <- tmp$det@estimates
p.start <- c(p.start[1:2], 0, p.start[3], 0)
starts <- c(tmp$lambda@estimates, tmp$phi@estimates, p.start, tmp$scale@estimates)

summary(fm4B <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
                          phiformula = ~(rep-1)+ POTATO+GRASS+LSCALE, 
                          pformula = ~ DATE + I(DATE^2) + HOUR + I(HOUR^2),
                          keyfun = "haz", output = "density", unitsOut = "ha",
                          mixture = "P", K = 100, control=list(trace=TRUE, REPORT=1),
                          se = TRUE, data = umf, starts = starts) )

starts <- coef(fm4B)[-16]   # Drop coef for HOUR^2
summary(fm4C <- gdistsamp(~ POTATO+GRASS+LSCALE, 
                          ~(rep-1)+ POTATO+GRASS+LSCALE, ~ DATE + I(DATE^2) + HOUR,
                          keyfun = "haz", output = "density", unitsOut = "ha",
                          mixture = "P", K = 100, control=list(trace=TRUE, REPORT=1),
                          se = TRUE, data = umf, starts = starts) )

starts <- c(coef(fm4C), 0)
summary(fm5 <- gdistsamp(~ POTATO+GRASS+LSCALE, 
                         ~(rep-1)+ POTATO+GRASS+LSCALE, ~ DATE + I(DATE^2) + HOUR,
                         keyfun = "haz", output = "density", unitsOut = "ha",
                         mixture = "NB", K = 100, control=list(trace=TRUE, REPORT=1),
                         se = TRUE, data = umf , starts = starts) )

# Now we create a model selection table of these various models
modSel(fitList(fm0.haz, fm1, fm2, fm3, fm4A, fm4B, fm4C, fm5) )

summary(fm5)

# Bootstrap Goodness-of-fit assessment: ART ~ 20 hours
set.seed(1234)
(pb <- parboot(fm5, fitstats, nsim=100, report=5))

# Compute magnitude of "overdispersion" c.hat as ratio of observed to expected 
#    chisquare test statistic
(c.hat <- pb@t0[2] / mean(pb@t.star[,2]))  # c-hat as ratio of observed/expected 


# Predictions of lambda for POTATO, GRASS and LSCALE
newdat1 <- data.frame(POTATO=0, GRASS=0, LSCALE = seq(-1.8,4.33,,100))
newdat2 <- data.frame(POTATO=seq(-0.75,3,,100), GRASS=0, LSCALE = 0)
newdat3 <- data.frame(POTATO=0, GRASS=seq(-0.4, 3.6,,100), LSCALE = 0)
pred1 <- predict(fm5, type="lambda", newdata=newdat1, append = T)
pred2 <- predict(fm5, type="lambda", newdata=newdat2, append = T)
pred3 <- predict(fm5, type="lambda", newdata=newdat3, append = T)

# Predictions of phi for POTATO, GRASS and LSCALE and for rep = 1
newdat4 <- data.frame(rep = factor('1', levels = c('1','2','3','4')), POTATO=0, GRASS=0, LSCALE = seq(-1.8,4.33,,100))
newdat5 <- data.frame(rep = factor('1', levels = c('1','2','3','4')), POTATO=seq(-0.75,3,,100), GRASS=0, LSCALE = 0)
newdat6 <- data.frame(rep = factor('1', levels = c('1','2','3','4')), POTATO=0, GRASS=seq(-0.4, 3.6,,100), LSCALE = 0)
pred4 <- predict(fm5, type="phi", newdata=newdat4, append = T)
pred5 <- predict(fm5, type="phi", newdata=newdat5, append = T)
pred6 <- predict(fm5, type="phi", newdata=newdat6, append = T)

# Predictions of detection function sigma for DATE and HOUR
newdat7 <- data.frame(DATE = seq(-1.51,1.69,,100), HOUR = 0)
newdat8 <- data.frame(DATE=0, HOUR = seq(-1.92,3.1,,100))
pred7 <- predict(fm5, type="det", newdata=newdat7, append = T)
pred8 <- predict(fm5, type="det", newdata=newdat8, append = T)

par(mfrow = c(1,3), mar = c(5,5,2,2), cex.lab = 1.5, cex.axis = 1.5)
plot(newdat1$LSCALE, pred1[,1], xlab="Standardized covariate", ylab="Density (birds/ha)", lwd=3,type="l", frame = F)
lines(newdat2$POTATO, pred2[,1], lwd=3, col="red")
lines(newdat3$GRASS, pred3[,1], lwd=3, col="blue")
legend(-1.6, 1.65, c("LSCALE", "POTATO", "GRASS"), col=c("black", "red", "blue"), lty=1, lwd=3, cex=1.2) 

plot(newdat4$LSCALE, pred4[,1], xlab="Standardized covariate", ylab="Availability (phi)", lwd=3,type="l", frame = F)
lines(newdat5$POTATO, pred5[,1], lwd=3, col="red")
lines(newdat6$GRASS, pred6[,1], lwd=3, col="blue")
legend(2, 0.65, c("LSCALE", "POTATO", "GRASS"), col=c("black", "red", "blue"), lty=1, lwd=3, cex=1.2) 

plot(newdat7$DATE, pred7[,1], xlab="Standardized covariate", ylab="Detection function (sigma)", lwd=3,type="l", frame = F, ylim = c(100, 200))
lines(newdat8$HOUR, pred8[,1], lwd=3, col="red")
legend(0.5, 140, c("DATE", "HOUR"), col=c("black", "red"), lty=1, lwd=3, cex=1.2) 



# 9.5.4 Fitting temporary emigration HDS models in BUGS
# ------------------------------------------------------------------------


# 9.5.4.1 Simulating a temporary emigration system
# ------------------------------------------------------------------------
simHDSopen(type="line", nsites = 100, mean.lam = 2, beta.lam = 0, mean.sig = 1, beta.sig = 0, B = 3, discard0=TRUE, nreps=2, phi=0.7, nyears=5, beta.trend = 0)

# Obtain a temporary emigration data set 
set.seed(1234)
str(tmp <- simHDSopen("point", nreps=7, nyears=5, nsites=100)  ) 
attach(tmp)


apply(tmp$M.true,2,sum)

# Define distance class information
delta <- 0.5
nD <- B%/%delta                 # Number of distance classes
midpt <- seq(delta/2, B, delta) # Mid-point of distance intervals

# Create the 4-d array
y4d <- array(0,dim=c(nsites, nD, K, nyears))
for(yr in 1:nyears){
  for(rep in 1:K){
    data <- tmp$data[[yr]][[rep]]
    site <- data[,1]
    dclass <- data[,"d"]%/%delta + 1
    ndclass <- B%/%delta
    dclass <- factor(dclass, levels= 1:ndclass)
    y4d[1:nsites,1:nD,rep,yr] <- table(site, dclass)
  }
}

y3d <- y4d[,,,1]


# Bundle and summarize the data set
nobs <- apply(y3d, c(1,3), sum)  # Total detections per site and occasion
str( data <- list(y3d=y3d, nsites=nsites, K=K, nD=nD, midpt=midpt, delta=delta, habitat=habitat, B=B, nobs = nobs) )


# Define model in BUGS
cat("
    model {
    # Prior distributions
    beta0 ~ dnorm(0, 0.01)  # Intercept for log(lambda)
    mean.lam <- exp(beta0)
    beta1 ~ dnorm(0, 0.01)  # Coefficient of lambda on habitat
    phi ~ dunif(0,1)        # Probability of availability
    sigma ~ dunif(0.01,5)   # Distance function parameter
    
    # Detection probs for each distance interval and related things
    for(b in 1:nD){
    log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma) # half-normal 
    f[b] <- (2*midpt[b]*delta)/(B*B)    # radial density function
    cellprobs[b] <- g[b]*f[b]
    cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
    }
    cellprobs[nD+1]<- 1-sum(cellprobs[1:nD])
    
    for (s in 1:nsites) {
    for (k in 1:K) {
    pdet[s,k] <- sum(cellprobs[1:nD])   # Distance class probabilities
    pmarg[s,k] <- pdet[s,k]*phi         # Marginal probability
    
    # Model part 4: distance class frequencies
    y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[1:nD], nobs[s,k])
    # Model part 3: total number of detections:
    nobs[s,k] ~ dbin(pmarg[s,k], M[s])
    # nobs[s,k] ~ dbin(pdet[s,k], Navail[s,k]) # Alternative formulation
    # Model part 2: Availability. Not used in this model but simulated.
    Navail[s,k] ~ dbin(phi, M[s]) 
    }  # end k loop
    # Model part 1: Abundance model
    M[s] ~ dpois(lambda[s])    
    log(lambda[s]) <- beta0 + beta1*habitat[s]
    }  # End s loop
    
    # Derived quantities
    Mtot <- sum(M[])
    for(k in 1:K){ 
    Ntot[k]<- sum(Navail[,k])
    }
    } # End model
    ",file="model.txt")


# Assemble the initial values and parameters to save for JAGS
Navail.st <- apply(y3d, c(1,3),sum)
Mst <- apply(Navail.st, c( 1), max)  +2
inits <- function(){
  list(M=Mst, sigma = 1.0, phi=0.9, beta0=log(2), beta1=0.5)
}
params <- c("sigma", "phi", "beta0", "mean.lam", "beta1", "Mtot", "Ntot")

# MCMC settings
ni <- 60000   ;   nb <- 10000   ;   nt <- 5   ;   nc <- 3

# Run WinBUGS or JAGS
library("R2WinBUGS")
library("jagsUI")  # JAGS works but WinBUGS does not!
# bd <- "c:/WinBUGS14/"
# out1 <- bugs(data, inits, parameters, "model.txt", n.thin=nthin, 
#       n.chains=nc, n.burnin=nb,n.iter=ni,debug=TRUE, bugs.dir = bd)
# We get this error: vector valued relation y3d must involve consecutive 
# elements of variable

# Run JAGS: This fails quite often ('invalid parent node'), just keep trying
outTE1 <- jags(data, inits, params, "model.txt", n.thin=nt,n.chains=nc,
               n.burnin=nb,n.iter=ni, parallel = TRUE)
traceplot(outTE1)   ;    print(outTE1, 3)            # ART 4 min


# Put true values into a vector
truth <- c(tmp$parms[c(1:3,5)], Mtot = sum(tmp$M[,1]), 
           Ntot = (apply(tmp$Na.real[,,1],2,sum))) 

# Get posterior means and 2.5% and 97.5% percentiles (95% CRI)
post <- outTE1$summary[c("mean.lam", "beta1", "sigma", "phi", "Mtot", "Ntot[1]", "Ntot[2]", "Ntot[3]" ,"Ntot[4]", "Ntot[5]", "Ntot[6]", "Ntot[7]"), c(1,3,7)]

# Table compares truth with posterior mean and 95% CRI from JAGS
cbind(truth, posterior = round(post, 3))


# 9.5.4.2 Bayesian analysis of the Wagtail data
# ------------------------------------------------------------------------
y3d <- array(NA,dim=c(nrow(Y), 6, 4) )          # Create 3d array
y3d[,,1] <- Y[,1:6]  ;  y3d[,,2] <- Y[,7:12]    # Fill the array
y3d[,,3] <- Y[,13:18]  ;  y3d[,,4] <- Y[,19:24]

K <- 4                          # Number of primary occasions
nsites <- nrow(Y)               # Number of sites
nD <- 6                         # Number of distance classes
midpt <- seq(25,275,50)         # Class midpoint distance
delta <- 50                     # Class width
B <- 300                        # Maximum distance
nobs <- apply(y3d, c(1,3), sum) # Total detections per site and occasion

# Bundle and summarize data set
area <- pi*(300^2)/10000
str(data <- list(y3d=y3d, nsites=nsites, K=K, nD=nD, midpt=midpt, delta=delta, 
                 B=B, nobs=nobs, area=area))

# Write out the BUGS model file
cat("
    model {
    
    # Priors
    # Abundance parameters
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    
    # Availability parameter
    phi ~ dunif(0,1)
    
    # Detection parameter
    sigma ~ dunif(0,500)  
    
    # Multinomial cell probabilities
    for(b in 1:nD){
    log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma)  # Half-normal model
    f[b] <- (2*midpt[b]*delta)/(B*B) # Scaled radial density function
    cellprobs[b] <- g[b]*f[b]
    cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
    }
    cellprobs[nD+1] <- 1-sum(cellprobs[1:nD])
    
    for (s in 1:nsites) {
    for (k in 1:K) {
    # Conditional 4-part version of the model
    pdet[s,k] <- sum(cellprobs[1:nD])
    pmarg[s,k] <- pdet[s,k]*phi
    y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[1:nD], nobs[s,k]) # Part 4: distance  
    nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: number of detected individuals
    Navail[s,k] ~ dbin(phi,M[s])        # Part 2: Number of available individuals
    }  # end k loop
    
    M[s] ~ dpois(lambda[s])    #  Part 1: Abundance model
    log(lambda[s]) <- beta0    #  Habitat variables would go here
    }  # end s loop
    
    # Derived quantities
    for(k in 1:K){
    Davail[k] <- phi*exp(beta0)/area
    }
    Mtotal <- sum(M[])
    Dtotal<- exp(beta0)/area
    } # end model
    ",fill=TRUE,file="wagtail.txt")

# Inits
Navail.st <- apply(y3d, c(1,3),sum)  
Mst <- apply(Navail.st, c( 1), max,na.rm=TRUE) + 2
inits <- function() list(M=Mst, sigma = 100.0)

# Parameters to save
params <- c("sigma", "phi", "beta0", "beta1", "Mtotal", "Davail", "Dtotal")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3

# Run JAGS (ART 3 min)
library("jagsUI")
wag1 <- jags(data, inits, params, "wagtail.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, parallel = TRUE)
par(mfrow = c(3,3))    ;   traceplot(wag1)
summary(wag1)

exp(-2.06)


# Bundle and summmarize data set for BUGS
rep <- matrix(as.numeric(rep), ncol=4)
area <- pi*(300^2)/10000
str(data <- list(y3d=y3d, nsites=nsites, K=K, nD=nD, midpt = midpt, delta=delta, B=B, nobs=nobs, POTATO=POTATO, GRASS=GRASS, LSCALE=LSCALE, rep=rep, DATE=DATE,
                 HOUR=HOUR, area=area))

# Define model in BUGS
cat("
    model {
    
    # Priors
    # Abundance parameters
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    beta2 ~ dnorm(0, 0.01)
    beta3 ~ dnorm(0, 0.01)
    
    # Availability parameters
    phi0 ~ dunif(0,1)
    logit.phi0 <- log(phi0/(1-phi0))
    for(k in 1:4){
    gamma1[k] ~ dunif(0, 1) # Availability effects of surveys 1 - 4
    logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k])) 
    }
    gamma2 ~ dnorm(0, 0.01)  
    gamma3 ~ dnorm(0, 0.01)
    gamma4 ~ dnorm(0, 0.01)
    
    # Detection parameters
    sigma0 ~ dunif(0.1,500)   # Intercept  
    alpha2 ~ dnorm(0, 0.01)   # effect of DATE (linear)
    alpha3 ~ dnorm(0, 0.01)   # effect of DATE (squared)
    alpha4 ~ dnorm(0, 0.01)   # effect of HOUR
    theta ~ dgamma(0.1, 0.1)
    r ~ dunif(0, 10)
    
    for (s in 1:nsites) {
    for (k in 1:K) {
    # Availability parameter
    logit.phi[s,k] <- logit.gamma1[k] + gamma2*POTATO[s] + gamma3*GRASS[s] + gamma4*LSCALE[s]
    phi[s,k] <- exp(logit.phi[s,k])/(1+ exp(logit.phi[s,k]))
    # Distance sampling parameter
    log(sigma[s,k]) <- log(sigma0) + alpha2*DATE[s,k] + alpha3*pow(HOUR[s,k],2) + alpha4*HOUR[s,k]
    # Multinomial cell probability construction
    for(b in 1:nD){
    #log(g[s,b,k]) <- -midpt[b]*midpt[b]/(2*sigma[s,k]*sigma[s,k]) # half-normal
    cloglog(g[s,b,k]) <- theta*log(sigma[s,k])  - theta*log(midpt[b])  # hazard
    f[s,b,k] <- (2*midpt[b]*delta)/(B*B)
    cellprobs[s,b,k] <- g[s,b,k]*f[s,b,k]
    cellprobs.cond[s,b,k] <- cellprobs[s,b,k]/sum(cellprobs[s,1:nD,k])
    }
    cellprobs[s,nD+1,k]<- 1-sum(cellprobs[s,1:nD,k])
    
    #  Conditional 4-part hierarchical model
    pdet[s,k] <- sum(cellprobs[s,1:nD,k])
    pmarg[s,k] <- pdet[s,k]*phi[s,k]
    y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k]) # Part 4
    nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: Number of detected individuals
    Navail[s,k] ~ dbin(phi[s,k],M[s])   # Part 2: Number of available individuals
    } # end k loop
    
    M[s] ~ dnegbin(prob[s], r)
    prob[s] <- r/(r+lambda[s])
    # M[s] ~ dpois(lambda[s])         # Part 1: Abundance model
    log(lambda[s]) <- beta0 + beta1*POTATO[s] + beta2*GRASS[s] + beta3*LSCALE[s]
    }  # end s loop
    
    # Derived quantities
    for(k in 1:K){
    Davail[k] <- mean(phi[,k])*exp(beta0)/area
    }
    Mtotal <- sum(M[])
    Dtotal <- exp(beta0)/area
    } # End model
    ",fill=TRUE,file="wagtail2.txt")

# Inits
Navail.st <- apply(y3d, c(1,3),sum)
Mst <- apply(Navail.st, c( 1), max,na.rm=TRUE) + 2
inits <- function() list(M=Mst, sigma0 = 100, alpha2=0, alpha3=0, alpha4=0, gamma2=0, gamma3=0, gamma4=0, beta1=0,beta2=0,beta3=0, r = 1)

# Parameters to save
params <- c("r","sigma0", "beta0", "beta1", "beta2", "beta3", "Mtotal", "alpha2", "alpha3",  "alpha4", "theta", "Dtotal", "Davail", "phi0", "gamma1", "gamma2", "gamma3", "gamma4" ,"logit.gamma1")

# MCMC settings
ni <- 32000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 5

# Run JAGS (ART 79 min), check convergence and summarize posteriors
wag2 <- jags(data, inits, params, "wagtail2.txt", n.thin=nt,n.chains=nc, 
             n.burnin=nb,n.iter=ni, parallel = TRUE)

# Compare posterior means to MLEs obtained from unmarked
mle <- coef(fm5)
mle[12] <- exp(mle[12]) # convert to distance units
mle[16] <- exp(mle[16]) # back-transform the hazard parameter
mle[17] <- exp(mle[17]) # back-transform the NB dispersion parameter
bayes <- wag2$summary[,1]
bayes <- c(bayes[3:6],bayes[c(25:28,22:24)],bayes[c(2,8:10)],bayes[11],bayes[1])
bayes[1] <- log(exp(bayes[1])/area)# Convert from N per site to log(density) per ha

round( cbind(mle,bayes), 3)


# 9.5.4.3 Robust Design: Replicates within and among years
# ------------------------------------------------------------------------



# 9.6 Open HDS models: Implicit Dynamics
# ------------------------------------------------------------------------


# Obtain a data set
set.seed(1236)
str(tmp <- simHDSopen("point", nreps=7, nyears=5, nsites=100, beta.trend=0.2) )
attach(tmp)

apply(tmp$M.true,2,sum)  # True population size per year

# Define distance class information
delta <- 0.5
nD <- B%/%delta                 # Number of distance classes
midpt <- seq(delta/2, B, delta) # Mid-point of distance intervals

# Create the 4-d array
y4d <- array(0, dim=c(nsites, nD, K, nyears))
for(yr in 1:nyears){
  for(rep in 1:K){
    data <- tmp$data[[yr]][[rep]]
    site <- data[,1]
    dclass <- data[,"d"]%/%delta + 1
    ndclass <- B%/%delta
    dclass <- factor(dclass, levels=  1:ndclass)
    y4d[1:nsites,1:nD,rep,yr] <- table(site, dclass)
  }
}


# Bundle and summarize the data set
nobs <- apply(y4d, c(1,3,4), sum)  # Total detections per site and occasion
str( data <- list(y4d=y4d, nsites=nsites, K=K, nD=nD, midpt=midpt, delta=delta, habitat=habitat, B=B, nobs = nobs, T=tmp$nyears) )

# Define model in BUGS
cat("
    model {
    
    # Prior distributions
    beta0 ~ dnorm(0, 0.01)  # Intercept for log(lambda)
    mean.lam <- exp(beta0)
    beta1 ~ dnorm(0, 0.01)  # Coefficient on habitat
    phi ~ dunif(0,1)        # Probability of availability
    sigma ~ dunif(0,5)      # Detection function parameter
    beta.trend ~ dnorm(0, 0.01)
    
    # Construct the multinomial cell probabilities
    for(b in 1:nD){
    log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma) # half-normal 
    f[b] <- (2*midpt[b]*delta)/(B*B)                # radial density function
    cellprobs[b] <- g[b]*f[b]
    cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
    }
    cellprobs[nD+1] <- 1-sum(cellprobs[1:nD])
    for (s in 1:nsites) {
    for (k in 1:K) {
    pdet[s,k] <- sum(cellprobs[1:nD]) # Distance class probabilities
    pmarg[s,k] <- pdet[s,k]*phi       # Marginal probability
    }
    }
    
    for(t in 1:T){                        # Years
    for (s in 1:nsites) {               # Sites
    for (k in 1:K) {                  # Replicates
    # Model part 4: distance class frequencies
    y4d[s,1:nD,k,t] ~ dmulti(cellprobs.cond[1:nD], nobs[s,k,t])
    # Model part 3: total number of detections:
    nobs[s,k,t] ~ dbin(pmarg[s,k], M[s,t])  
    # Model part 2: Availability. Not used in this model but simulated.
    Navail[s,k,t] ~ dbin(phi, M[s,t]) 
    }  # end k loop
    # Model part 1: Abundance model
    M[s,t] ~ dpois(lambda[s,t])
    log(lambda[s,t]) <- beta0 + beta1*habitat[s] + beta.trend*(t-2.5)
    }  # end s loop
    } # end t loop
    
    # Derived quantities
    for(t in 1:T){
    Mtot[t] <- sum(M[,t])
    for(k in 1:K){ 
    Ntot[k,t] <- sum(Navail[,k,t])
    }
    }
    } # End model
    ",file="tempemig4d.txt")

# Inits and parameters to save
Navail.st <- apply(y4d, c(1,3,4),sum)
Mst <- apply(Navail.st, c( 1,3), max) +2
inits <- function(){
  list(M=Mst, Navail = Navail.st, sigma = 1.0, phi=.9,beta0=log(2),beta1=.5)
}
params <- c("sigma", "phi", "beta0", "mean.lam", "beta.trend",
            "beta1", "Mtot", "Ntot")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 5   ;   nc <- 3

# Run JAGS (ART 9 min), look at trace plots and summarize
outRD <- jags(data, inits, params, "tempemig4d.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, parallel = FALSE)
par(mfrow = c(3,3))   ;   traceplot(outRD)
summary(outRD)



# 9.7 Open HDS models: modelling population dynamics
# ------------------------------------------------------------------------


# 9.7.1 Simulating the ISSJ data over multiple years
# ------------------------------------------------------------------------
# We load the ISSJ data analyzed in chapter 8, package into an unmarked frame  
library(unmarked)
library(rjags)
data(issj)
covs <- issj[,c("elevation","forest","chaparral")]
area <- pi*300^2 / 100^2             # Area in ha
jayumf <- unmarkedFrameGDS(y=as.matrix(issj[,1:3]),
                           siteCovs=data.frame(covs, area), numPrimary=1,
                           dist.breaks=c(0, 100, 200, 300),
                           unitsIn="m", survey="point")
sc <- siteCovs(jayumf)
sc.s <- scale(sc)
sc.s[,"area"] <- pi*300^2 / 10000  # Don't standardize area
covs<- siteCovs(jayumf) <- sc.s
summary(jayumf)

# Fit the model using gdistsamp and look at the fit summary
(nb.C2E.C <- gdistsamp( ~chaparral + I(chaparral^2) + elevation , ~1, ~chaparral,
                        data =jayumf, output="abund", mixture="NB", K = 150))

# Get coefficient estimates to be used in data simulation
beta <- coef(nb.C2E.C)
betaFall <- beta[c("lambda(Int)", "lambda(chaparral)",
                   "lambda(elevation)", "lambda(I(chaparral^2))")]

# Predict expected abundance per point count on log-scale for simulation
Xmat <- cbind(rep(1,307),covs[,3],covs[,3]^2,covs[,1]) # Order: chap, chap^2, elev
loglam <- Xmat%*%(betaFall)
lamnew <- exp(loglam)

# Parameters of the detection function
dparm <- beta[c("p(Int)", "p(chaparral)")]
sigma <- exp(Xmat[,c(1, 2)]%*%dparm)
J <- nsites <- 307 # number of sampling points

# Number of years
nyrs <- 6

# Set dynamics parameters to achieve a target growth rate of 0.95
phi <- 0.6       # Survival probability
gamma <- 0.35    # Recruitment rate

# Distance category info
db <- c(0,50, 100, 150, 200, 250, 300)
midpt <- c(25, 75, 125, 175, 225, 275)
nD <- length(midpt)
delta <- 50       # Distance interval width
B <- 300

# Simulate an ISSJ data set and harvest the data objects
set.seed(2015)
dat <- issj.sim(B=300, db = db, lam=lamnew, sigma=sigma, phi=phi, gamma=gamma, npoints=nsites, nyrs=nyrs)

y <- dat$y
dclass <- dat$dclass
site <- dat$site

# Bundle and summarize the data set
str(data1<-list(nsites=nsites, chap=as.vector(covs[,"chaparral"])[dat$cell], 
                chap2=as.vector(covs[,"chaparral"]^2)[dat$cell],
                elev=as.vector(covs[,"elevation"])[dat$cell], T=nyrs, nD=nD, midpt=midpt, 
                B=B, delta=delta, y=y, dclass=dclass, site=site, nind=sum(y)) )



# 9.7.2 Fitting a slurry of open population models  
# ------------------------------------------------------------------------


# 9.7.2.1 The independence model
# ------------------------------------------------------------------------
# Write out the BUGS model file
cat("
    model{
    
    # Prior distributions
    # Regression parameters
    alpha0 ~ dunif(0,20)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-20,20)
    beta1 ~ dunif(-20,20)
    beta2 ~ dunif(-20,20)
    beta3 ~ dunif(-20,20)
    beta4 ~ dunif(-20,20) # Population trend parameter
    r ~ dunif(0,5)        # NegBin dispersion parameter
    rout <- log(r)
    
    # 'Likelihood'
    for (s in 1:nsites){
    # Linear model for detection function scale
    log(sigma[s]) <- alpha0+alpha1*chap[s]
    # Compute detection probability
    for(k in 1:nD){
    pi[k,s] <- (2*midpt[k]*delta )/(B*B) 
    log(p[k,s]) <- -midpt[k]*midpt[k]/(2*sigma[s]*sigma[s])
    f[k,s] <- p[k,s]*pi[k,s]
    fc[k,s] <- f[k,s]/pcap[s]
    fct[k,s] <- fc[k,s]/sum(fc[1:nD,s])
    }
    pcap[s]<-sum(f[1:nD,s])  # Overall detection probability
    
    # Process model
    for (t in 1:T){
    log(lambda[s,t]) <- beta0 + beta1*chap[s] + beta2*chap2[s] + beta3*elev[s] + beta4*(t - t/2)  # Note trend parameter here
    y[s,t] ~ dbin(pcap[s], N[s,t])
    N[s,t] ~ dnegbin(prob[s,t], r)
    prob[s,t] <- r/(r+lambda[s,t])
    } # End loop over years
    } # End loop over sites
    
    # Distance sampling observation model for observed (binned) distance data
    for(i in 1:nind){
    dclass[i] ~ dcat(fct[1:nD,site[i]]) 
    }
    # Derived parameters
    for(t in 1:6){
    Ntot[t] <- sum(N[,t])
    D[t] <- Ntot[t] / (28.27*nsites)   # 300 m point = 28.27 ha
    }
    }
    ", file="Sollmann1.txt")

# Set up initial values, parameters vector and MCMC settings
Nst <- y+1  # this is for trend model
inits <- function(){list(N=Nst, beta0=runif(1), beta1=runif(1), beta2=runif(1), 
                         beta3=runif(1), beta4=runif(1), alpha0=runif(1,3,5), alpha1=runif(1), r = 1)}
params <-c('beta0', 'beta1', 'beta2', 'beta3', 'beta4', 'alpha0', 'alpha1', 'Ntot', 'D', 'r')
ni <- 22000   ;   nb <- 2000   ;   nt <- 1   ;   nc <- 3

# JAGS setting b/c otherwise JAGS cannot build a sampler, rec. by M. Plummer
set.factory("bugs::Conjugate", FALSE, type="sampler")

# Execute JAGS, look at convergence and summarize the results
library(jagsUI)
open1 <- jags (data1, inits, params, "Sollmann1.txt", n.thin=nt, n.chains=nc,
               n.burnin=nb, n.iter=ni)
par(mfrow = c(3,3))   ;   traceplot(open1)   ;   print(open1, 2)




# 9.7.2.2 The reduced-dynamics model
# ------------------------------------------------------------------------
# Write out the BUGS model file
cat("
    model{
    
    # Prior distributions
    # Regression parameters
    alpha0 ~ dunif(0,20)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-20,20)
    beta1 ~ dunif(-20,20)
    beta2 ~ dunif(-20,20)
    beta3 ~ dunif(-20,20)
    theta ~ dunif(0,5)
    # NegBin dispersion parameter
    r ~ dunif(0,5)
    rout <- log(r)
    
    # 'Likelihood'
    for (s in 1:nsites){
    # Linear model for detection function scale
    log(sigma[s]) <- alpha0+alpha1*chap[s]
    # Compute detection probability
    for(k in 1:nD){
    log(p[k,s]) <- -midpt[k]*midpt[k]/(2*sigma[s]*sigma[s])
    f[k,s] <- p[k,s]*pi[k,s]
    fc[k,s] <- f[k,s]/pcap[s]
    fct[k,s] <- fc[k,s]/sum(fc[1:nD,s])
    pi[k,s] <- (2*midpt[k]*delta )/(B*B)
    }
    pcap[s]<-sum(f[1:nD,s])  # Overall detection probability
    
    # Process model
    # Abundance model for Yr1 as in Sillett et al 2012
    log(lambda[s,1]) <- beta0 + beta1*chap[s] + beta2*chap2[s] + beta3*elev[s]
    y[s,1] ~ dbin(pcap[s], N[s,1])
    N[s,1] ~ dnegbin(prob[s,1], r)
    prob[s,1] <- r/(r+lambda[s,1])
    
    # Population dynamics model for subsequent years
    for (t in 2:T){
    N[s,t] ~ dpois(N[s, t-1] * theta)
    y[s,t] ~ dbin(pcap[s], N[s,t])
    }
    }
    # Distance sampling observation model for observed (binned) distance data
    for(i in 1:nind){
    dclass[i] ~ dcat(fct[1:nD,site[i]]) 
    }
    
    # Derived parameters
    for(t in 1:6){
    Ntot[t] <- sum(N[,t])
    D[t] <- Ntot[t] / (28.27*nsites)  # 300 m point = 28.27 ha
    }
    }
    ", file="Sollmann2.txt")

# Set up initial values, parameters vector and MCMC settings
Nst <- y+1 # this is for trend model
inits <- function(){list(N=Nst, beta0=runif(1), beta1=runif(1), beta2=runif(1), 
                         beta3=runif(1), alpha0=runif(1,3,5), alpha1=runif(1), theta=runif(1,0.6,0.99))} 
params <- c('beta0', 'beta1', 'beta2', 'beta3', 'alpha0', 'alpha1', 'theta', 
            'rout', 'Ntot', 'D', 'r')
ni <- 22000   ;   nb <- 2000   ;   nt <- 1   ;   nc <- 3

# Execute JAGS, look at convergence and summarize the results
open2 <- jags (data1, inits, params, "Sollmann2.txt", n.thin=nt, n.chains=nc, 
               n.burnin=nb, n.iter=ni)
par(mfrow = c(3,3))   ;   traceplot(open2)   ;   print(open2, 2)



# 9.7.2.3 The Glorious Integrated HDS/Dail-Madsen Model
# ------------------------------------------------------------------------
# Write out the BUGS model file
cat("
    model{
    # Prior distributions
    # Regression parameters
    alpha0 ~ dunif(0,20)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-20,20)
    beta1 ~ dunif(-20,20)
    beta2 ~ dunif(-20,20)
    beta3 ~ dunif(-20,20)
    
    # Priors for dynamics parameters: here they are constant across years
    # We could add covariate models for logit(phi) and log(gamma)
    phi ~ dunif(0,1)
    gamma ~ dunif(0,5)
    
    # NegBin dispersion parameter
    r ~ dunif(0,5)
    rout <- log(r)
    
    # 'Likelihood'
    for (s in 1:nsites){
    # Linear model for detection function scale
    log(sigma[s]) <- alpha0+alpha1*chap[s]
    
    # Compute detection probability
    for(k in 1:nD){
    log(p[k,s]) <- -midpt[k]*midpt[k]/(2*sigma[s]*sigma[s])
    f[k,s] <- p[k,s]*pi[k,s]
    fc[k,s] <- f[k,s]/pcap[s]
    fct[k,s] <- fc[k,s]/sum(fc[1:nD,s])
    pi[k,s] <- (2*midpt[k]*delta )/(B*B)
    }
    pcap[s]<-sum(f[1:nD,s])  # Overall detection probability
    
    # Process model
    # Abundance model for year 1
    log(lambda[s,1]) <- beta0 + beta1*chap[s] + beta2*chap2[s] + beta3*elev[s]
    y[s,1] ~ dbin(pcap[s], N[s,1])
    N[s,1] ~ dnegbin(prob[s,1], r)
    prob[s,1] <- r/(r+lambda[s,1])
    
    # Population dynamics model for subsequent years
    for (t in 2:T){                      # Loop over years
    S[s,t] ~ dbinom(phi, N[s, t-1])   # Survivors
    R[s,t] ~ dpois(gamma * N[s, t-1]) # Recruits
    N[s,t] <- S[s,t] + R[s,t]         # N = Survivors + Recruits
    y[s,t]~ dbin(pcap[s],N[s,t])       # Measurement error
    }
    }
    
    # Distance sampling observation model for observed (binned) distance data
    for(i in 1:nind){
    dclass[i] ~ dcat(fct[1:nD,site[i]])
    }
    
    # Derived parameters
    for(t in 1:6){
    Ntot[t] <- sum(N[,t])
    D[t] <- Ntot[t] / (28.27*nsites)     # 300 m point = 28.27 ha
    } 
    }
    ", file="Sollmann3.txt")


# Set up some sensible starting values for S and R
yin <- y+1
yin[,2:6] <- NA
Sin <- Rin <- matrix(NA, nrow=nsites, ncol=nyrs)
y1 <- y + 1
for(s in 1:nsites){
  for (t in 2:6){
    Sin[s,t] <- rbinom(1,y1[s,t-1], phi )
    Rin[s,t] <- ifelse((y1[s,t]-Sin[s,t])>0, y1[s,t]-Sin[s,t], 0)
  }
}

# Set up initial values, parameters vector and MCMC settings
inits <-function(){list(N=yin, beta0=runif(1), beta1=runif(1), beta2=runif(1),
                        beta3=runif(1), alpha0=runif(1,3,5), alpha1=runif(1), phi=0.6, gamma=0.3, 
                        R=Rin, S=Sin) }
params <- c('beta0', 'beta1', 'beta2', 'beta3', 'alpha0', 'alpha1', 'phi', 
            'gamma', 'Ntot', 'D', 'r')
ni <- 152000   ;   nb <- 2000   ;   nt <- 10   ;   nc <- 5

# Run JAGS, look at convergence and summarize the results
library(jagsUI)
open3  <- jags (data1, inits, params, "Sollmann3.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, parallel=TRUE)
par(mfrow = c(3,3))   ;   traceplot(open3)   ;   print(open3, 2)


# Compare inferences in graph .... (Fig. 9-6)
plot(apply(dat$N,2,sum),ylim=c(600,1300),xlab="Year",ylab="Population size (307 sample units)")
lines(apply(open1$sims.list$Ntot,2,mean), lty=1, col="blue", lwd=2)
lines(apply(open2$sims.list$Ntot,2,mean), lty=1, col="red", lwd=2)
lines(apply(open3$sims.list$Ntot,2,mean), lty=1, col="green", lwd=2)
open1.95cri<- apply(open1$sims.list$N,2,function(x) quantile(x, c(0.025,0.975)))
open2.95cri<- apply(open2$sims.list$N,2,function(x) quantile(x, c(0.025,0.975)))
open3.95cri<- apply(open3$sims.list$N,2,function(x) quantile(x, c(0.025,0.975)))
legend(1,750,legend=c("Independence","Reduced dynamics","Full dynamics"), lty=1, col=c("blue","red","green"))

matlines(1:6, t(open1.95cri), type="l", lty=2, lwd=2, col="blue")
matlines(1:6, t(open2.95cri), type="l", lty=2, lwd=2, col="red")
matlines(1:6, t(open3.95cri), type="l", lty=2, lwd=2, col="green")

# .... and table
parms <- c(betaFall, dparm)
round(post <- cbind(parms, Independent=open1$summary[c(1:4,6,7),1], 
                    Partial=open2$summary[1:6,1], Full=open3$summary[1:6,1]), 3)

parms Independent Partial   Full
lambda(Int)             0.827       0.990   0.971  0.949
lambda(chaparral)       1.432       1.755   1.783  1.806
lambda(elevation)      -0.227      -0.509  -0.519 -0.523
lambda(I(chaparral^2)) -0.376      -0.343  -0.317 -0.323
p(Int)                  4.679       4.682   4.672  4.680
p(chaparral)           -0.199      -0.194  -0.187 -0.192



# 9.7.2.4 Summary remarks on modelling populations over time
# ------------------------------------------------------------------------



# 9.8 Spatial Distance Sampling: Modelling within-unit variation in density
# ------------------------------------------------------------------------


# 9.8.1 Distance sampling with location of encounter
# ------------------------------------------------------------------------
# Simulate a data set and harvest the output
set.seed(1234)
str(tmp <- sim.pdata(N=200,sigma=1,keep.all=FALSE,B=3))

List of 8
$ N     : num 200
$ sigma : num 1
$ B     : num 3
$ u1    : num [1:37] 3.73 3.66 3.27 1.72 3.32 ...
$ u2    : num [1:37] 3.17 1.9 1.68 3.83 2.54 ...
$ d     : num [1:37] 0.753 1.276 1.345 1.529 0.557 ...
$ y     : int [1:200] 0 1 1 0 0 0 0 0 0 0 ...
$ N.real: int 152

# Harvest some data objects
B <- tmp$B
d <- tmp$d
u1 <- tmp$u1
u2 <- tmp$u2
nind <- length(d)

# Data augmentation
M <- 400                       # max of 400 individuals
nz <- M-nind                   # augment by nz individuals
y <- c(rep(1,nind), rep(0,nz)) # augmented data augmentation variable
u <- cbind(u1,u2)               # Locations of individuals detected
u <- rbind(u, matrix(NA, nrow=nz, ncol=2))


# Bundle and summarize the data set
str(data <- list (B=B, nind=nind, u=u, y=y, nz=nz))


# Write out the BUGS model file
cat("
    model{ 
    
    # Priors
    sigma ~ dunif(0,10)
    psi ~ dunif(0,1)
    
    # Categorical observation model
    for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)
    u[i,1] ~ dunif(0, 2*B)  # Here is the uniformity assumption made explicit
    u[i,2] ~ dunif(0, 2*B)
    # Compute distance as a derived quantity
    d[i] <- pow( pow( u[i,1]-B,2) + pow(u[i,2]-B,2), 0.5) # Pythagoras
    p[i] <- exp(-d[i]*d[i] / (2*sigma*sigma))
    mu[i] <- p[i] * z[i]
    y[i] ~ dbern(mu[i])
    }
    # Other derived quantity
    N <- sum(z[])
    D <- N / (B*B)
    }
    ",fill=TRUE,file="model1.txt")


# Load libraries and specify MCMC settings
library("R2WinBUGS")   ;   library(jagsUI)
ni <- 22000   ;   nb <- 2000   ;   nthin <- 2   ;   nc <- 3

# Inits and parameters
inits <- function(){
  list(sigma=runif(1,1,10), psi=runif(1),z = c(rep(1,nind),rep(0,nz)) ) }
params <- c("sigma", "N", "psi")

# Execute jags and summarize the posterior distributions
out1 <- jags (data, inits, parameters, "model1.txt", n.thin=nthin, 
              n.chains=nc, n.burnin=nb,n.iter=ni, parallel = FALSE)
par(mfrow = c(2,2))   ;   traceplot(out1)
print(out1, 2)

# 9.8.2 The line transect case
# ------------------------------------------------------------------------


# 9.8.3 Modelling spatial covariates
# ------------------------------------------------------------------------
# Simulator function for spatial distance sampling data
sim.spatialDS <- 
  function(N=1000, beta = 1, sigma=1, keep.all=FALSE, B=B, model="halfnorm"){
    # Function simulates coordinates of individuals on a square
    # Square is [0,2B] x [0,2B], with a count location on the point (B, B)
    #   N: total population size in the square
    #   beta: coefficient of SOEMTHING on spatial covariate x
    #   sigma: scale of half-normal detection function
    #   B: circle radius
    #   keep.all: return the data for y=0 individuals or not
    library(raster)      # Load required packages
    library(plotrix)
    
    # Create coordinates for 30 x 30 grid
    delta <- (2*B-0)/30                # '2D bin width'
    grx <- seq(delta/2, 2*B - delta/2, delta) # mid-point coordinates
    gr <- expand.grid(grx,grx)         # Create grid coordinates
    
    # Create spatially correlated covariate x and plot it
    V <- exp(-e2dist(gr,gr)/1)
    x <- t(chol(V))%*%rnorm(900)
    par(mar=c(3,3,3,6))
    image(rasterFromXYZ(cbind(gr,x)), col=topo.colors(10))
    draw.circle(3, 3, B)
    points(3, 3, pch="+", cex=3)
    image.scale(x, col=topo.colors(10))
    
    # Simulate point locations as function of habitat covariate x
    probs <- exp(beta*x)/sum(exp(beta*x)) # probability of point in pixel (sum = 1)
    pixel.id <- sample(1:900, N, replace=TRUE, prob=probs)
    # could simulate randomly within the pixel but it won't matter so place centrally
    u1 <- gr[pixel.id,1]
    u2 <- gr[pixel.id,2]
    points(u1, u2, pch=20, col='black', cex = 0.8)  # plot points
    title("This is so cool !")         # express your appreciation of all this
    
    d <- sqrt((u1 - B)^2 + (u2-B)^2)   # distance to center point of square
    #plot(u1, u2, pch = 1, main = "Point transect")
    N.real <- sum(d<= B)               # Population size inside of count circle
    
    # Can only count individuals in the circle, so set to zero detection probability of individuals in the corners (thereby truncating them)
    # p <- ifelse(d< B, 1, 0) * exp(-d*d/(2*(sigma^2)))
    # We do away with the circle constraint here.   
    if(model=="hazard")
      p <- 1-exp(-exp(-d*d/(2*sigma*sigma)))
    if(model=="halfnorm")
      p <- exp(-d*d/(2*sigma*sigma))
    # Now we decide whether each individual is detected or not
    y <- rbinom(N, 1, p)                                           # detected or not
    points(u1[d<= B], u2[d<= B], pch = 16, col = "black", cex = 1) # not detected
    points(u1[y==1], u2[y==1], pch = 16, col = "red", cex = 1)     # detected
    
    # Put all of the data in a matrix
    if(!keep.all){
      u1 <- u1[y==1]
      u2 <- u2[y==1]
      d <- d[y==1]   }
    # Output
    return(list(model=model, N=N, beta=beta, B=B, u1=u1, u2=u2, d=d, y=y, N.real=N.real, Habitat=x, grid=gr))
  }

# Generate one data set and harvest the output
set.seed(1234)
str(tmp <- sim.spatialDS(N=200, beta=1, sigma=1.5, keep.all=FALSE, B=3)) # Fig. 9-7


# Harvest data
B <- tmp$B
d <- tmp$d
u1 <- tmp$u1
u2 <- tmp$u2
Habitat <- as.vector(tmp$Habitat)
Habitat <- Habitat - mean(Habitat)
Habgrid <- tmp$grid
nind <- length(d)
G <- nrow(Habgrid)

# Do data augmentation, including for pixel ID
M <- 400
nz <- M-nind
pixel <- rep(NA, M)   # We use discrete "pixel ID" here instead of "s"
y <- c(rep(1,nind), rep(0,nz))

# Pick some starting values and figure out the pixel of each observation
s <- cbind(u1,u2)
s <- rbind(s, matrix(NA,nrow=nz,ncol=2))
D <- e2dist(s[1:nind,], Habgrid)
for(i in 1:nind){
  pixel[i] <- (1:ncol(D))[D[i,]==min(D[i,])]
}

# Bundle and summarize the data for BUGS
str(data <- list (B=B, nind=nind, y=y, nz=nz, Habitat=Habitat, Habgrid=Habgrid, G=G, pixel=pixel))


# Write BUGS model
cat("
    model{ 
    
    # Prior distributions
    sigma ~ dunif(0,10)
    psi ~ dunif(0,1)
    beta ~ dnorm(0,0.01)
    
    for(g in 1:G){   # g is the pixel index, there are G total pixels
    probs.num[g] <- exp(beta*Habitat[g])
    probs[g] <- probs.num[g]/sum(probs.num[])
    }
    
    # Models for DA variables and location (pixel)
    for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)
    pixel[i] ~ dcat(probs[])
    s[i,1:2] <- Habgrid[pixel[i],]   # location = derived quantity  
    # compute distance = derived quantity
    d[i] <- pow(   pow( s[i,1]-B,2) + pow(s[i,2]-B,2), 0.5)
    p[i] <- exp(-d[i]*d[i]/(2*sigma*sigma))  # Half-normal detetion function
    mu[i]<- p[i]*z[i]
    y[i] ~ dbern(mu[i])                      # Observation model
    }
    # Derived parameters
    N <- sum(z[])                      # N is a derived parameter
    D <- N/9                           # area = 9 ha 
    }
    ",fill=TRUE, file="spatialDS.txt")


# Load libraries and specify MCMC settings
library("R2WinBUGS")
library(jagsUI)
ni <- 12000   ;   nb <- 2000   ;   nthin <- 2   ;   nc <- 3

# Create inits and define parameters to monitor
inits <- function(){  list (sigma=runif(1,1,10),psi=runif(1), 
                            z = c(rep(1,nind), rep(0,nz)) ) }
params <- c("sigma", "N", "psi", "beta", "D")

# Run JAGS, check convergence and summarize posteriors
out2 <- jags (data, inits, params, "spatialDS.txt", n.thin=nthin,
              n.chains=nc, n.burnin=nb, n.iter=ni, parallel = TRUE)
par(mfrow = c(2,3)   ;   traceplot(out2)
    print(out2, 2)
    
    
    # Add pixel in order to make a density map
    params <- c("sigma", "N", "psi", "beta", "D", "pixel")
    
    # Run JAGS, check convergence and summarize posteriors
    out2 <- jags (data, inits, params, "spatialDS.txt", n.thin=nthin,
                  n.chains=nc, n.burnin=nb, n.iter=ni, parallel = FALSE)
    
    # Plot density maps
    library(raster) 
    par(mfrow=c(1,2))
    pixel <- out2$sims.list$pixel
    post <- table(pixel)/nrow(pixel)   # Average number of locations in each pixel
    prior.mean <- mean(out2$sims.list$beta)*as.vector(Habitat)
    prior.mean <- mean(out2$sims.list$psi)*M*exp(prior.mean)/sum(exp(prior.mean))
    plot(rast.data <- rasterFromXYZ(cbind(Habgrid,prior.mean)), axes=FALSE, 
         col=topo.colors(10) )
    title("Prior mean density (estimated)")
    plot(rast.post <- rasterFromXYZ(cbind(Habgrid,as.vector(post))),axes=FALSE, 
         col=topo.colors(10) )
    title("Posterior mean density")
    
    
    # 9.8.4 Spatial HDS models in unmarked using the pcount function
    # ------------------------------------------------------------------------
    # Simulate a data set, N = 600 for the population size
    set.seed(1234)
    tmp <-sim.spatialDS(N=600, sigma=1.5, keep.all=FALSE, B=3, model= "hazard")
    
    # Harvest stuff
    B <- tmp$B
    d <- tmp$d
    u1 <- tmp$u1
    u2 <- tmp$u2
    Habitat <- as.vector(tmp$Habitat)
    Habitat <- Habitat - mean(Habitat)
    Habgrid <- tmp$grid
    nind <- length(d)
    G <- nrow(Habgrid)
    
    # Find which pixel each observation belongs to
    s <- cbind(u1,u2)
    D <- e2dist(s[1:nind,], Habgrid)
    pixel <-rep(NA,nind)
    for(i in 1:nind){
      pixel[i] <- (1:ncol(D))[D[i,]==min(D[i,])]
    }
    
    # Create a vector of counts in each pixel and pad it with zeros
    pixel.count <- rep(0, G)
    names(pixel.count) <- 1:G
    pixel.count[names(table(pixel))] <- table(pixel)
    # Create a covariate: distance between observer and pixel center
    dist <- sqrt( (Habgrid[,1]-3)^2 + (Habgrid[,2]-3)^2  )
    # Construct an unmarkedFrame
    umf <- unmarkedFramePCount(y=matrix(pixel.count,ncol=1), 
                               siteCovs=data.frame(dist=dist,Habitat=Habitat))
    summary(umf)
    
    
    # Fit an N-mixture model with no intercept and distance squared using
    #   the hacked function pcount.hds
    
    (fm1 <- pcount.spHDS(~ -1 + I(dist^2) ~ Habitat, umf, K=20))
    
    lam <- exp( coef(fm1)[1] + coef(fm1)[2]*Habitat )
    pred <- predict(fm1, type='state')
    sum(lam)
    sum(pred[,1])  # Same
    
    
    # 9.8.5 Hierarchical spatial distance sampling
    # ------------------------------------------------------------------------
    sim.spatialHDS(lam0 = 4, sigma = 1.5, B = 3, nsites = 100)
    # lam0 = expected population size per site
    # nsites = number of point count locations
    # B = count radius. Function simulates coordinates of individuals on a square
    #       [0,2*B] x[0,2*B], with a count location on the point (B,B)
    # sigma = scale of half-normal detection function
    
    library(raster)
    
    # Simulate a data set and harvest the output
    set.seed(1234)
    str(tmp <-sim.spatialHDS(lam0 = 3, sigma = 1.5, B = 3, nsites = 100))
    
    # Process the simulated data set
    data <- tmp$data
    # To make it a ërealí data set:
    data <- data[!is.na(data[,2]),] # Get rid of the 0 sites 
    data <- data[data[,"y"]==1,]    # Only keep detected individuals
    
    # Now zero-pad the observed counts
    nsites <- tmp$nsites
    nobs <- rep(0, nsites)
    names(nobs) <- 1:nsites
    nobs[names(table(data[,1]))] <- table(data[,1])
    
    # Extract data elements that we need
    site <- data[,"site"]
    s <- data[,c("u1","u2")]
    B <- tmp$B
    Habitat <- (tmp$Habitat) # Raster values
    Habgrid <- tmp$grid   # Raster coordinates
    nind <- nrow(data)
    G <- nrow(Habgrid)
    
    # We have to convert observed locations to pixels
    pixel <- rep(NA,nrow(data))
    D <- e2dist(s[1:nind,], Habgrid)
    for(i in 1:nind){
      pixel[i] <- (1:ncol(D))[D[i,]==min(D[i,])]
    }
    
    # Do data augmentation of data from each site ìS-foldî DA. Three objects need
    # to have DA applied to them: Ymat, umat, pixmat
    Msite <- 2*max(nobs)  # Perhaps use a larger value 
    Ymat <- matrix(0,nrow=Msite,ncol=nsites)
    umat <- array(NA, dim=c(Msite,nsites,2))
    pixmat <- matrix(NA,nrow=Msite,ncol=nsites)
    for(i in 1:nsites){
      if(nobs[i]==0) next
      Ymat[1:nobs[i],i]<- data[data[,1]==i,"y"]
      umat[1:nobs[i],i,1:2]<- data[data[,1]==i,c("u1","u2")]
      pixmat[1:nobs[i],i]<- pixel[data[,1]==i]
    }
    
    # Bundle the data for BUGS
    str(data <- list (y=Ymat, pixel=pixmat, Habitat=Habitat, Habgrid=Habgrid, G = G, 
                      nsites=nsites, M = Msite, B = B))
    
    # Write out the BUGS model file
    cat("
        model{ 
        
        # Prior distributions
        sigma ~ dunif(0,10)
        beta1 ~ dnorm(0,0.01)
        beta0 ~ dnorm(0,0.01)
        lam0 <- exp(beta0)*G           # Baseline lambda in terms of E(N) per sample unit
        
        # For each site, construct the DA parameter as a function of lambda
        for(s in 1:nsites){
        lamT[s] <- sum(lambda[,s])   # total abundance at a site
        psi[s] <- lamT[s]/M
        for(g in 1:G){             # g is the pixel index, there are G total pixels
        lambda[g,s] <- exp(beta0 + beta1*Habitat[g,s])
        probs[g,s] <- lambda[g,s]/sum(lambda[,s])
        }
        }
        
        # DA variables and spatial location variables:
        for(s in 1:nsites){
        for(i in 1:M){
        z[i,s] ~ dbern(psi[s])
        pixel[i,s] ~ dcat(probs[,s])
        u[i,s,1:2] <- Habgrid[pixel[i,s],]   # location = derived quantity  
        # distance = derived quantity
        d[i,s] <- pow(   pow( u[i,s,1]-B,2) + pow(u[i,s,2]-B,2), 0.5)
        p[i,s] <- exp(-d[i,s]*d[i,s]/(2*sigma*sigma))  # Half-normal model
        mu[i,s] <- p[i,s]*z[i,s]
        y[i,s] ~ dbern(mu[i,s])    # Observation model
        }
        # Derived parameters
        N[s]<- sum(z[,s])            # Site specific abundance
        }
        Ntotal <- sum(N[])             # Total across all sites
        D <- Ntotal/(9*nsites)         # Density: point area = 9 ha 
        }
        ",fill=TRUE, file="spatialHDS.txt")
    
    
    # Inits and parameters saved
    zst <- Ymat
    inits <- function(){ list (sigma=1.5, z = zst, beta0 = -5, beta1=1 ) }
    params <- c("sigma", "Ntotal", "beta1", "beta0", "D", "lam0")
    
    # MCMC settings
    ni <- 12000   ;   nb <- 2000   ;   nthin <- 2   ;   nc <- 3
    
    # Call JAGS, check convergence and summarize the results
    out3 <- jags(data, inits, params, "spatialHDS.txt", n.thin=nthin, n.chains=nc, n.burnin=nb, n.iter=ni, parallel = TRUE)
    par(mfrow = c(2,3))   ;   traceplot(out3)
    print(out3, 3)
    
    Ntrue <- tmp$N
    Nhat <- out3$summary[7:106,1]  
    plot(Ntrue, Nhat, xlab="Local population size, N, for each site", ylab="Posterior mean N for each site", pch=20)
    abline(0, 1, lwd=2)
    
    
    
    # 9.9 Summary
    # ------------------------------------------------------------------------