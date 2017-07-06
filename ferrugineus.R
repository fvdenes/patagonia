

## ferrugineus
site.data <- read.csv("C:/Users/voeroesd/Dropbox/EBD/Loros Patagonia/pat_site.csv")
obs_fer <- read.csv("C:/Users/voeroesd/Dropbox/EBD/Loros Patagonia/pat_obs_fer.csv")

# Site-level covariates
elevation <- site.data$elevation
agroflorestal <- site.data$Agroforestal
agroganadero <- site.data$Agroganadero
#aramix <- site.data$Aramix
#estepa <- site.data$Estepa
#matorral <- site.data$Matorral
#noto <- site.data$Noto
#otros <- site.data$Otros
#urbano <- site.data$Urbano

nsites <- nrow(site.data)

# Observation-level covariates
data <- obs_fer

groupsize <- data[,"count"] -1    # Input groupsize-1 as data

B <- 510 # strip half-width

delta <- 10 # width of distance bin
d <- data$distance
dclass <- d%/%delta+1 # assign observed distances to distance intervals
xg <- seq(0,B,delta)
midpt <- xg[-1] - delta/2 # get midpoint of distance intervals
nD <- length(midpt) # how many intervals

M <- 4000                        # Size of augmented data set is M
nz <- M-nrow(data)              # Number of "pseudo-groups" added
y <- c(rep(1,nrow(data)),rep(0,nz))      # Indicator of capture (== 1 for all obs. groups)
nind <- nrow(data)              # Number of observed groups
site <- c(data[,1], rep(NA,nz)) # Site they belong to is unknown 

groupsize <- c(groupsize, rep(NA,nz)) # .... as is their size
zst <- y                        # Starting values for data augmentation variable

dclass <- c(dclass,rep(NA,nz))


# (few habitats, just to see if code runs)

# Bundle data and produce summary
str(bugs.data <- list (y=y, B=B, nind=nind, nsites=nsites, midpt=midpt, delta=delta, dclass=dclass, site=site, nz=nz, nD=nD, groupsize=groupsize,elevation=elevation,agroflorestal=agroflorestal,agroganadero=agroganadero))


# Define model in BUGS language (few habitats, just to see if code runs)
cat("
    model{ 
    
    # Prior distributions for model parameters
    alpha0 ~ dunif(-2,2)              # groupsize intercept
    alpha1 ~ dunif(-5,5)              # groupsize slope
    beta0 ~ dunif(-10,10)             # intercept for log(lambda)
    beta1 ~ dunif(-10,10)             # elevation slope
    beta2 ~ dunif(-10,10)             # agroflorestal slope
    beta3 ~ dunif(-10,10)             # agroganadero slope
  
    lambda.group ~ dgamma(0.1, 0.1)
    
    
    # psi is a derived parameter for sites
    psi <- sum(lambda[])/(nind+nz)
    
    # Individual level model: observations and process
    
    for(i in 1:(nind+nz)){
    
    z[i] ~ dbern(psi)                    # Data augmentation variables
    dclass[i] ~ dcat(pi.probs[site[i],]) # Population distribution of distance class
    mu[i] <- z[i]*p[site[i],dclass[i]]   # p depends on site and distance class
    groupsize[i] ~ dpois(lambda.group)   # Group size is Poisson
    
    log(sigma[i]) <- alpha0 +  alpha1*groupsize[i]
    y[i] ~ dbern(mu[i])
    site[i] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
    zg[i]<- z[i]*(groupsize[i] + 1)      # Number of individuals in that group
    }
    
    for(s in 1:nsites){                  # Loop over sites
    # Construct cell probabilities for nD cells (distance classes)
    for (g in 1:nD){                     # midpt[g] = midpoint of each cell; g is each distance class  
    log(p[s,g]) <- -midpt[g]*midpt[g]/(2*sigma[s]*sigma[s])
    pi[s,g] <- ((2*midpt[g])/(B*B))*delta
    pi.probs[s,g] <- pi[s,g]/norm[s]
    f[s,g] <- p[s,g] * pi[s,g]
    fc[s,g] <- f[s,g]/pcap[s]          # Conditional probabilities
    }
    pcap[s] <- sum(f[s,])                # Capture prob. is the sum of all rectangular areas
    norm[s] <- sum(pi[s,])
    
    # Model for population size of groups
    N[s] ~ dpois(lambda[s])
    log(lambda[s])<- beta0 + beta1*elevation[s] + beta2*agroflorestal[s] + beta3*agroganadero[s] 
    site.probs[s]<- lambda[s]/sum(lambda[])
    }
    
    # Derived quantities
    G <- sum(z[])        # Total number of groups
    Ntotal <- sum(zg[])  # Total population size (all groups combined)
    }
    ",fill=TRUE, file="model1.txt")


# Load some libraries, define MCMC settings, inits function and parameters to save

library("jagsUI")  # 
ni <- 3000   ;   nb <- 500   ;   nt <- 2   ;   nc <- 3
inits <- function(){list(alpha0=1, alpha1=0.5, beta0=1, beta1=-0.5,beta2=-0.5, z=zst)}
params <- c("alpha0", "alpha1", "beta0", "beta1", "beta2", "psi", "Ntotal", "G", "lambda.group")


# Call JAGS, check convergence and summarize posterior distributions
out1 <- jags(bugs.data, inits, params, "model1.txt", n.thin=nt,  n.chains=nc, n.burnin=nb,n.iter=ni, parallel=FALSE)
traceplot(out1)    
print(out1, 3)



# (all habitat types)
# Bundle data and produce summary
str(bugs.data <- list (y=y, B=B, nind=nind, nsites=nsites, midpt=midpt, delta=delta, dclass=dclass, site=site, nz=nz, nD=nD, groupsize=groupsize,elevation=elevation,agroflorestal=agroflorestal,agroganadero=agroganadero,aramix=aramix,estepa=estepa,matorral=matorral,noto=noto,otros=otros,urbano=urbano ))

# Define model in BUGS language (all habitat types)
cat("
    model{ 
    
    # Prior distributions for model parameters
    alpha0 ~ dunif(-2,2)              # groupsize intercept
    alpha1 ~ dunif(-5,5)              # groupsize slope
    beta0 ~ dunif(-10,10)             # intercept for log(lambda)
    beta1 ~ dunif(-10,10)             # elevation slope
    beta2 ~ dunif(-10,10)             # agroflorestal slope
    beta3 ~ dunif(-10,10)             # agroganadero slope
    beta4 ~ dunif(-10,10)             # aramix slope 
    beta5 ~ dunif(-10,10)             # estepa slope
    beta6 ~ dunif(-10,10)             # matorral slope
    beta7 ~ dunif(-10,10)             # noto slope
    beta8 ~ dunif(-10,10)             # otros slope
    beta9 ~ dunif(-10,10)             # urbano slope
    lambda.group ~ dgamma(0.1, 0.1)
    

    # psi is a derived parameter for sites
    psi <- sum(lambda[])/(nind+nz)
    
    # Individual level model: observations and process
    
    for(i in 1:(nind+nz)){

      z[i] ~ dbern(psi)                    # Data augmentation variables
      dclass[i] ~ dcat(pi.probs[site[i],]) # Population distribution of distance class
      mu[i] <- z[i]*p[site[i],dclass[i]]   # p depends on site and distance class
      groupsize[i] ~ dpois(lambda.group)   # Group size is Poisson
    
      log(sigma[i]) <- alpha0 +  alpha1*groupsize[i]
      y[i] ~ dbern(mu[i])
      site[i] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
      zg[i]<- z[i]*(groupsize[i] + 1)      # Number of individuals in that group
    }
    
    for(s in 1:nsites){                  # Loop over sites
     # Construct cell probabilities for nD cells (distance classes)
     for (g in 1:nD){                     # midpt[g] = midpoint of each cell; g is each distance class  
      log(p[s,g]) <- -midpt[g]*midpt[g]/(2*sigma[s]*sigma[s])
      pi[s,g] <- ((2*midpt[g])/(B*B))*delta
      pi.probs[s,g] <- pi[s,g]/norm[s]
      f[s,g] <- p[s,g] * pi[s,g]
      fc[s,g] <- f[s,g]/pcap[s]          # Conditional probabilities
     }
     pcap[s] <- sum(f[s,])                # Capture prob. is the sum of all rectangular areas
     norm[s] <- sum(pi[s,])
    
     # Model for population size of groups
     N[s] ~ dpois(lambda[s])
     log(lambda[s])<- beta0 + beta1*elevation[s] + beta2*agroflorestal[s] + beta3*agroganadero[s] + beta4*aramix[s] + beta5*estepa[s] + beta6*matorral[s] + beta7*noto[s] + beta8*otros[s] + beta9*urbano[s]
     site.probs[s]<- lambda[s]/sum(lambda[])
    }
    
    # Derived quantities
    G <- sum(z[])        # Total number of groups
    Ntotal <- sum(zg[])  # Total population size (all groups combined)
    }
    ",fill=TRUE, file="model1.txt")


# Load some libraries, define MCMC settings, inits function and parameters to save

library("jagsUI")  # 
ni <- 3000   ;   nb <- 500   ;   nt <- 2   ;   nc <- 3
inits <- function(){list(alpha0=1, alpha1=0.5, beta0=1, beta1=-0.5,beta2=-0.5,beta3=-0.5,beta4=-0.5,beta5=-0.5,beta6=-0.5,beta7=-0.5,beta8=-0.5,beta9=-0.5, z=zst)}
params <- c("alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "beta8", "beta9", "psi", "Ntotal", "G", "lambda.group")

# Call JAGS, check convergence and summarize posterior distributions
out1 <- jags(bugs.data, inits, params, "model1.txt", n.thin=nt,  n.chains=nc, n.burnin=nb,n.iter=ni, parallel=FALSE)
traceplot(out1)    
print(out1, 3)

sum(temp$N.true) # = number of groups (G)
sum(temp$groupsize) # = number of individuals (Ntotal)
