

#### E. ferrugineus #### 
obs_fer <- read.csv("C:/Users/voeroesd/Dropbox/EBD/Loros Patagonia/pat_obs_fer.csv") # Francisco Windows
#obs_fer <- read.csv("~/Dropbox/EBD/Loros Patagonia/pat_obs_fer.csv") # Francisco Mac

x <- as.data.frame(obs_fer)
### Exploratory plots ####
head(x)

## Some notes on dataset
# each row is a group
# count is group size
# season is breeding and non-breeding
# three habitat types: agricultural (pastures, crops), urban and others (includes different forest types)
# transect length (t.length) in km
# detection distance in m

plot(table(x$count))
plot(table(subset(x,count==2)$distance)) # detection distance for groups=2
plot(table(subset(x,count!=2)$distance)) # detection distance for groups other than 2


# new variable identifying pairs from other groups
x$pair<-0
x$pair[which(x$count==2)]<-1

## Estimating effective detection radius (EDR) ####

x$distance[x$distance == 500] <- 300 # correct outlier
library(mefa4)
library(detect)

br <- seq(0, 320, by=20)
x$bin <- cut(x$distance, br,  include.lowest = TRUE)
table(x$bin)
x$pres <- 1

x$groupid <- 0
for (i in unique(x$site)) {
    ii <- x$site == i
    x[ii, "groupid"] <- seq_len(sum(ii))
}
x$site_g <- paste0(x$site, "_", x$groupid)
head(x[order(x$site),], 25)

# Matrix with binned detection distances in each site
Y <- as.matrix(Xtab(pres ~site + bin, x))

# Dataframe for covariate vectors
## sum of all inds in all groups
X <- data.frame(ntot=rowSums(Xtab(count ~site + bin, x)))
## number of pairs / site
tmp <- aggregate(x[,c("pair"),drop=FALSE], list(Site=x$site), sum)
stopifnot(all(rownames(X) == as.character(tmp$Site)))
X$pair <- tmp$pair
## number of groups / site
tmp <- aggregate(x[,c("groupid"),drop=FALSE], list(Site=x$site), max)
X$ngroups <- tmp$groupid
## average group size = ntot / ngroups
X$gavg <- X$ntot / X$ngroups
## the rest is just repeated, so we take the unique values
tmp <- nonDuplicated(x, site, TRUE)
X <- data.frame(X, tmp[rownames(X), c("jdate", "Urban", "Agropastoral","site")])

head(X, n=25)

D <- matrix(br[-1], nrow(Y), length(br)-1, byrow=TRUE)
rownames(X) <- rownames(D) <- rownames(Y)
colnames(D) <- colnames(Y)


# Null model
m1 <- cmulti(Y | D ~1, type="dis")
summary(m1)

# Average group size
m2 <- cmulti(Y | D ~ gavg, X, type="dis")
summary(m2)

# Number of pairs
m3 <- cmulti(Y | D ~ pair, X, type = "dis")
summary(m3)

m4 <- cmulti(Y | D ~ ngroups, X, type = "dis")
summary(m4)

m5 <- cmulti(Y | D ~ Urban + Agropastoral, X, type = "dis")
summary(m5)

# There is a significant effect of habitat type on EDR

AIC(m1,m2,m3,m4,m5)[order(AIC(m1,m2,m3,m4,m5)$AIC),]


# Others
exp(coef(m5))[1]

# Urban
exp(sum(coef(m5)[1:2]))

# Agropastoral
exp(sum(coef(m5)[c(1,3)]))


## Estimating area surveyed for each site, in km^2
sites <- read.csv("C:/Users/voeroesd/Dropbox/EBD/Loros Patagonia/pat_site.csv")
#sites <- read.csv("~/Dropbox/EBD/Loros Patagonia/pat_site.csv")

sites$A <- sites$habitat.length.km*(2*exp(coef(m5))[1]/1000)
for (i in 1:nrow(sites)){
  if (sites$habitat[i]=="Urban"){
    sites$A[i] <- sites$habitat.length.km[i]*(2*exp(sum(coef(m5)[1:2]))/1000)
  }
  if (sites$habitat[i]=="Agropastoral"){
    sites$A[i] <- sites$habitat.length.km[i]*(2*exp(sum(coef(m5)[c(1,3)]))/1000)
  }
}

head(sites,n=50)


# Model for number of groups: Gi ~ Poisson(DiAi), where Di = covariates and Ai = area sampled in site
str(sites)
summary(sites)
sites$ngroups<- 0
sites$ngroups[which(sites$site%in%X$site)]<-X$ngroups

table(sites$ngroups)

boxplot(ngroups~habitat,data=sites)

# some exploratory plots
test<- with(sites,table(ngroups,habitat))
test<-as.matrix(as.data.frame.matrix(test))
plot( col(test),
      row(test),
      cex=log(test),#log-transformed sample size
      xlim=c(0.5,ncol(test)+0.5),
      ylim=c(0.5,nrow(test)+0.5),
      axes=FALSE,
      ann=FALSE
)
axis(1,at=1:ncol(test),labels=colnames(test),cex.axis=0.8)
axis(2,at=1:nrow(test),labels=rownames(test),cex.axis=0.8)
title(xlab="Habitat",ylab="Number of groups")

test2<- with(sites,table(ngroups,season))
test2<-as.matrix(as.data.frame.matrix(test2))
plot( col(test2),
      row(test2),
      cex=log(test2),#log-transformed sample size
      xlim=c(0.5,ncol(test2)+0.5),
      ylim=c(0.5,nrow(test2)+0.5),
      axes=FALSE,
      ann=FALSE
)
axis(1,at=1:ncol(test2),labels=colnames(test2),cex.axis=0.8)
axis(2,at=1:nrow(test2),labels=rownames(test2),cex.axis=0.8)
title(xlab="Season",ylab="Number of groups")

test3<- with(sites,table(ngroups,year))
test3<-as.matrix(as.data.frame.matrix(test3))
plot( col(test3),
      row(test3),
      cex=log(test3),#log-transformed sample size
      xlim=c(0.5,ncol(test3)+0.5),
      ylim=c(0.5,nrow(test3)+0.5),
      axes=FALSE,
      ann=FALSE
)
axis(1,at=1:ncol(test3),labels=colnames(test3),cex.axis=0.8)
axis(2,at=1:nrow(test3),labels=rownames(test3),cex.axis=0.8)
title(xlab="Year",ylab="Number of groups")


plot(ngroups~elevation, data=sites)
plot(ngroups~A, data=sites)
plot(ngroups~jdate, data=sites)

sites <- sites[!is.na(sites$A) & sites$habitat.length.km>0,]

## source diagnostic functions from Peter's code
source("diagnostics-functions.R")
#siplot(ngroups ~ ., sites)



# linear effects of Area, habitat and elevation
glm1 <- glm(ngroups~habitat+elevation, family=poisson, data=sites, offset=log(sites$A))
summary(glm1)
anova(glm1)

# quadratic effect of elevation
glm2 <- glm(ngroups~habitat+elevation+I(elevation^2), family=poisson, data=sites, offset=log(sites$A))
summary(glm2)
anova(glm2)

AIC(glm1,glm2) # model with quadratic elevation has better fit

# test only habitat vs. only elevation
glm3 <- glm(ngroups~habitat, family=poisson, data=sites, offset=log(sites$A))
summary(glm3)
anova(glm3)

glm4 <- glm(ngroups~elevation+I(elevation^2), family=poisson, data=sites, offset=log(sites$A))
summary(glm4)
anova(glm4)

AIC(glm2,glm3,glm4)

# The model with both covariates is best

# enter within-year temporal covariates
glm3 <- glm(ngroups~habitat+elevation+I(elevation^2)+season, family=poisson, data=sites, offset=log(sites$A))
summary(glm3)

glm4 <- glm(ngroups~habitat+elevation+I(elevation^2)+jdate, family=poisson, data=sites, offset=log(sites$A))
summary(glm4)

# test interaction between season and habitat and jdate and habitat
glm5 <- glm(ngroups~elevation+I(elevation^2)+habitat*season, family=poisson, data=sites, offset=log(sites$A))
summary(glm5)

glm6 <- glm(ngroups~elevation+I(elevation^2)+habitat*jdate, family=poisson, data=sites, offset=log(sites$A))
summary(glm6)

AIC(glm2,glm3,glm4,glm5,glm6)[order(AIC(glm2,glm3,glm4,glm5,glm6)$AIC),] # although the model with season*habitat interaction has a slightly better fit, the interaction is not significant, so continue with model with season and no interaction (glm3)

# enter year covariate
glm7 <- glm(ngroups~habitat+elevation+I(elevation^2)+season+year, family=poisson, data=sites, offset=log(sites$A))
summary(glm7)

AIC(glm3,glm7)

# model with year has better fit
anova(glm7)

mep(glm7) # some diagnostic plots


#### Models for group size ####
plot(table(x$count))
x$count
str(x)
x$A <- NA
  
  sites$A[which(sites$site%in%x$site)]

for(i in 1:nrow(x)){
  x$A[i] <- sites$A[which(sites$site==x$site[i])]
}

X <- model.matrix(~Urban+Agropastoral,x)
Z <- model.matrix(~season, x)


mod. <- vip(Y=x$count, X=X, offsetx=log(x$A), Z=Z, V=2, truncate=TRUE, hessian=TRUE, method="SANN")
 
summary(mod.)


#### # E. leptorhynchus ####
obs_lep <- read.csv("C:/Users/voeroesd/Dropbox/EBD/Loros Patagonia/pat_obs_lep.csv")
x <- as.data.frame(obs_lep)
### Exploratory plots ####
head(x)

## Some notes on dataset
# each row is a group
# count is group size
# season is breeding and non-breeding
# three habitat types: agricultural (pastures, crops), urban and others (includes different forest types)
# transect length (t.length) in km
# detection distance in m

plot(table(x$count))
plot(table(subset(x,count==2)$distance)) # detection distance for groups=2
plot(table(subset(x,count!=2)$distance)) # detection distance for groups other than 2


# new variable identifying pairs from other groups
x$pair<-0
x$pair[which(x$count==2)]<-1

## Estimating effective detection radius (EDR) ####

x$distance[x$distance == 2000] <- 1000 # correct outlier
library(mefa4)
library(detect)

br <- seq(0, 1020, by=20)
x$bin <- cut(x$distance, br,  include.lowest = TRUE)
table(x$bin)
x$pres <- 1

x$groupid <- 0
for (i in unique(x$site)) {
  ii <- x$site == i
  x[ii, "groupid"] <- seq_len(sum(ii))
}
x$site_g <- paste0(x$site, "_", x$groupid)
head(x[order(x$site),], 25)

# Matrix with binned detection distances in each site
Y <- as.matrix(Xtab(pres ~site + bin, x))

# Dataframe for covariate vectors
## sum of all inds in all groups
X <- data.frame(ntot=rowSums(Xtab(count ~site + bin, x)))
## number of pairs / site
tmp <- aggregate(x[,c("pair"),drop=FALSE], list(Site=x$site), sum)
stopifnot(all(rownames(X) == as.character(tmp$Site)))
X$pair <- tmp$pair
## number of groups / site
tmp <- aggregate(x[,c("groupid"),drop=FALSE], list(Site=x$site), max)
X$ngroups <- tmp$groupid
## average group size = ntot / ngroups
X$gavg <- X$ntot / X$ngroups
## the rest is just repeated, so we take the unique values
tmp <- nonDuplicated(x, site, TRUE)
X <- data.frame(X, tmp[rownames(X), c("jdate", "Urban", "Agropastoral","site")])

head(X, n=25)

D <- matrix(br[-1], nrow(Y), length(br)-1, byrow=TRUE)
rownames(X) <- rownames(D) <- rownames(Y)
colnames(D) <- colnames(Y)


# Null model
m1 <- cmulti(Y | D ~1, type="dis")
summary(m1)

m2 <- cmulti(Y | D ~ gavg, X, type="dis")
summary(m2)

m3 <- cmulti(Y | D ~ pair, X, type = "dis")
summary(m3)

m4 <- cmulti(Y | D ~ ngroups, X, type = "dis")
summary(m4)

m5 <- cmulti(Y | D ~ Urban + Agropastoral, X, type = "dis")
summary(m5)

AIC(m1,m2,m3,m4,m5)[order(AIC(m1,m2,m3,m4,m5)$AIC),]

# Others
exp(coef(m5))[1]

# Urban
exp(sum(coef(m5)[1:2]))

# Agropastoral
exp(sum(coef(m5)[c(1,3)]))



## Estimating area surveyed for each site, in km^2
sites <- read.csv("C:/Users/voeroesd/Dropbox/EBD/Loros Patagonia/pat_site.csv")
#sites <- read.csv("~/Dropbox/EBD/Loros Patagonia/pat_site.csv")

sites$A <- sites$habitat.length.km*(2*exp(coef(m5))[1]/1000)
for (i in 1:nrow(sites)){
  if (sites$habitat[i]=="Urban"){
    sites$A[i] <- sites$habitat.length.km[i]*(2*exp(sum(coef(m5)[1:2]))/1000)
  }
  if (sites$habitat[i]=="Agropastoral"){
    sites$A[i] <- sites$habitat.length.km[i]*(2*exp(sum(coef(m5)[c(1,3)]))/1000)
  }
}

head(sites,n=50)


# Model for number of groups: Gi ~ Poisson(DiAi), where Di = covariates and Ai = area sampled in site
str(sites)
summary(sites)
sites$ngroups<- 0
sites$ngroups[which(sites$site%in%X$site)]<-X$ngroups

table(sites$ngroups)

boxplot(ngroups~habitat,data=sites)

# some exploratory plots
test<- with(sites,table(ngroups,habitat))
test<-as.matrix(as.data.frame.matrix(test))
plot( col(test),
      row(test),
      cex=log(test),#log-transformed sample size
      xlim=c(0.5,ncol(test)+0.5),
      ylim=c(0.5,nrow(test)+0.5),
      axes=FALSE,
      ann=FALSE
)
axis(1,at=1:ncol(test),labels=colnames(test),cex.axis=0.8)
axis(2,at=1:nrow(test),labels=rownames(test),cex.axis=0.8)
title(xlab="Habitat",ylab="Number of groups")

test2<- with(sites,table(ngroups,season))
test2<-as.matrix(as.data.frame.matrix(test2))
plot( col(test2),
      row(test2),
      cex=log(test2),#log-transformed sample size
      xlim=c(0.5,ncol(test2)+0.5),
      ylim=c(0.5,nrow(test2)+0.5),
      axes=FALSE,
      ann=FALSE
)
axis(1,at=1:ncol(test2),labels=colnames(test2),cex.axis=0.8)
axis(2,at=1:nrow(test2),labels=rownames(test2),cex.axis=0.8)
title(xlab="Season",ylab="Number of groups")

test3<- with(sites,table(ngroups,year))
test3<-as.matrix(as.data.frame.matrix(test3))
plot( col(test3),
      row(test3),
      cex=log(test3),#log-transformed sample size
      xlim=c(0.5,ncol(test3)+0.5),
      ylim=c(0.5,nrow(test3)+0.5),
      axes=FALSE,
      ann=FALSE
)
axis(1,at=1:ncol(test3),labels=colnames(test3),cex.axis=0.8)
axis(2,at=1:nrow(test3),labels=rownames(test3),cex.axis=0.8)
title(xlab="Year",ylab="Number of groups")


plot(ngroups~elevation, data=sites)
plot(ngroups~A, data=sites)
plot(ngroups~jdate, data=sites)

sites <- sites[!is.na(sites$A) & sites$habitat.length.km>0,]

## source diagnostic functions from Peter's code
source("diagnostics-functions.R")
#siplot(ngroups ~ ., sites)



# linear effects of Area, habitat and elevation
glm1 <- glm(ngroups~habitat+elevation, family=poisson, data=sites, offset=log(sites$A))
summary(glm1)
anova(glm1)

# quadratic effect of elevation
glm2 <- glm(ngroups~habitat+elevation+I(elevation^2), family=poisson, data=sites, offset=log(sites$A))
summary(glm2)
anova(glm2)

AIC(glm1,glm2) # model with quadratic elevation has equal fit

# test only habitat vs. only elevation
glm3 <- glm(ngroups~habitat, family=poisson, data=sites, offset=log(sites$A))
summary(glm3)
anova(glm3)

glm4 <- glm(ngroups~elevation+I(elevation^2), family=poisson, data=sites, offset=log(sites$A))
summary(glm4)
anova(glm4)

AIC(glm2,glm3,glm4)

# The model with both covariates is best

# enter within-year temporal covariates
glm3 <- glm(ngroups~habitat+elevation+I(elevation^2)+season, family=poisson, data=sites, offset=log(sites$A))
summary(glm3)

glm4 <- glm(ngroups~habitat+elevation+I(elevation^2)+jdate, family=poisson, data=sites, offset=log(sites$A))
summary(glm4)

# test interaction between season and habitat and jdate and habitat
glm5 <- glm(ngroups~elevation+I(elevation^2)+habitat*season, family=poisson, data=sites, offset=log(sites$A))
summary(glm5)

glm6 <- glm(ngroups~elevation+I(elevation^2)+habitat*jdate, family=poisson, data=sites, offset=log(sites$A))
summary(glm6)

AIC(glm2,glm3,glm4,glm5,glm6)[order(AIC(glm2,glm3,glm4,glm5,glm6)$AIC),] # model with jdate*habitat interaction has a better fit

# enter year covariate
glm7 <- glm(ngroups~elevation+I(elevation^2)+habitat*jdate+year, family=poisson, data=sites, offset=log(sites$A))
summary(glm7)

AIC(glm6,glm7)

# model with year has better fit
anova(glm7)

mep(glm7) # some diagnostic plots


#### Models for group size ####
plot(table(x$count))
x$count
str(x)
x$A <- NA

sites$A[which(sites$site%in%x$site)]

for(i in 1:nrow(x)){
  x$A[i] <- sites$A[which(sites$site==x$site[i])]
}

X <- model.matrix(~Urban+Agropastoral+A,x)
Z <- model.matrix(~jdate, x)


mod. <- vip(Y=x$count, X=X, Z=Z, V=2, truncate=TRUE)


