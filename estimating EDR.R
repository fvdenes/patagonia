obs_fer <- read.csv("~/Dropbox/collaborations/patagonia/pat_obs_fer.csv")
obs_fer <- read.csv("C:/Users/voeroesd/Dropbox/EBD/Loros Patagonia/pat_obs_fer.csv")

#obs_lep <- read.csv("C:/Users/voeroesd/Dropbox/EBD/Loros Patagonia/pat_obs_lep.csv")


# E. ferrugineus
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

x$distance[x$distance == 500] <- 300
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
X <- data.frame(X, tmp[rownames(X), c("jdate", "urban", "others")])

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


m4 <- cmulti(Y | D ~ jdate , X, type="dis")
summary(m4)

m5 <- cmulti(Y | D ~ urban + others, X, type = "dis")
summary(m5)

m6 <- cmulti(Y | D ~ ngroups, X, type = "dis")
summary(m6)

# There is a significant effect of habitat type on EDR

AIC(m1,m2,m3,m4,m5,m6)[order(AIC(m1,m2,m3,m4,m5,m6)$AIC),]


# Agricultural
exp(coef(m5))[1]

# Urban
exp(sum(coef(m5)[1:2]))

# Others
exp(sum(coef(m5)[c(1,3)]))


## Estimating area surveyed for each site, in km^2
sites <- read.csv("C:/Users/voeroesd/Dropbox/EBD/Loros Patagonia/pat_site.csv")

sites$A <- sites$habitat.length.km*(2*exp(coef(m5))[1]/1000)
for (i in 1:nrow(sites)){
  if (sites$habitat[i]=="Urbano"){
    sites$A[i] <- sites$habitat.length.km[i]*(2*exp(sum(coef(m5)[1:2]))/1000)
  }
  if (sites$habitat[i]=="Otros"){
    sites$A[i] <- sites$habitat.length.km[i]*(2*exp(sum(coef(m5)[c(1,3)]))/1000)
  }
}

head(sites,n=50)


