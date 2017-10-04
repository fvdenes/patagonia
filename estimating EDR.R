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

# Matrix with binned detection distances in each site
Y <- as.matrix(Xtab(pres ~site + bin, x))

# Dataframe for covariate vectors
X <- data.frame(
  gsize=rowSums(Xtab(count ~site + bin, x)), # what about sites with more than one group detected (e.g. 226)? this command sums the counts from all groups, irrespective of the distance at which each group was detected. Should they not be separated?
  pair=rowSums(Xtab(pair ~site + bin, x)), # not sure if this and lines below are correct...
  jdate=rowSums(Xtab(jdate ~site + bin, x)),
  urban=rowSums(Xtab(urban ~site + bin, x)), 
  others=rowSums(Xtab(others ~site + bin, x))
  )


D <- matrix(br[-1], nrow(Y), length(br)-1, byrow=TRUE)
rownames(X) <- rownames(D) <- rownames(Y)
colnames(D) <- colnames(Y)


m1 <- cmulti(Y | D ~1, type="dis")
summary(m1)
exp(coef(m1))
m2 <- cmulti(Y | D ~ log(gsize), X, type="dis")
summary(m2)

m3 <- cmulti(Y | D ~ pair, X, type = "dis")
summary(m3)

m4 <- cmulti(Y | D ~ log(gsize)*pair, X, type="dis")
summary(m4)

m5 <- cmulti(Y | D ~ urban + others, X, type = "dis")
summary(m5)

m6 <- cmulti(Y | D ~ log(jdate), X, type= "dis")
summary(m6)
