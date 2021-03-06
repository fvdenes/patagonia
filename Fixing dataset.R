##### Analysis of patagonia parrot count data #####

## Load excel data data ####
#mac
#pat_obs <- read.csv("~/Dropbox/EBD/Loros Patagonia/patagonia_obs.csv")
#pat_site <- read.csv("~/Dropbox/EBD/Loros Patagonia/patagonia_site.csv")
#windows
pat_obs <- read.csv("C:/Users/voeroesd/Dropbox/EBD/Loros Patagonia/patagonia_obs.csv")
pat_site <- read.csv("C:/Users/voeroesd/Dropbox/EBD/Loros Patagonia/patagonia_site.csv")


### Fixing datasets ####
str(pat_obs)
str(pat_site)

pat_obs$distance<-as.numeric(pat_obs$distance)

pat_obs$perpendicular[which(pat_obs$perpendicular=="")]<-NA
pat_obs$perpendicular<-as.numeric(pat_obs$perpendicular)

levels(pat_obs$species)
pat_obs$species[which(pat_obs$species=="E. Ferrugineus"|pat_obs$species=="E. ferrugineus ")]<-"E. ferrugineus"
pat_obs$species[which(pat_obs$species=="E. Leptorhynchus")]<-"E. leptorhynchus"
pat_obs$species<-factor(pat_obs$species)

pat_obs$count<-as.numeric(pat_obs$count)

#pat_obs$sample.unit<-as.character(pat_obs$sample.unit)

#### Assign dates to pat_obs ####
pat_obs$date<-as.character(pat_obs$date)
pat_obs$sample.unit<-as.character(pat_obs$sample.unit)
pat_site$date<-as.character(pat_site$date)
pat_obs$sample.unit<-as.character(pat_obs$sample.unit)

for(i in 1:nrow(pat_obs)){
  pat_obs$date[i] <-pat_site$date[which(pat_site$sample.unit==pat_obs$sample.unit[i])]
} 

head(pat_obs,n=100)

pat_obs$jdate<-strptime(pat_obs$date,format="%Y-%m-%d")$yday+1
pat_obs$year<-strptime(pat_obs$date,format="%Y-%m-%d")$year+1900
pat_obs$month<-as.numeric(strftime(pat_obs$date,format="%m"))

pat_site$jdate<-strptime(pat_site$date,format="%Y-%m-%d")$yday+1
pat_site$year<-strptime(pat_site$date,format="%Y-%m-%d")$year+1900
pat_site$month<-as.numeric(strftime(pat_site$date,format="%m"))

# assign season: breeding (11,12) and non-breeding (rest)
pat_obs$season<-"non-breeding"
pat_obs$season[which(pat_obs$month==11&12)]<-"breeding"

pat_site$season<-"non-breeding"
pat_site$season[which(pat_site$month==11&12)]<-"breeding"

pat_obs$time[which(pat_obs$time=="")]<-NA
pat_obs$time <- as.character(pat_obs$time)
pat_obs$time <- sapply(strsplit(pat_obs$time,":"),
       function(x) {
         x <- as.numeric(x)
         x[1]+x[2]/60
       }
)

pat_obs$time[which(pat_obs$year==2013&2014)]<-NA # Fix time variable (counts from 2013 and 2014 have weird times)

pat_obs



### Fixing site-level dataset

pat_site$elevation[which(pat_site$elevation=="#DIV/0!"|pat_site$elevation=="")]<-NA
pat_site$elevation[which(pat_site$elevation=="")]<-NA
pat_site$elevation<-as.numeric(pat_site$elevation)
pat_site$elevation[is.na(pat_site$elevation)] <-mean(pat_site$elevation,na.rm=T) # fill NA values with mean elevation
pat_site$elevation<-round(pat_site$elevation,0)

pat_site$habitat<-as.factor(pat_site$habitat)
#pat_site$sample.unit<-as.character(pat_site$sample.unit)
pat_site$transect<-as.factor(pat_site$transect)

head(pat_site)

# remove sites with habitat.length.km=0
which(pat_site$habitat.length.km==0)
pat_site<-pat_site[-which(pat_site$habitat.length.km==0),]

# removing updating sample unit identification
str(pat_site)
str(pat_obs)

(pat_site$site <- 1:nrow(pat_site))
(pat_obs$site <- pat_site$site[match(pat_obs$sample.unit,pat_site$sample.unit)])

head(pat_site,n=20)
head(pat_obs,n=20)

# remove observations of pat_obs without group size data, distance data, site
(pat_obs<- pat_obs[is.na(pat_obs$site)==FALSE,] )
(pat_obs<- pat_obs[is.na(pat_obs$count)==FALSE,] )
(pat_obs<- pat_obs[is.na(pat_obs$distance)==FALSE,] )
(pat_obs <- pat_obs[-which(pat_obs$count==0),]) # remove obs with count=0

# Convert habitat covariate into model matrix
levels(pat_site$habitat)
habitat <- as.character(pat_site$habitat)
pat_site$habitat[which(pat_site$habitat=="agroganadero")]<- "Agroganadero"


pat_site$habitat2<-pat_site$habitat # backup original habitat classification before simplification

# Reducing habitat levels:
# Araucaria + Aramix + Estepa + Matorral + Otros = Other
pat_site$habitat<- as.character(pat_site$habitat)
pat_site$habitat[which(pat_site$habitat=="Estepa")]<-"Other"
pat_site$habitat[which(pat_site$habitat=="Matorral")]<-"Other"
pat_site$habitat[which(pat_site$habitat=="Araucaria")]<-"Other"
pat_site$habitat[which(pat_site$habitat=="Aramix")]<-"Other"
pat_site$habitat[which(pat_site$habitat=="Noto")]<-"Other"
pat_site$habitat[which(pat_site$habitat=="Otros")]<-"Other"

# Agroganadero + Agroforestal = Agropastoral
pat_site$habitat[which(pat_site$habitat=="Agroforestal")]<-"Agropastoral"
pat_site$habitat[which(pat_site$habitat=="Agroganadero")]<-"Agropastoral"

# urbano = Urban
pat_site$habitat[which(pat_site$habitat=="Urbano")]<-"Urban"

pat_site$habitat <- as.factor(pat_site$habitat)
pat_site$habitat<-relevel(pat_site$habitat,"Other")


my_matrix <- model.matrix(~ pat_site$habitat)
head(my_matrix)
my_matrix <- data.frame(my_matrix[,-1])

colnames(my_matrix) <- unlist(strsplit(colnames(my_matrix),"habitat"))[seq(2,4,2)]

pat_site<-cbind(pat_site,my_matrix)

head(pat_site,n=20)


# put site-level data in observation dataset
head(pat_obs,n=100)
pat_obs$t.length<-NA
pat_obs$elevation<- NA
pat_obs$Agropastoral<-NA
pat_obs$Urban<-NA

for(i in 1:nrow(pat_obs)){
  pat_obs$t.length[i]<-pat_site$habitat.length.km[which(pat_site$site==pat_obs$site[i])]
  pat_obs$elevation[i]<- pat_site$elevation[which(pat_site$site==pat_obs$site[i])]
  pat_obs$Agropastoral[i]<- pat_site$Agropastoral[which(pat_site$site==pat_obs$site[i])]
  pat_obs$Urban[i]<-pat_site$Urban[which(pat_site$site==pat_obs$site[i])]
}

str(pat_obs)
str(pat_site)




# Write data files
write.csv(pat_obs[which(pat_obs$species=="E. leptorhynchus"),c(14,3,11,12,6:9,15:18)],"pat_obs_lep.csv",row.names = FALSE)
write.csv(pat_obs[which(pat_obs$species=="E. ferrugineus"),c(14,3,11,12,6:9,15:18)],"pat_obs_fer.csv",row.names = FALSE)
write.csv(pat_site[,c(12,4:6,8:9,11,14:15)],"pat_site.csv",row.names = FALSE)

rm(pat_obs,pat_site,my_matrix,habitat,i)

#### Load count and site data ####

site.data <- read.csv("~/Dropbox/EBD/Loros Patagonia/pat_site.csv")
obs_fer <- read.csv("~/Dropbox/EBD/Loros Patagonia/pat_obs_fer.csv")
obs_lep <- read.csv("~/Dropbox/EBD/Loros Patagonia/pat_obs_lep.csv")
