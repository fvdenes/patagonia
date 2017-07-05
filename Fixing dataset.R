##### Analysis of patagonia parrot count data #####

## Load excel data data ####
pat_obs <- read.csv("~/Dropbox/Pos-doc/Loros Patagonia/patagonia_obs.csv")
pat_site <- read.csv("~/Dropbox/Pos-doc/Loros Patagonia/patagonia_site.csv")

### Fixing datasets ####
str(pat_obs)

pat_obs$distance<-as.numeric(pat_obs$distance)

pat_obs$perpendicular[which(pat_obs$perpendicular=="")]<-NA
pat_obs$perpendicular<-as.numeric(pat_obs$perpendicular)

levels(pat_obs$species)
pat_obs$species[which(pat_obs$species=="E. Ferrugineus"|pat_obs$species=="E. ferrugineus ")]<-"E. ferrugineus"
pat_obs$species[which(pat_obs$species=="E. Leptorhynchus")]<-"E. leptorhynchus"
pat_obs$species<-factor(pat_obs$species)

pat_obs$count<-as.numeric(pat_obs$count)

#pat_obs$sample.unit<-as.character(pat_obs$sample.unit)

pat_obs$jdate<-strptime(pat_obs$date,format="%e/%m/%y")$yday+1
pat_obs$year<-strptime(pat_obs$date,format="%e/%m/%y")$year+1900

pat_obs$time[which(pat_obs$time=="")]<-NA
pat_obs$time <- as.character(pat_obs$time)
pat_obs$time <- sapply(strsplit(pat_obs$time,":"),
       function(x) {
         x <- as.numeric(x)
         x[1]+x[2]/60
       }
)

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


# removing updating sample unit identification
str(pat_site)
str(pat_obs)

(pat_site$site <- 1:nrow(pat_site))
(pat_obs$site <- pat_site$site[match(pat_obs$sample.unit,pat_site$sample.unit)])

head(pat_site,n=20)
head(pat_obs,n=20)

# remove observations of pat_obs without group size data, distance data

(pat_obs<- pat_obs[is.na(pat_obs$count)==FALSE,] )
(pat_obs<- pat_obs[is.na(pat_obs$distance)==FALSE,] )
(pat_obs <- pat_obs[-which(pat_obs$count==0),]) # remove obs with count=0

# Convert habitat covariate into model matrix
levels(pat_site$habitat)
habitat <- as.character(pat_site$habitat)
habitat[which(habitat=="agroganadero")]<- "Agroganadero"
pat_site$habitat<- as.factor(habitat)
pat_site$habitat<-relevel(pat_site$habitat,"Araucaria")

my_matrix <- model.matrix(~ pat_site$habitat)
head(my_matrix)
my_matrix <- data.frame(my_matrix[,-1])

colnames(my_matrix) <- unlist(strsplit(colnames(my_matrix),"habitat"))[seq(2,16,by=2)]

pat_site<-cbind(pat_site,my_matrix)

head(pat_site,n=20)



# Write data files
write.csv(pat_obs[which(pat_obs$species=="E. leptorhynchus"),c(13,3,11,12,6:10)],"pat_obs_lep.csv",row.names = FALSE)
write.csv(pat_obs[which(pat_obs$species=="E. ferrugineus"),c(13,3,11,12,6:10)],"pat_obs_fer.csv",row.names = FALSE)
write.csv(pat_site[,c(7,2,4:6,8:15)],"pat_site.csv",row.names = FALSE)

rm(pat_obs,pat_site,my_matrix,habitat)

#### Load count and site data ####

site.data <- read.csv("~/Dropbox/Pos-doc/Loros Patagonia/pat_site.csv")
obs_fer <- read.csv("~/Dropbox/Pos-doc/Loros Patagonia/pat_obs_fer.csv")
obs_lep <- read.csv("~/Dropbox/Pos-doc/Loros Patagonia/pat_obs_lep.csv")