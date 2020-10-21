#Libraries
######################
# Charge libraries
######################
library(sp)
library(car)
library(hSDM)
library(nlme)
library(dismo)
library(gstat)
library(spdep)
library(Hmisc)
library(raster)
library(GGally)
library(fields)
library(ggplot2)
library(maptools)
library(corrplot)


######################################
#Read data of one species
data <- read.delim("Hemiramphus_brasiliensis.txt")

#OR

data<- read.delim("Hyporhamphus_unifasciatus.txt")
str(data)

#Select only longitude, latitude and presences
data=as.data.frame(cbind(data$Longitude,data$Latitude,data$Ano))
colnames(data)=c("lon","lat","year")

#Remove NAs
summary(data)
to.remove <- which(!complete.cases(data))
data <- data[-to.remove,]
summary(data)
dim(data)

################################################################################
#which records are duplicates?
################################################################################
dups2 <- duplicated(data[, c("lon","lat")]);dups2

################################################################################
#Number of duplicates
################################################################################
sum(dups2)

################################################################################
#Remove duplicates
################################################################################
data <- data[!dups2, ]

#################################################################################
#               Environmental variables                              ##
##############################################################################
files<-(list.files(".", full.names=T, pattern=".asc"))#change directory
predictors <- stack(files);predictors

######################################
### --- Standardize predictors --- ###
######################################
predictors2<-scale(predictors)

###################################################################################3
#Suitability Prediction
Bathy<- raster(predictors2, 'Bathy')
PP<- raster(predictors2, 'PP')
Rugosity<- raster(predictors2, 'Rugosity')
SSS<- raster(predictors2, 'SSS')
SST<-raster(predictors2, "SST")
SST26<-raster(predictors2, "SST26")
SST85<-raster(predictors2, "SST85")

#convert in data frame
bathy<-as.data.frame(Bathy)
pp<-as.data.frame(PP)
rugosity<-as.data.frame(Rugosity)
sss<-as.data.frame(SSS)
sst<-as.data.frame(SST)
sst26<-as.data.frame(SST26)
sst85<-as.data.frame(SST85)

#join in a unique database
x<-coordinates(Bathy)[,1]
y<-coordinates(Bathy)[,2]
env<-cbind(x,y,bathy,pp,rugosity,sss,sst,sst26,sst85)

#############################################################
#Extract environmental values from locations
##############################################################
coords=cbind(data$lon, data$lat);coords
pres_vals <- extract(predictors2, coords)
head(pres_vals)

set.seed(500)
backgr <- randomPoints(predictors2,500)#generate psuedo-absence 

abs_vals <- extract(predictors2, backgr)
head(abs_vals)

#Make a new database
coords=rbind(coords, backgr)
vals=rbind(pres_vals,abs_vals)
abs=rep(0, nrow(abs_vals))
species<- c(rep(1, nrow(pres_vals)), rep(0, nrow(abs_vals)));species
sdmdata <- data.frame(cbind(species, vals, coords))
colnames(sdmdata)=c("pb", "Bathy","PP" ,"Rugosity", "SSS" , "SST" ,"SST26","SST85","lon", "lat" )
head(sdmdata)

#Remove NAs
to.remove <- which(!complete.cases(sdmdata))
sdmdata <- sdmdata[-to.remove,]
summary(sdmdata)

################################################################################
#               Exploration of the dataset
################################################################################

################################################################################
#Step 1: Check correlation among explicative variables. 
#Variables highly correlated (>0.80) can NOT be used together in the model.
################################################################################

matrix<-rcorr(as.matrix(sdmdata[,c(1:6)]), type = "pearson")

# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

jpeg("Corr_matrix_1.jpeg", width = 1462, height = 1694, res = 300)
corrplot(matrix$r, type="lower", tl.col = "black",method="number",
         p.mat = matrix$P, sig.level = 0.05)
# in Corrplot "X" are no significant variables. We look at correlation among variables

################################################################################
#Step 2: Check multicollinearity among variables
################################################################################
source("HighstatLib.r")
corvif(sdmdata[,c(2:6)])#Remove step by step variable with high VIF (>3)

#Make both data sets spatial objects
data<-sdmdata[,c(1,9,10)]
coordinates(data) <- c(2,3)
cfr.env <- SpatialPixelsDataFrame(points=env[c("x","y")],
                                  tol=0.94,
                                  data=env[,-c(1,2)])
fullgrid(cfr.env) <- TRUE

#Get the indices of cells where presences and absences have been observed.
cfr.env.rast <- stack(cfr.env)
pres <- extract(cfr.env.rast, SpatialPoints(data[data$pb==1,]),
                cellnumbers=TRUE)[,1]
abs <- extract(cfr.env.rast, SpatialPoints(data[data$pb==0,]),
               cellnumbers=TRUE)[,1]

# Make the data frame used in regressions
ncelltot <- length(cfr.env) # Including NULL cells
d <- data.frame(lon=coordinates(cfr.env)[,1],lat=coordinates(cfr.env)[,2],
                Y=rep(0,ncelltot),
                trials=rep(0,ncelltot),
                cell.orig=1:ncelltot,
                cfr.env@data)
d$Y[pres] <- 1
d$trials[c(pres,abs)] <- 1

#Remove NAs
to.remove <- which(!complete.cases(d))
d <- d[-to.remove,]
summary(d)

# Find cells' neighborhood with function 'adjacent' from the 'raster' package
nb <- dnearneigh(as.matrix(d[,1:2]), 0, 1.5)
table(card(nb))
winnb <- nb2WB(nb)

#Find cells' neighborhood with function 'adjacent' from the 'raster' package
sel.cell <- d$cell.orig
s.cell <- sort(unique(d$cell.orig))
d$cell <- match(d$cell.orig,s.cell)

#Make d a spatial object for later use
coordinates(d) <- c(1,2)

mod.hSDM.binomial.iCAR<-hSDM.binomial.iCAR(presences=d$Y[d$trials>0],
                                           trials=d$trials[d$trials>0],
                                           suitability=~Bathy+SST+SSS*Rugosity,#change variables
                                           spatial.entity=d$cell[d$trials>0],
                                           data=d[d$trials>0,],
                                           n.neighbors = winnb$num,
                                           neighbors = winnb$adj,
                                           suitability.pred=d,
                                           spatial.entity.pred=d$cell,
                                           burnin=5000,
                                           mcmc=50000, thin=5,
                                           beta.start=0,
                                           Vrho.start=0.1,
                                           priorVrho=10,
                                           mubeta=0, Vbeta=1.0E6,
                                           #shape=0.1, rate=0.1,
                                           #Vrho.max=1000,
                                           seed=1234, verbose=1, save.rho=0)

#Outputs
summary(mod.hSDM.binomial.iCAR$mcmc)
#plot(mod.hSDM.binomial.iCAR$mcmc)

#Put output together
out <- data.frame(d,pred=mod.hSDM.binomial.iCAR$theta.pred,
                  sp.ef=mod.hSDM.binomial.iCAR$rho.pred)
#Plot results
#Attribute predicted values to raster cells
out <- data.frame(mod.hSDM.binomial.iCAR$theta.pred)
colnames(out) <- c("pred")

#correlations between predicted and observed with the same set:(neg bin)
range(out$pred)
cor(out$pred,d$Y) 

#Plot results
coordinates(out) <- coordinates(d)
proj4string(out)=CRS("+init=epsg:4724")
out = spTransform(out,CRS("+init=epsg:4724"))

#raster to plot
gridded(out) = TRUE
out <- raster(out)

###################################################################################3
#####################################################
#  CUT WITH BATHY
##################################################
depth<- raster::raster(predictors, "Bathy")
depth=abs(depth)

matrix<- cbind(coordinates(depth), depth=getValues(depth))
I <- is.na(matrix[,3])
matrix<- matrix[!I,]
matrix<-as.data.frame(matrix)
new<- subset(matrix, matrix$depth <=1500)

xy <- cbind(new$x, new$y)
rast<- raster(xmn=-118, xmx=-27, ymn=-40, ymx=43, nrows=996, ncols=1092)#2  
p<- rasterize(xy, rast, new$depth, fun=max,na.rm=F)
p<-resample(p,predictors)
e<-extent(p)

sp<-crop(out,e)
sp=resample(sp,p)
sp<-raster::mask(sp,p)

#read shapefile coast
coast=readOGR(".", "america_simp")
ext=c(-118, -27, -40, 43)
coast=crop(coast,ext)

plot(y,col=tim.colors(100)[1:100],main=" ",zlim=c(0,1),axes=T)
plot(coast, add=TRUE,col='dark grey')
box()

#########
#Save prediction
matrix<- cbind(coordinates(out), pred=getValues(out))
I <- is.na(matrix[,3])
matrix<- matrix[!I,]
colnames(matrix) <- c("Lon", "Lat", "p")
