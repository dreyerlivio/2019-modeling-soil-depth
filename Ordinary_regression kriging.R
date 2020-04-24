#Section 1: Load packages#############################################################################
library(gstat)
library(maptools)
library(rgdal)
library(sp)
library(e1071)
library(raster)
library(ithir)

#Section 2: Read depthfile and log-transform if not normaly distributed###############################
depth <- read.delim("Soil_depth.csv",sep=";", header=TRUE)#read depth file

#Make depth data a spatial Data Frame
coordinates(depth)=~POINT_X+POINT_Y

#Plot histogramm of Depth
hist(depth$Soil_depth)

#Log tranform data so data is normally distributed
depth$Soil_depth.log<-log(depth$Soil_depth)
hist(depth$Soil_depth.log)#Plot histogramm

#Section 3: Read all input rasters and set categorical rasters as factor##############################
#read rasters
slope=read.asciigrid("slope.txt")
geotech=read.asciigrid("geotech.txt")
twi=read.asciigrid("twi.txt")
geology=read.asciigrid("geology.txt")
geomorphon=read.asciigrid("geomorphon.txt")
tpi=read.asciigrid("tpi.txt")
curv_plan=read.asciigrid("curv_plan.txt")
curv_prof=read.asciigrid("curv_prof.txt")
overlandflow=read.asciigrid("overlandflow.txt")
aspect=read.asciigrid("aspect.txt")
dem=read.asciigrid("dem.txt")
streampower=read.asciigrid("streampower.txt")
valleydepth=read.asciigrid("valleydepth.txt")
landform=read.asciigrid("landform.txt")
rockdistance=read.asciigrid("rockdistance.txt")
normalizedheight=read.asciigrid("normalizedheight.txt")

#create factor for rasters containing categorical variables
geotech$geotech.txt<-as.factor(geotech$geotech.txt)
geology$geology.txt<-as.factor(geology$geology.txt)
geomorphon$geomorphon.txt<-as.factor(geomorphon$geomorphon.txt)
landform$landform.txt<-as.factor(landform$landform.txt)

#Section 4: Copy values of auxilary rasters at soil depth location###################################
#Overlay rasters
slope.ov = over(depth,slope)  
geotech.ov = over( depth,geotech)
twi.ov=over(depth,twi)
geology.ov=over(depth,geology)
geomorphon.ov=over(depth,geomorphon)
tpi.ov=over(depth,tpi)
curv_plan.ov=over(depth,curv_plan)
curv_prof.ov=over(depth,curv_prof)
overlandflow.ov=over(depth,overlandflow)
aspect.ov=over(depth,aspect)
dem.ov=over(depth,dem)
streampower.ov=over(depth,streampower)
valleydepth.ov=over(depth,valleydepth)
dem.ov=over(depth,dem)
landform.ov=over(depth,landform)
normalizedheight.ov=over(depth,normalizedheight)
rockdistance.ov=over(depth,rockdistance)

#Copy auxilary varialbles to depth file
depth$slope.txt =slope.ov$slope.txt  
depth$geotech.txt=geotech.ov$geotech.txt
depth$twi.txt=twi.ov$twi.txt
depth$geomorphon.txt=geomorphon.ov$geomorphon.txt
depth$geology.txt=geology.ov$geology.txt
depth$tpi.txt=tpi.ov$tpi.txt
depth$curv_plan.txt=curv_plan.ov$curv_plan.txt
depth$curv_prof.txt=curv_prof.ov$curv_prof.txt
depth$overlandflow.txt=overlandflow.ov$overlandflow.txt
depth$aspect.txt=aspect.ov$aspect.txt
depth$streampower.txt=streampower.ov$streampower.txt
depth$valleydepth.txt=valleydepth.ov$valleydepth.txt
depth$dem.txt=dem.ov$dem.txt
depth$landform.txt=landform.ov$landform.txt
depth$normalizedheight.txt=normalizedheight.ov$normalizedheight.txt
depth$rockdistance.txt=rockdistance.ov$rockdistance.txt


#Section 5:  Create dataframe and add ESPG projection code####################################################
#Create dataframe
depthframe<-as.data.frame(depth)
str(depthframe)

#Set categorigal variables in dataframe
depthframe$geology.txt<-as.factor(depthframe$geology.txt)
depthframe$geomorphon.txt<-as.factor(depthframe$geomorphon.txt)
depthframe$geotech.txt<-as.factor(depthframe$geotech.txt)
depthframe$landform.txt<-as.factor(depthframe$landform.txt )

#Ckeck for NA Values
which(!complete.cases(depthframe))
##integer(0)
depthframe<-depthframe[complete.cases(depthframe),]#remove NA values

#Make depthframe a spatial dataframe and add LV95 ESPG projection
coordinates(depthframe)=~POINT_X+POINT_Y
crs(depthframe) <- "+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=2600000 +y_0=1200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs "

# Section 6: Build the regression kriging model#########################################################
#Fit multiple linear regression model
lm.depth <- lm(Soil_depth.log~geology.txt+geotech.txt+aspect.txt
               +twi.txt*tpi.txt+curv_plan.txt+geomorphon.txt
               +streampower.txt+curv_prof.txt+normalizedheight.txt
               +valleydepth.txt+streampower.txt+dem.txt+rockdistance.txt, as.data.frame(depthframe))
names(lm.depth)
summary(lm.depth)#summary of model

#Plot histogramms of residuals if values are normally distributed
hist(lm.depth$residuals)#Histogramm of residuals

#Section 7: Check correlation between soildepth and coovariates
#Plot selected covariate and soil depth
plot(Soil_depth.log ~ slope.txt, data=depthframe)
legend(x=0.7, y=3.2, legend=c("1","2","3"), pch=1, col=1:3)
abline(lm(Soil_depth.log~slope.txt, as.data.frame(depthframe)))
#Check corelation of selected covariate and soildepth
x <- cor(depthframe$Soil_depth.log,depthframe$slope.txt)
print(x)

# Section 8: Automatically choose predictors using the step function##################################
#Select predictors step wise
sel_predict=step(lm.depth)
summary(sel_predict)
sel_predict$call$formula # display variables used for kriging

#Model diagnostics
par(mfrow=c(1,3))
plot(lm.depth, which=c(1,2,5))#Display model diagnostics

# Section 9: Model the variagramm for ordinary kriging and make prediction####################################################################
#Modeling of the variogramm
vg<-variogram(Soil_depth.log~1,location=depthframe)
plot(vg)#Display variogramm
modvg<-vgm(nugget=0.2,model="Pen",range=150,sill=0.1)#Fit by eye
plot(vg,mode=modvg,main="fitted by eye")# Display variogram fitted by eye
model_vg <- fit.variogram(vg, modvg)#Fit variogram automatically
plot(vg,mode=model_vg,main="fitted by gstat")# Plot fitted variogram

# Depth prediction by ordinary kriging:
depth.ok <- krige(Soil_depth.log~1, depthframe,grids, model=model_vg)#Make ordinary kriging prediction
summary(depth.ok)

#Section 9: Model the variogramm for regression kriging and make prediction###########################
#Modeling of the variogramm
vgresidual<-variogram(sel_predict$call$formula,location=depthframe)
plot(vgresidual)#Display variogram of residuals
modres=vgm(psill =0.4,"Pen", range =180,nugget=0.3)#Fit by eye
plot(vgresidual,mode=modres,main="fitted by eye")# Display variogram fitted by eye
model_residual <- fit.variogram(vgresidual, modres)#Fit variogram automatically
plot(vgresidual,mode=model_residual,main="fitted by gstat")# Plot fitted variogram
#Universal Kriging
depth.uk<-krige(sel_predict$call$formula,depthframe,grids,model=model_residual)#Make regression kriging prediciton
summary(depth.uk)

# Section 10: Backtransform prediction of ordinary and regression kriging##############################
#Backtransform values
depth.ok$pred<-exp(depth.ok$var1.pred)
depth.uk$pred<-exp(depth.uk$var1.pred)
#OLS Kriging
dept.ols<-krige(sel_predict$call$formula,depthframe,grids,model=NULL)

# Section 11: Use cross validation to check model performance##########################################
#Cross Validation
depth.ok.cv <- krige.cv(Soil_depth.log~1, depthframe,grids, model=model_vg)#Perform cross-validation
depth.ok.RMSE<-sqrt(sum(depth.ok.cv$residual^2)/length(depth.ok.cv$residual))#Calculate RMSE
str(depth.ok.RMSE)#Display RMSE
depth.uk.cv <- krige.cv(sel_predict$call$formula,depthframe,grids,model=model_residual)#Perform cross-validation
depth.uk.RMSE<-sqrt(sum(depth.uk.cv$residual^2)/length(depth.uk.cv$residual))#Calculate RMSE
str(depth.uk.RMSE)#Display RMSE

Section 12: #Export maps##############################################################################
write.asciigrid(depth.ok["pred"],"Soil_depth_OK.asc")#Export prediction of depth
write.asciigrid(depth.uks["sim1"],"UKs_kriging.asc")#Export prediction of depth



