#Section 1:Load packages#########################################################################################
library(ithir)
library(raster)
library(rgdal)
library(sp)
library(randomForest)

#Section 2:#Read depth file and make it a spatial grid frame##############################################
depth <- read.delim("Soil_depth.csv",sep=";", header=TRUE)
coordinates(depth)=~POINT_X+POINT_Y

#Section 3: Read Coovariates to create rasterstack########################################################
slope=raster("slope.txt")
geotech=as.factor(raster("geotech.txt"))
twi=raster("twi.txt")
geology=as.factor(raster("geology.txt"))
geomorphon=as.factor(raster("geomorphon.txt"))
tpi=raster("tpi.txt")
curv_plan=raster("curv_plan.txt")
curv_prof=raster("curv_prof.txt")
overlandflow=raster("overlandflow.txt")
aspect=raster("aspect.txt")
dem=raster("dem.txt")
streampower=raster("streampower.txt")
valleydepth=raster("valleydepth.txt")
landform=raster("landform.txt")
rockdistance=raster("rockdistance.txt")
normalizedheight=raster("normalizedheight.txt")
midslopeposition=raster("midslopeposition.txt")
tci=raster("tci.txt")
stack=stack(slope,geotech,twi,geology,geomorphon,tpi,
            curv_plan,curv_prof,aspect,
            normalizedheight,rockdistance,dem,streampower,valleydepth,midslopeposition,tci)

#Section 4: Extract covariate value at depth location##################################################
#Extract coovariate data 
depthframe <- extract(stack, depth, sp = 1, method = "simple")
#Make depthframe a dataframe 
depthframe<-as.data.frame(depthframe)
str(depthframe)

#Section 5: Check for NA Values and set categorical variables##########################################
#Ckeck for NA Values in Depthframe
which(!complete.cases(depthframe))
##integer(0)
depthframe<-depthframe[complete.cases(depthframe),]#Remove NA values
#set categorical variables
depthframe$geology<-as.factor(depthframe$geology)
depthframe$geomorphon<-as.factor(depthframe$geomorphon)
depthframe$geotech<-as.factor(depthframe$geotech)

#Section 6: #Build Random Forest Model#################################################################
#Create calibration data
set.seed(123)# Shuffle dataset randomly
training <- sample(nrow(depthframe), 0.7 * nrow(depthframe))# Select 70% of dataset as calibration data

# Create the  model
soilddepth.RF <- randomForest(Soil_depth ~ slope + geotech + twi + geology +
                            geotech + geomorphon + tpi + curv_plan+midslopeposition+
                           aspect+dem+streampower+valleydepth+rockdistance+normalizedheight+tci,
                         data = depthframe[training, ],
                          importance = TRUE, ntree = 500)
#Model summary
print(soildepth.RF)
varImpPlot(soilddepth.RF)# Display covariate importance

# Section 7: Model validation
#Validation of calibration dataset
soildepth.RF.cal <- predict(soildepth.RF, newdata = depthframe[training, ])
goof(observed = depthframe$Soil_depth[training], predicted = soildepth.RF.cal)

#Validation of validation dataset
soilddepth.RF.val <- predict(soildepth.RF, newdata = depthframe[-training, ])
goof(observed = depthframe$Soil_depth[-training], predicted = soildepth.RF.val)

#Section 8: Apply Random Forest model spatially and plot the result####################################
#Apply the RF model spatially
soildepth.RF.map <- predict(stack, soildepth.RF, "soildepth_RF.tif",
                     format = "GTiff", datatype = "FLT4S", overwrite = TRUE)
#Plot the map
plot(soildepth.RF.map)
