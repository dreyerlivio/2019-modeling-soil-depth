#Section 1: Import modules and set working directory
from pcraster import *
from pcraster.framework import *
import os
import numpy as np
from numpy import empty
os.chdir("U:/Documents/Model/GIST_Sachseln")# Set your working directory
print(os.getcwd())
#%% Section 2: Import maps#####################################################
dem=readmap("dem.map")
TPI=readmap("TPI.map")
C=readmap("Cs.map")
S=readmap("slope.map")
LU=readmap("geology.map")
depth=readmap("depth.map")
threshold=readmap("threshold.map")
#%% Section 3: Normalize data and create massmovement map
TPI= (TPI-mapmaximum(TPI))/(mapminimum(TPI)-mapmaximum(TPI))#Normalize data
C= (C-mapminimum(C))/(mapmaximum(C)-mapminimum(C))#Normalize data
C= (C-mapmaximum(C))/(mapminimum(C)-mapmaximum(C))#Invert raster
P= (16.381*TPI**4) - (33.232*TPI**3) + (22.261*TPI**2) - (5.4301*TPI) + 0.4522# Insert soildepth-TPI function
P= (P-mapmaximum(P))/(mapminimum(P)-mapmaximum(P))#Normalize data
P= (P-mapmaximum(P))/(mapminimum(P)-mapmaximum(P))#Invert raster
M= ifthenelse(S > threshold,((1+tan(threshold))**-1),1)#Create massmovement map based on threshold


#%% Section 4: Prepare definitions ###########################################
def prepcal(x):
    """ prepares dataset for calibration of the geological units
    x=nominal number of geological unit"""
    u=LU==x
    u=boolean(u)
    t=ifthen(u,scalar(x))
    t=cover(t,0)
    tar=ifthen(u,depth)
    dbol=boolean(tar)
    cal=1
    pred=ifthen(dbol,scalar(cal))
    return pred,tar

def cal():
    """ calibration of soil depth for a geological unit
    """
    target=pcr2numpy(tar,1000)
    b=target<1000
    target=target[b]
    i=i=np.arange(0,100,0.1).tolist()
    o= empty([1,2])
    for t in i:
        predi=pred*M*C*P*t
        prediction=pcr2numpy(predi,1000)
        c=prediction<1000
        prediction=prediction[c]
        o=np.append(o,[(rmse(prediction,target),t)],axis=0)
    o = numpy.delete(o, (0), axis=0)
    o=o[np.argmin(o[:, 0]), 1] 
    print (o)
    return o

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())
#%% Section 5: Calibrate each geological unit#################################
pred,tar=prepcal(1)
geo=cover(scalar(LU==1),0)*cal()
pred,tar=prepcal(2)
geo=cover(scalar(LU==2),0)*cal()+geo
pred,tar=prepcal(3)
geo=cover(scalar(LU==3),0)*cal()+geo
pred,tar=prepcal(4)
geo=cover(scalar(LU==4),0)*cal()+geo
pred,tar=prepcal(6)
geo=cover(scalar(LU==6),0)*cal()+geo
pred,tar=prepcal(7)
geo=cover(scalar(LU==7),0)*cal()+geo
pred,tar=prepcal(8)
geo=cover(scalar(LU==8),0)*cal()+geo
pred,tar=prepcal(11)
geo=cover(scalar(LU==11),0)*cal()+geo
pred,tar=prepcal(13)
geo=cover(scalar(LU==13),0)*cal()+geo
pred,tar=prepcal(14)
geo=cover(scalar(LU==14),0)*cal()+geo
pred,tar=prepcal(16)
geo=cover(scalar(LU==16),0)*cal()+geo
pred,tar=prepcal(22)
geo=cover(scalar(LU==22),0)*cal()+geo
#%% Section 6: Calculate final soil depth, plot it and writte map to disk######
soildepth=geo*M*C*P#Calculate soildepth using the GIST function
soildepth=windowaverage(soildepth,3.0*3.0)# smooth out DEM errors
aguila(soildepth)# Plot the map
report(soildepth,"soildepth_GIST.map")#Writte map to disk

#%% Section 7: Calculate RMSE of soildepth prediction##########################
target=pcr2numpy(depth,1000)
prediction=pcr2numpy(soildepth*(depth*0+1),1000)
b=target<1000
target=target[b]
c=prediction<1000
prediction=prediction[c]
rmse_value=rmse(prediction,target)
print(rmse_value)


