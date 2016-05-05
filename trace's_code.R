rm(list = ls())
#setwd("Documents/Spatial Stats/final_proj")
#setwd("Z:/Spatial Stats/final_proj")
load("predictors.Rdata")
library(gstat)
library(fields)
library(maps)
library(mvtnorm)
library(base)
library(sp)
library(geoR)
library(maptools)
library(rgdal)

#------maps for each month & ptype (36 individual plots)------#
#first obtain prior probs based on observed frequency

months = as.numeric(substr(dates,5,6))
pi.snow = matrix(nrow = 551, ncol = 9)
pi.rain = matrix(nrow = 551, ncol = 9)
pi.ip = matrix(nrow = 551, ncol = 9)
pi.fzra = matrix(nrow = 551, ncol = 9)

for(i in 1:551){ 
  #i = station
  station.indicies = which(station.ind == i) #indices of station i
  
  for(j in 1:length(unique(months))){ 
    #j = month
    
    #There are some missing probabilities for May and September, so we use the observations from April and October to record frequency:
    if(j == 5){m = sort(unique(months))[j-1]}else{if(j == 6){m = sort(unique(months))[j+1]}else{m = sort(unique(months))[j]}}  
    
    month.ind = which(months == m) 
    month.indicies = which(date.ind%in%month.ind)
    total.month.station = intersect(month.indicies, station.indicies)
    
    snow.indicies = which(ptype[total.month.station] == "SN")
    rain.indicies = which(ptype[total.month.station] == "RA")
    ip.indicies = which(ptype[total.month.station] == "IP")
    fzra.indicies = which(ptype[total.month.station] == "FZRA")
    
    pi.snow[i,j] = length(snow.indicies)/length(total.month.station)
    pi.rain[i,j] = length(rain.indicies)/length(total.month.station)
    pi.ip[i,j] = length(ip.indicies)/length(total.month.station)
    pi.fzra[i,j] = length(fzra.indicies)/length(total.month.station)
    
  }
}

lon.new = lon-360
#alaska.row=which(lat>50)
#us.row=which(lat<=50)
#xy.alaska=cbind(lon.new[alaska.row],lat[alaska.row])
#colnames(xy.alaska)
#xy.us=cbind(lon.new[us.row],lat[us.row])
#colnames(xy.us)
xy=cbind(lon.new,lat)
colnames(xy)
utm=project(xy,proj="+proj=utm +zone=15:19 ellps=WGS84")
#utm.alaska=project(xy.alaska,proj="+proj=utm +zone=01:0 ellps=WGS84")
#utm.us=project(xy.us,proj="+proj=utm +zone=15:19 ellps=WGS84")
#utm=rbind(utm.us,utm.alaska);dim(utm)

utms=utm/1000
colnames(utms)=c("easting","northing")
#plot(utms[,1],utms[,2],pch=20,col=4);title('Utm Locations')
#dim(utm)
d=as.matrix(dist(utms))
max(d)/3
#par(mfrow=c(1,2))
#plot(lon.lat[,1],lon.lat[,2],pch=20);title('Longitute/Latitude Locations')

#Euclid.utm=as.matrix(dist(utm))					#Euclidean distances between utm coordinates
# GCD.latlon=rdist.earth(xy,miles=FALSE)		#Great circle distances between lat/long coordinates
# max(GCD.latlon)/3

# par(mfrow=c(1,2))
# Euclid.utm=as.matrix(dist(utms))
# plot(Euclid.utm);title('Euclidean Distance')
# plot(GCD.latlon):title('Great Circle Distance')

month.names=c("January","February","March","April","May","September","October","November","December")


no.nug.snow = array()
no.nug.rain = array()
no.nug.ip = array()
no.nug.fzra = array()

val.snow=matrix(0,9,3)
val.rain=matrix(0,9,3)
val.ip=matrix(0,9,3)
val.fzra=matrix(0,9,3)

mod.snow=array()
mod.rain=array()
mod.ip=array()
mod.fzra=array()


mod.snow.k=array()
mod.rain.k=array()
mod.ip.k=array()
mod.fzra.k=array()

for(i in 1:9){
  
#   m = sort(unique(months))[i]
#   pdf(file = paste("figures/Priors/prior", m,".pdf", sep="_"))
#   par(mfrow = c(2,2))
#   quilt.plot(lon.new, lat, pi.snow[,i],zlim = c(0,1),main = paste("Snow Probabilities, in ",month.names[i] , sep = ""))
#   map("world", add = T)
#   map("state", add = T)
#   
#   quilt.plot(lon.new, lat, pi.rain[,i], zlim = c(0,1),main = paste("Rain Probabilities, in ",month.names[i] , sep = ""))
#   map("world", add = T)
#   map("state", add = T)
#   
#   quilt.plot(lon.new, lat, pi.ip[,i], zlim = c(0,0.3),main = paste("Ice Pellet Probabilities, in ",month.names[i] , sep = ""))
#   map("world", add = T)
#   map("state", add = T)
#   
#   quilt.plot(lon.new, lat, pi.fzra[,i], zlim = c(0,0.3),main = paste("Freezing Rain Probabilities, in ",month.names[i], sep = ""))
#   map("world", add = T) 
#   map("state", add = T)
#   dev.off()
#   
  new.h=seq(0,max(d)/3,len=1000)
  
  
  
  ####################### SNOW #######################
  
  #### Estimating Model ####
  probs.snow = as.data.frame(cbind(utms, pi.snow[,i]))
  coordinates(probs.snow) = ~easting+northing
  vg.snow=variogram(pi.snow[,i]~1, data=probs.snow,cutoff=max(d)/3,width=100)
  invals = c(0.15,20)
  nugget=0.001
  chosen=choose.model(vg.snow,invals,nugget);print(chosen)
  mod.snow[i]=chosen[1]
  mod.snow.k[i]=chosen[7]
  no.nug.snow[i]=as.logical(chosen[6])
  
  #### Plotting SMVG Model ####
  #pdf(file = paste("figures/Semivariograms/snow_", m,".pdf", sep=""))
  if(no.nug.snow[i]==TRUE){
      values = c(as.numeric(chosen[3]),as.numeric(chosen[4]))
      #plot.smvg.model(vg.snow,values,model=mod.snow[i],new.h,main = paste("Snow in",month.names[i]))
      val.snow[i,2:3]=values
    }else {
      values = c(as.numeric(chosen[2]),as.numeric(chosen[3]),as.numeric(chosen[4]))
      #plot.smvg.nug.model(vg.snow,values,model=mod.snow[i],new.h, main = paste("Snow in",month.names[i])) #xlim = c(0,50),ylim=c(0,0.18),
      val.snow[i,]=values
    }
  #dev.off()
  
  
  ####################### RAIN ####################### 
  
  #### Estimating Model ####
  probs.rain = as.data.frame(cbind(utms, pi.rain[,i]))
  coordinates(probs.rain) = ~easting+northing#~utms[,1]+utms[,2]
  vg.rain=variogram(pi.rain[,i]~1, data=probs.rain,cutoff=max(d)/3,width=100)
  invals = c(0.15,20);nugget=0.001
  chosen=choose.model(vg.rain,invals,nugget);print(chosen)
  mod.rain[i]=chosen[1]
  mod.rain.k[i]=chosen[7]
  no.nug.rain[i]=as.logical(chosen[6])
  
  #### Plotting SMVG Model ####
  #pdf(file = paste("figures/Semivariograms/rain_", m,".pdf", sep=""))
  if(no.nug.rain[i]==TRUE){
    values = c(as.numeric(chosen[3]),as.numeric(chosen[4]))
    #plot.smvg.model(vg.rain,values,model=mod.rain[i],new.h,main = paste("Rain in ",month.names[i]))
    val.rain[i,2:3]=values  
  }else {
    values = c(as.numeric(chosen[2]),as.numeric(chosen[3]),as.numeric(chosen[4]))
   # plot.smvg.nug.model(vg.rain,values,model=mod.rain[i],new.h,main=paste("Rain in ",month.names[i]))
    val.rain[i,]=values  
  }
  #plot(vg.rain,pch=19,col=1,ylab=expression(paste("Estimated ",gamma(h))), main = paste("Semivariogram for Rain month",m, sep = " "))
  #dev.off()

  
  ####################### ICE PELLETS #######################  
  
  #### Estimating Model ####
  probs.ip = as.data.frame(cbind(utms, pi.ip[,i]))
  coordinates(probs.ip) = ~easting+northing#~utms[,1]+utms[,2]
  vg.ip=variogram(pi.ip[,i]~1, data=probs.ip,cutoff=max(d)/3,width=100)
  invals = c(0.0018,20);nugget=0.001
  chosen=choose.model(vg.ip,invals,nugget);print(chosen)
  mod.ip[i]=chosen[1]
  mod.ip.k[i]=chosen[7]
  no.nug.ip[i]=as.logical(chosen[6])
  
  #### Plotting SMVG Model ####
  #pdf(file = paste("figures/Semivariograms/ip_", m,".pdf", sep=""))
  if(no.nug.ip[i]==TRUE){
    values = c(as.numeric(chosen[3]),as.numeric(chosen[4]))
    #plot.smvg.model(vg.ip,values,model=mod.ip[i],new.h,main = paste("Pellets in ",month.names[i]))
    val.ip[i,2:3]=values
  }else {
    values = c(as.numeric(chosen[2]),as.numeric(chosen[3]),as.numeric(chosen[4]))
    #plot.smvg.nug.model(vg.ip,values,model=mod.ip[i],new.h, main=paste("Pellets in ",month.names[i]))
    val.ip[i,]=values
  }
  #plot(vg.ip,pch=19,col=1,ylab=expression(paste("Estimated ",gamma(h))), main = paste("Semivariogram for Ice Pellets month",m, sep = " "))
  #dev.off()

  
  ####################### FREEZING RAIN #######################  
  
  #### Estimating Model ####
  probs.fzra = as.data.frame(cbind(utms, pi.fzra[,i]))
  coordinates(probs.fzra) = ~easting+northing#~utms[,1] +utms[,2]
  vg.fzra=variogram(pi.fzra[,i]~1, data=probs.fzra,cutoff=max(d)/3,width=100)
  invals = c(0.0018,20);nugget=0.001
  chosen=choose.model(vg.fzra,invals,nugget);print(chosen)
  mod.fzra[i]=chosen[1]
  mod.fzra.k[i]=chosen[7]
  no.nug.fzra[i]=as.logical(chosen[6])
  
  #### Plotting SMVG Model ####
  #pdf(file = paste("figures/Semivariograms/fzra_", m,".pdf", sep=""))
  if(no.nug.fzra[i]==TRUE){
    values = c(as.numeric(chosen[3]),as.numeric(chosen[4]))
    #plot.smvg.model(vg.fzra,values,model=mod.fzra[i],new.h,main = paste("Freeze in ",month.names[i])) #c(0,50),c(0,0.0025),
    val.fzra[i,2:3]=values
  }else {
    values = c(as.numeric(chosen[2]),as.numeric(chosen[3]),as.numeric(chosen[4]))
    #plot.smvg.nug.model(vg.fzra,values,model=mod.fzra[i],new.h,main=paste("Freeze in ",month.names[i])) #c(0,50),c(0,0.0025), 
    val.fzra[i,]=values
  }
  #dev.off()
  
  print(i)
}

#using the smvg model: krig at each station location (leave one out X val)
#Thin plate spline 
#compare krig & thin plate spline
#go through and check probabilities (> 0, sum =1)

#############################################################   
####################### LOOCV Kriging ####################### 
#############################################################

fix.probs = function(probs){
  new = probs
  if(max(probs)>=1){new=probs/max(probs)}
  if(min(probs)<=0){new=probs-min(probs)}
  new = new/sum(new) 
  return(new)
}

#sum(fix.probs(probs))




prob.pred.rain=matrix(0,551,9)
prob.pred.snow=matrix(0,551,9)
prob.pred.ip=matrix(0,551,9)
prob.pred.fzra=matrix(0,551,9)

for(i in 1:9){
  m = sort(unique(months))[i]
  
####################### SNOW #######################
  prb.snow= as.data.frame(cbind(utms,pi.snow[,i]))
  colnames(prb.snow)=c("x","y","probs")
  coordinates(prb.snow)=~x+y
  
  #vgm1=variogram(probs~1,prb.snow)
  vg.snow=variogram(probs~1, data=prb.snow,cutoff=max(d)/3,width=100)
  #### Kriging Using GSTAT package ####
  initl = vgm(psill=0.15,model="Exp",range=20)
  fit.vg.snow=fit.variogram(vg.snow,initl,fit.method=2)
  par(mfrow=c(1,2))
  plot(vg.snow,model=fit.vg.snow)
  #fit.vg
  snow.krig=krige.cv(probs~1,locations=prb.snow,model=fit.vg.snow)$var1.pred
  
  #quilt.plot(cbind(lon.new,lat),snow.krig$var1.pred,main="Snow Kringing Predictions")

  ####################### RAIN #######################
  prb.rain= as.data.frame(cbind(utms,pi.rain[,i]))
  colnames(prb.rain)=c("x","y","probs")
  coordinates(prb.rain)=~x+y
  
  vg.rain=variogram(probs~1, data=prb.rain,cutoff=max(d)/3,width=100)
  #### Kriging Using GSTAT package ####
  initl = vgm(psill=0.15,model="Exp",range=20)
  fit.vg.rain=fit.variogram(vg.rain,initl,fit.method=2)
  plot(vg.rain,model=fit.vg.rain)
  #fit.vg
 rain.krig=krige.cv(probs~1,locations=prb.rain,model=fit.vg.rain)$var1.pred

  #quilt.plot(cbind(lon.new,lat),rain.krig$var1.pred,main="Rain Kringing Predictions")
  
#   ####################### PELLETS #######################
#   prb.ip= as.data.frame(cbind(utms,pi.ip[,i]))
#   colnames(prb.ip)=c("x","y","probs")
#   coordinates(prb.ip)=~x+y
#   
#   vg.ip=variogram(pi.ip[,i]~1, data=probs.ip,cutoff=max(d)/3,width=100)
#   
#   #### Kriging Using GSTAT package ####
#   initl = vgm(psill=0.00008,model="Sph",range=20)
#   fit.vg.ip=fit.variogram(vg.ip,initl,fit.method=2)
#   plot(vg.ip,model=fit.vg)
#   #fit.vg
#   ip.krig=krige.cv(probs~1,locations=prb.ip,model=fit.vg.ip)$var1.pred
  ip.krig=pi.ip[,i]
  #quilt.plot(cbind(lon.new,lat),ip.krig$var1.pred,main="Pellets Kringing Predictions")
  
  
  ####################### FREEZING #######################
#   prb.fzra= as.data.frame(cbind(utms,pi.fzra))
#   colnames(prb.fzra)=c("x","y","probs")
#   coordinates(prb.fzra)=~x+y
#   
#   vg.fzra=variogram(pi.fzra[,i]~1, data=probs.fzra,cutoff=max(d)/3,width=100)
#   
#   #### Kriging Using GSTAT package ####
#   initl = vgm(psill=0.0018,model="Exp",range=20)
#   fit.vg.fzra=fit.variogram(vg.fzra,initl,fit.method=2)
  #fit.vg
  #fzra.krig=krige.cv(probs~1,locations=prb.fzra,model=fit.vg.fzra)$var1.pred
  fzra.krig=pi.fzra[,i]
  #quilt.plot(cbind(lon.new,lat),fzra.krig$var1.pred,main="Freezing Kringing Predictions")
  
  #norm.probs=matrix(0,551,4)
  #normalizing the probs
  norm.probs=apply(cbind(snow.krig,rain.krig,ip.krig,fzra.krig),1,fix.probs)
  
  prob.pred.snow[,i] = norm.probs[1,]
  prob.pred.rain[,i] = norm.probs[2,]
  prob.pred.ip[,i] = norm.probs[3,]
  prob.pred.fzra[,i] = norm.probs[4,]
  
  ###### Plotting SMVG Models ######
  #pdf(file = paste("figures/Semivariograms/snow_", m,".pdf", sep=""))
    plot(vg.snow,model=fit.vg.snow, main=paste("Snow SMVG Model in ", month.names[i]))
  #dev.off()
  
  #pdf(file = paste("figures/Semivariograms/rain_", m,".pdf", sep=""))
    plot(vg.rain,model=fit.vg.rain, main=paste("Rain SMVG Model in ", month.names[i]))
  #dev.off()
  
  #pdf(file = paste("figures/Semivariograms/ip_", m,".pdf", sep=""))
    #plot(vg.ip,model=fit.vg.ip, main=paste("Pellets SMVG Model in ", month.names[i]))
  #dev.off()
  
 # pdf(file = paste("figures/Semivariograms/fzra_", m,".pdf", sep=""))
    #plot(vg.fzra,model=fit.vg.fzra, main=paste("Freezing SMVG Model in ", month.names[i]))
  #dev.off()
  
  
  ###### Plotting Predicted Values ######
  pdf(file = paste("figures/Priors/Kriging_preds_snra", m,".pdf", sep="_"),width = 15,height = 8)
  par(mfrow=c(1,2))
  
  quilt.plot(cbind(lon.new,lat),prob.pred.snow[,i],main=paste("Snow Kringing Predictions in ",month.names[i]))
  map("world", add = T)
  map("state", add = T)
  
  quilt.plot(cbind(lon.new,lat),prob.pred.rain[,i],main=paste("Rain Kringing Predictions in ",month.names[i]))
  map("world", add = T)
  map("state", add = T)
  
#   quilt.plot(cbind(lon.new,lat),  prob.pred.ip[,i],main=paste("Pellets Kringing Predictions in ",month.names[i]))
#   map("world", add = T)
#   map("state", add = T)
  
#   quilt.plot(cbind(lon.new,lat),prob.pred.fzra[,i],main=paste("Freezing Rain Kringing Predictions in ",month.names[i]))
#   map("world", add = T)
#   map("state", add = T)
  
  dev.off()
  
  print(i)
  
}



write.table(prob.pred.snow, file="LOOCV_snow.txt")
write.table(prob.pred.rain, file="LOOCV_rain.txt")
write.table(prob.pred.ip, file="LOOCV_ip.txt")
write.table(prob.pred.fzra, file="LOOCV_fzra.txt")


par(mfrow=c(1,1))
############################################################   
#################### Thin Plate Splines #################### 
############################################################   
tps.pred.snow=matrix(0,551,9)
tps.pred.rain=matrix(0,551,9)
tps.pred.ip  =matrix(0,551,9)
tps.pred.fzra=matrix(0,551,9)

for(i in 1:9){
print(i)
####################### SNOW #######################
spline.snow=Tps(xy,pi.snow[,i],lon.lat=TRUE,miles=FALSE)
tps.snow=predict(spline.snow)

####################### RAIN #######################
spline.rain=Tps(xy,pi.rain[,i],lon.lat=TRUE,miles=FALSE)
tps.rain=predict(spline.rain)

####################### PELLETS #######################
#spline.ip=Tps(xy,pi.ip[,i],lon.lat=TRUE,miles=FALSE)
#tps.ip=predict(spline.ip)
tps.ip = pi.ip[,i]

####################### FREEZE #######################
#spline.fzra=Tps(xy,pi.fzra[,i],lon.lat=TRUE,miles=FALSE)
#tps.fzra=predict(spline.fzra)
tps.fzra=pi.fzra[,i]
#normalizing the probs
norm.probs=apply(cbind(tps.snow,tps.rain,tps.ip,tps.fzra),1,fix.probs)

tps.pred.snow[,i] = norm.probs[1,]
tps.pred.rain[,i] = norm.probs[2,]
tps.pred.ip[,i] = norm.probs[3,]
tps.pred.fzra[,i] = norm.probs[4,]

#pdf(file = paste("figures/TPS/smoothing", m,".pdf", sep="_"))
# par(mfrow=c(1,1))
# surface(spline.snow,main=paste("Snow Spline in ",month.names[i]),zlim=c(0,1))
# map("world", add = T)
# map("state", add = T)
# surface(spline.rain,main=paste("Rain Spline in ",month.names[i]),zlim=c(0,1))
# map("world", add = T)
# map("state", add = T)
#surface(spline.ip,main=paste("Pellets Spline in ",month.names[i]))
#map("world", add = T)
#map("state", add = T)
#surface(spline.fzra,main=paste("Freezing Spline in ",month.names[i]))
#map("world", add = T)
#map("state", add = T)
#dev.off()

#pdf(file = paste("figures/TPS/tps_pred_", m,".pdf", sep=""))
par(mfrow=c(2,2))

quilt.plot(lon.new,lat,tps.pred.snow[,i],main="Snow TPS Predictions",zlim=c(0,1))
map("world", add = T)
map("state", add = T)

quilt.plot(lon.new,lat,tps.pred.rain[,i],main="Rain TPS Predictions",zlim=c(0,1))
map("world", add = T)
map("state", add = T)

quilt.plot(lon.new,lat,tps.pred.ip[,i],main="Pellets (No Predictions)",zlim=c(0,max(tps.pred.ip[,i])))
map("world", add = T)
map("state", add = T)

quilt.plot(lon.new,lat,tps.pred.fzra[,i],main="Freezing (No Predictions)",zlim=c(0,max(tps.pred.fzra[,i])))
map("world", add = T)
map("state", add = T)
#dev.off()

}

write.table(tps.pred.snow, file="TPS_snow.txt")
write.table(tps.pred.rain, file="TPS_rain.txt")
write.table(tps.pred.ip, file="TPS_ip.txt")
write.table(tps.pred.fzra, file="TPS_fzra.txt")


############# MSE For Kriging Data ###############

MSE=matrix(nrow=2,ncol=4)
colnames(MSE)=c("snow","rain","ip","fzra")
rownames(MSE)= c("krige","TPS")

error.snow.k = (pi.snow-prob.pred.snow)^2
MSE.k.snow.mon=apply(error.snow.k,2,mean)
MSE[1,1] = mean(MSE.k.snow.mon)

error.rain.k = (pi.rain-prob.pred.rain)^2
MSE.k.rain.mon=apply(error.rain.k,2,mean)
MSE[1,2] = mean(MSE.k.rain.mon)

error.ip.k = (pi.ip-prob.pred.ip)^2
MSE.k.ip.mon=apply(error.ip.k,2,mean)
MSE[1,3] = mean(MSE.k.ip.mon)

error.fzra.k = (pi.fzra-prob.pred.fzra)^2
MSE.k.fzra.mon=apply(error.fzra.k,2,mean)
MSE[1,4] = mean(MSE.k.fzra.mon)

############# MSE For TPS ###############
error.snow.t = (pi.snow-tps.pred.snow)^2
MSE.t.snow.mon=apply(error.snow.t,2,mean)
MSE[2,1] = mean(MSE.t.snow.mon)

error.rain.t = (pi.rain-tps.pred.rain)^2
MSE.t.rain.mon=apply(error.rain.t,2,mean)
MSE[2,2] = mean(MSE.t.rain.mon)

error.ip.t = (pi.ip-tps.pred.ip)^2
MSE.t.ip.mon=apply(error.ip.t,2,mean)
MSE[2,3]= mean(MSE.t.ip.mon)

error.fzra.t = (pi.fzra-tps.pred.fzra)^2
MSE.t.fzra.mon=apply(error.fzra.t,2,mean)
MSE[2,4] = mean(MSE.t.fzra.mon)

MSE

###could be cool to plot all 4 ptype smvg's on one plot for each month:
#https://stat.ethz.ch/pipermail/r-sig-geo/2009-September/006404.html