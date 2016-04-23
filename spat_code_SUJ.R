rm(list = ls())
setwd("Documents/Spatial Stats/final_proj")
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
xy=cbind(lon.new,lat)
d=as.matrix(dist(xy))
max(d)/3

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

for(i in 1:9){
  
  m = sort(unique(months))[i]
  pdf(file = paste("figures/Priors/prior", m,".pdf", sep="_"))
  par(mfrow = c(2,2))
  quilt.plot(lon.new, lat, pi.snow[,i],zlim = c(0,1),main = paste("Snow Probabilities, Month",m , sep = ""))
  map("world", add = T)
  map("state", add = T)
  
  quilt.plot(lon.new, lat, pi.rain[,i], zlim = c(0,1),main = paste("Rain Probabilities, Month",m , sep = ""))
  map("world", add = T)
  map("state", add = T)
  
  quilt.plot(lon.new, lat, pi.ip[,i], zlim = c(0,0.3),main = paste("Ice Pellet Probabilities, Month",m , sep = ""))
  map("world", add = T)
  map("state", add = T)
  
  quilt.plot(lon.new, lat, pi.fzra[,i], zlim = c(0,0.3),main = paste("Freezing Rain Probabilities, Month",m , sep = ""))
  map("world", add = T) 
  map("state", add = T)
  dev.off()
  
  new.h=seq(0,50,len=1000)
  pdf(file = paste("figures/Semivariograms/snow_", m,".pdf", sep=""))
  probs.snow = as.data.frame(cbind(lon.new, lat, pi.snow[,i]))
  coordinates(probs.snow) = ~lon.new +lat
  vg.snow=variogram(pi.snow[,i]~1, data=probs.snow,cutoff=max(d)/3,width=1)
  invals = c(0.15,20)
  nugget=0.001
  chosen=choose.model(vg.snow,invals,nugget)
  mod=chosen[1]
  no.nug.snow[i]=as.logical(chosen[6])
  if(no.nug.snow[i]==TRUE){
    values = c(as.numeric(chosen[3]),as.numeric(chosen[4]))
    plot.smvg.model(vg.snow,values,model=mod,new.h,main = "Snow")
    val.snow[i,2:3]=values
  }else {
    values = c(as.numeric(chosen[2]),as.numeric(chosen[3]),as.numeric(chosen[4]))
    plot.smvg.nug.model(vg.snow,values,model=mod,new.h, main="Snow") #xlim = c(0,50),ylim=c(0,0.18),
    val.snow[i,]=values
  }
  #plot(vg.snow,pch=19,col=1,ylab=expression(paste("Estimated ",gamma(h))), main = paste("Semivariogram for Snow month",m, sep = " ")) 
  dev.off()
  
  pdf(file = paste("figures/Semivariograms/rain_", m,".pdf", sep=""))
  probs.rain = as.data.frame(cbind(lon.new, lat, pi.rain[,i]))
  coordinates(probs.rain) = ~lon.new +lat
  vg.rain=variogram(pi.rain[,i]~1, data=probs.rain,cutoff=max(d)/3,width=1)
  invals = c(0.15,20);nugget=0.001
  chosen=choose.model(vg.rain,invals,nugget)
  mod=chosen[1]
  no.nug.rain[i]=as.logical(chosen[6])
  if(no.nug.rain[i]==TRUE){
    values = c(as.numeric(chosen[3]),as.numeric(chosen[4]))
    plot.smvg.model(vg.rain,values,model=mod,new.h,main = "Rain")
    val.rain[i,2:3]=values
  }else {
    values = c(as.numeric(chosen[2]),as.numeric(chosen[3]),as.numeric(chosen[4]))
    plot.smvg.nug.model(vg.rain,values,model=mod,new.h,main="Rain")
    val.rain[i,]=values
  }
  #plot(vg.rain,pch=19,col=1,ylab=expression(paste("Estimated ",gamma(h))), main = paste("Semivariogram for Rain month",m, sep = " "))
  dev.off()
  
  pdf(file = paste("figures/Semivariograms/ip_", m,".pdf", sep=""))
  probs.ip = as.data.frame(cbind(lon.new, lat, pi.ip[,i]))
  coordinates(probs.ip) = ~lon.new +lat
  vg.ip=variogram(pi.ip[,i]~1, data=probs.ip,cutoff=max(d)/3,width=1)
  invals = c(0.15,20);nugget=0.001
  chosen=choose.model(vg.ip,invals,nugget)
  mod=chosen[1]
  no.nug.ip[i]=as.logical(chosen[6])
  if(no.nug.ip[i]==TRUE){
    values = c(as.numeric(chosen[3]),as.numeric(chosen[4]))
    plot.smvg.model(vg.ip,values,model=mod,new.h,main = "Pellets")
    val.ip[i,2:3]=values
  }else {
    values = c(as.numeric(chosen[2]),as.numeric(chosen[3]),as.numeric(chosen[4]))
    plot.smvg.nug.model(vg.ip,values,model=mod,new.h, main="Pellets")
    val.ip[i,]=values
  }
  #plot(vg.ip,pch=19,col=1,ylab=expression(paste("Estimated ",gamma(h))), main = paste("Semivariogram for Ice Pellets month",m, sep = " "))
  dev.off()
  
  pdf(file = paste("figures/Semivariograms/fzra_", m,".pdf", sep=""))
  probs.fzra = as.data.frame(cbind(lon.new, lat, pi.fzra[,i]))
  coordinates(probs.fzra) = ~lon.new +lat
  vg.fzra=variogram(pi.fzra[,i]~1, data=probs.fzra,cutoff=max(d)/3,width=1)
  invals = c(0.15,20);nugget=0.001
  chosen=choose.model(vg.fzra,invals,nugget)
  mod=chosen[1]
  no.nug.fzra[i]=as.logical(chosen[6])
  if(no.nug.fzra[i]==TRUE){
    values = c(as.numeric(chosen[3]),as.numeric(chosen[4]))
    plot.smvg.model(vg.fzra,values,model=mod,new.h,main = "Freeze") #c(0,50),c(0,0.0025),
    val.fzra[i,2:3]=values
  }else {
    values = c(as.numeric(chosen[2]),as.numeric(chosen[3]),as.numeric(chosen[4]))
    plot.smvg.nug.model(vg.fzra,values,model=mod,new.h,main="Freeze") #c(0,50),c(0,0.0025), 
    val.fzra[i,]=values
  }
  dev.off()
  
  print(i)
}

#predictions: 4 matrices, same size as pi.snow...
pred.snow = matrix(nrow= 551, ncol=9)
pred.rain = matrix(nrow= 551, ncol=9)
pred.ip = matrix(nrow= 551, ncol=9)
pred.fzra = matrix(nrow= 551, ncol=9)
i = 1 # do this in the loop for each ptype:
#no nug:
for(j in 1:551){
  j=1
  obs = pi.snow[j,i]
  obs = as.data.frame(c(lon.new[j], lat[j], obs))
  coordinates(obs) = ~lon.new+lat #having trouble with this, might work with euclid distances
  
  pi.snow.j = pi.snow[-(j),i] 
  attach(pi.snow.j)
  
  #Fitted Spherical Variogram
  coordinates(pi.snow.j) = ~lon.new+lat
  initial=vgm(psill=val.snow[i,2],model= mod.snow,range=val.snow[i,3]) 
  fit.vg=fit.variogram(vg,initial,fit.method=2)

#Kriging Predicion for data point i
pred[j, i]=krige(pi.snow.j~1, pi.snow.j,obs, model=fit.vg)$var1.pred
}

#Thin plate spline 
#do this in loop for each ptype
xy = cbind(lon.new,lat) 
spline=Tps(xy,array(pi.snow[,i]),lambda=0.15)
pred.snow = predict(spline)
abs(sum(pred.snow - pi.snow[,i])) #just checking how far off it is


#compare krig & thin plate spline: use predictions in main code
#go through and check probabilities (> 0, sum =1)


###could be cool to plot all 4 ptype smvg's on one plot for each month:
#https://stat.ethz.ch/pipermail/r-sig-geo/2009-September/006404.html