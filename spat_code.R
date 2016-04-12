rm(list = ls())
#setwd("Z:\Spatial_Project_2016\Spatial_Project_2016")
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
    
    m = sort(unique(months))[j]
    month.ind = which(months == m) 
    month.indicies = which(date.ind%in%month.ind)
    total.month.station = intersect(month.indicies, station.indicies)
    
    
    snow.indicies = which(ptype[total.month.station] == "SN")
    rain.indicies = which(ptype[total.month.station] == "RA")
    ip.indicies = which(ptype[total.month.station] == "IP")
    fzra.indicies = which(ptype[total.month.station] == "FZRA")
    if(length(total.month.station)==0){
      pi.snow[i,j] = 0
      pi.rain[i,j] = 0
      pi.ip[i,j] = 0
      pi.fzra[i,j] = 0
      
    }
    else{
    pi.snow[i,j] = (length(snow.indicies))/length(total.month.station)
    pi.rain[i,j] = (length(rain.indicies))/length(total.month.station)
    pi.ip[i,j] = (length(ip.indicies))/length(total.month.station)
    pi.fzra[i,j] = (length(fzra.indicies))/length(total.month.station)
    }
    
    }
}
lon.new = lon-360
xy=cbind(lon.new,lat)
d=as.matrix(dist(xy))
max(d)/3

for(i in 1:9){
  m = sort(unique(months))[i]
  pdf(file = paste("figures/Prior_plot", m,".pdf", sep=""))
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
  
  pdf(file = paste("figures/Semivariograms", m,".pdf", sep=""))
  par(mfrow = c(2,2))
  probs.snow = as.data.frame(cbind(lon.new, lat, pi.snow[,i]))
  coordinates(probs.snow) = ~lon.new +lat 
  vg.snow=variogram(pi.snow[,i]~1, data=probs.snow,cutoff=max(d)/2,width=4)
  plot(vg.snow,pch=19,col=1,ylab=expression(paste("Estimated ",gamma(h)))) #change titles
  
  probs.rain = as.data.frame(cbind(lon.new, lat, pi.rain[,i]))
  coordinates(probs.rain) = ~lon.new +lat 
  vg.rain=variogram(pi.rain[,i]~1, data=probs.rain,cutoff=max(d)/2,width=4)
  plot(vg.rain,pch=19,col=1,ylab=expression(paste("Estimated ",gamma(h))))
  
  probs.ip = as.data.frame(cbind(lon.new, lat, pi.ip[,i]))
  coordinates(probs.ip) = ~lon.new +lat 
  vg.ip=variogram(pi.ip[,i]~1, data=probs.ip,cutoff=max(d)/2,width=4)
  plot(vg.ip,pch=19,col=1,ylab=expression(paste("Estimated ",gamma(h))))
  
  probs.fzra = as.data.frame(cbind(lon.new, lat, pi.fzra[,i]))
  coordinates(probs.fzra) = ~lon.new +lat 
  vg.fzra=variogram(pi.fzra[,i]~1, data=probs.fzra,cutoff=max(d)/2,width=4)
  plot(vg.fzra,pch=19,col=1,ylab=expression(paste("Estimated ",gamma(h))))
  dev.off()
}

#try to build a semvg model for each
#using the smvg model: krig at each station location (leave one out X val)
#Thin plate spline 

#compare krig & thin plate spline

#go through and check probabilities (> 0, sum =1)



