rm(list = ls())
setwd("Z:/Spatial Stats/final_proj")
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
  
  #set.panel
  #somethng is wrong with this dev.off() thing. this works when we do it one by one, maybe we should just do it that way.
  pdf(file = paste("figures/Semivariograms/snow_", m,".pdf", sep=""))
  probs.snow = as.data.frame(cbind(lon.new, lat, pi.snow[,i]))
  coordinates(probs.snow) = ~lon.new +lat
  vg.snow=variogram(pi.snow[,i]~1, data=probs.snow,cutoff=max(d)/2,width=2)
  plot(vg.snow,pch=19,col=1,ylab=expression(paste("Estimated ",gamma(h))), main = paste("Semivariogram for Snow month",m, sep = " ")) 
  dev.off()
  
  pdf(file = paste("figures/Semivariograms/rain_", m,".pdf", sep=""))
  probs.rain = as.data.frame(cbind(lon.new, lat, pi.rain[,i]))
  coordinates(probs.rain) = ~lon.new +lat
  vg.rain=variogram(pi.rain[,i]~1, data=probs.rain,cutoff=max(d)/2,width=2)
  plot(vg.rain,pch=19,col=1,ylab=expression(paste("Estimated ",gamma(h))), main = paste("Semivariogram for Rain month",m, sep = " "))
  dev.off()
  
  pdf(file = paste("figures/Semivariograms/ip_", m,".pdf", sep=""))
  probs.ip = as.data.frame(cbind(lon.new, lat, pi.ip[,i]))
  coordinates(probs.ip) = ~lon.new +lat
  vg.ip=variogram(pi.ip[,i]~1, data=probs.ip,cutoff=max(d)/2,width=4)
  plot(vg.ip,pch=19,col=1,ylab=expression(paste("Estimated ",gamma(h))), main = paste("Semivariogram for Ice Pellets month",m, sep = " "))
  dev.off()
  
  pdf(file = paste("figures/Semivariograms/fzra_", m,".pdf", sep=""))
  probs.fzra = as.data.frame(cbind(lon.new, lat, pi.fzra[,i]))
  coordinates(probs.fzra) = ~lon.new +lat
  vg.fzra=variogram(pi.fzra[,i]~1, data=probs.fzra,cutoff=max(d)/2,width=4)
  plot(vg.fzra,pch=19,col=1,ylab=expression(paste("Estimated ",gamma(h))), main = paste("Semivariogram for Freezing Rain month",m, sep = " "))
  dev.off()
  
  print(i)
}


i = 1 #just trying this outside the loop

#######Fit a wave model
probs.snow = as.data.frame(cbind(lon.new, lat, pi.snow[,i]))
coordinates(probs.snow) = ~lon.new +lat
vg.snow=variogram(pi.snow[,i]~1, data=probs.snow,cutoff=max(d)/2,width=4)
initial=vgm(psill=2.4,model="Wav",range=11)
fit.vg=fit.variogram(vg.snow,initial,fit.method=2)
plot(vg.snow, fit.vg, col = "black", pch = 19, lwd = 3, ylab=expression(paste("Estimated ",gamma(h))), main = "Wave Fit w/o Nugget")
#so basically the plots aren't showing up ! but we would theoretically do this for every ptype and month

######Building a function
vg = vg.snow
#Wave W/out a Nugget
wave=function(t2,t3,h){
  t2*(1-(t3/h)*sin(h/t3))
}

spherical=function(t2,t3,h){
  less.than=which(h<=t3)
  more.than=which(h>t3)
  fitted1=t2*(1.5*h[less.than]/t3-0.5*(h[less.than]/t3)^3)
  fitted2=(t2)*rep(1,length(more.than))
  return(c(fitted1,fitted2))
}	

WRSS.wave=function(thetas){
  #thetas = parameter vector
  
  gam.theta=wave(thetas[1],thetas[2],vg[,2])
  sum((vg[,1]/(gam.theta^2))*((vg[,3]-gam.theta)^2))
}

WRSS.sph=function(thetas){
  #thetas = parameter vector
  
  gam.theta=spherical(thetas[1],thetas[2],vg[,2])
  sum((vg[,1]/(gam.theta^2))*((vg[,3]-gam.theta)^2))
}

sph.initial=c(.15,20)	
param.sph=optim(sph.initial, WRSS.sph)$par 

new.h=seq(0,50,len=1000)
plot(vg.snow[,2],vg.snow[,3],pch=19,col=1,xlab = "distance", ylab=expression(paste("Estimated ",gamma(h),)),main="Wave Fit",ylim=c(0,0.2), xlim = c(0,50))
new.sph=spherical(param.sph[1],param.sph[2],new.h)
lines(new.h,new.sph,lwd=3)

######Look more at the fit before the sill, other models may work better



#using the smvg model: krig at each station location (leave one out X val)
#Thin plate spline 
#compare krig & thin plate spline
#go through and check probabilities (> 0, sum =1)


###could be cool to plot all 4 ptype smvg's on one plot for each month:
#https://stat.ethz.ch/pipermail/r-sig-geo/2009-September/006404.html