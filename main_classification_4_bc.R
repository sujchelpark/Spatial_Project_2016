rm(list = ls())
setwd("Z:/Practicum/HW4")

load("predictors.Rdata")
library(gstat)
library(fields)
library(maps)
library(mvtnorm)
library(base)
library(car)
#time: hours = 6.35
#Only looking at the first 16 temperature levels, i.e. Twb.prof[,1:16]


#########################################
#Find the prior probabilities for all data: just using frequency in data
##Split up by ptype, station, and month
#######################################
months = as.numeric(substr(dates,5,6))
pi.snow = matrix(nrow = 551, ncol = 9)
pi.rain = matrix(nrow = 551, ncol = 9)
pi.ip = matrix(nrow = 551, ncol = 9)
pi.fzra = matrix(nrow = 551, ncol = 9)
time <- proc.time()
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
    
    pi.snow[i,j] = (length(snow.indicies))/length(total.month.station)
    pi.rain[i,j] = (length(rain.indicies))/length(total.month.station)
    pi.ip[i,j] = (length(ip.indicies))/length(total.month.station)
    pi.fzra[i,j] = (length(fzra.indicies))/length(total.month.station)
  }
}

#############################
#plot the probabilities for each type and month
#############################
# lon.new = lon-360
# for(i in 1:9){
#   m = sort(unique(months))[i]
#   pdf(file = paste("figures/Prior_plot", m,".pdf", sep=""))
#   par(mfrow = c(2,2))
#   quilt.plot(lon.new, lat, pi.snow[,i],zlim = c(0,1),main = paste("Snow Probabilities, Month",m , sep = ""))
#   map("world", add = T)
#   map("state", add = T)
#   
#   quilt.plot(lon.new, lat, pi.rain[,i], zlim = c(0,1),main = paste("Rain Probabilities, Month",m , sep = ""))
#   map("world", add = T)
#   map("state", add = T)
#   
#   quilt.plot(lon.new, lat, pi.ip[,i], zlim = c(0,0.3),main = paste("Ice Pellet Probabilities, Month",m , sep = ""))
#   map("world", add = T)
#   map("state", add = T)
#   
#   quilt.plot(lon.new, lat, pi.fzra[,i], zlim = c(0,0.3),main = paste("Freezing Rain Probabilities, Month",m , sep = ""))
#   map("world", add = T)
#   map("state", add = T)
#   dev.off()
# }


#################
#Regularize Sigma
#################
ab.BSS = function(param){
  a = param[1]
  b = param[2]
  if(a < 0 | b < 0){return(0)}
  
  sigma.reg.snow[[i]] = a*sigma.snow[[i]] + b*diag(1,16,16)
  sigma.reg.rain[[i]] = a*sigma.rain[[i]] + b*diag(1,16,16)
  sigma.reg.ip[[i]] = a*sigma.ip[[i]] + b*diag(1,16,16)
  sigma.reg.fzra[[i]] = a*sigma.fzra[[i]] + b*diag(1,16,16)
  
  #preallocate space:
  o.snow.train = array(dim = n.train[i])
  o.rain.train = array(dim = n.train[i])
  o.ip.train = array(dim = n.train[i])
  o.fzra.train = array(dim = n.train[i])
  BS.train.indv = array(dim = n.train[i])
  BS.ref.train.indv = array(dim = n.train[i])
  p.snow.train[[i]] = array()
  p.rain.train[[i]] = array()
  p.ip.train[[i]] = array()
  p.fzra.train[[i]] = array()
  
  for(j in 1:n.train[i]){
    
    phi.snow = dmvnorm(snow.train[j,], mean = xbar.snow[,i], sigma = sigma.reg.snow[[i]])
    phi.rain = dmvnorm(train.temp[j,], mean = xbar.rain[,i], sigma = sigma.reg.rain[[i]])
    phi.ip = dmvnorm(ip.train[j,], mean = xbar.ip[,i], sigma = sigma.reg.ip[[i]])
    phi.fzra = dmvnorm(train.temp[j,], mean = xbar.fzra[,i], sigma = sigma.reg.fzra[[i]])
    
    #get station and month for this observation
    obs.ind = obs.indicies.train[j]
    obs.month = months[date.ind[obs.ind]]
    obs.station = stations[station.ind[obs.ind]]
    obs.month.col =  which(sort(unique(months)) == obs.month)
    obs.station.row = which(stations == obs.station)
    
    snow.num = (pi.snow[obs.station.row, obs.month.col]*phi.snow)
    rain.num = (pi.rain[obs.station.row, obs.month.col]*phi.rain)
    ip.num   = (pi.ip[obs.station.row, obs.month.col]*phi.ip)
    fzra.num = (pi.fzra[obs.station.row, obs.month.col]*phi.fzra)
    
    denom = snow.num + rain.num + ip.num + fzra.num 
    
    ###Rename these for training & create them outside of this function
    p.snow.train[[i]][j] = snow.num/denom #pi.snow as list. each element is an array of probs for each obs in testing set i
    p.rain.train[[i]][j] = rain.num/denom
    p.ip.train[[i]][j] = ip.num/denom
    p.fzra.train[[i]][j] = fzra.num/denom
    
    #calculate BS for each testing set, combine at end
    o.snow.train[j] = ifelse(train.type[j] == "SN", 1, 0)
    o.rain.train[j] = ifelse(train.type[j] == "RA", 1, 0)
    o.ip.train[j]   = ifelse(train.type[j] == "IP", 1, 0)
    o.fzra.train[j] = ifelse(train.type[j] == "FZRA", 1, 0)
    
    pi.denom = pi.snow[obs.station.row,obs.month.col] + pi.rain[obs.station.row,obs.month.col] + pi.ip[obs.station.row,obs.month.col] + pi.fzra[obs.station.row,obs.month.col]
    
    BS.ref.train.indv[j] = ((pi.snow[obs.station.row,obs.month.col]/pi.denom) - o.snow.train[j])^2 + ((pi.rain[obs.station.row,obs.month.col]/pi.denom) - o.rain.train[j])^2 + ((pi.ip[obs.station.row,obs.month.col]/pi.denom) - o.ip.train[j])^2 + ((pi.fzra[obs.station.row,obs.month.col]/pi.denom) - o.fzra.train[j])^2
    BS.train.indv[j] =  (p.snow.train[[i]][j] - o.snow.train[j])^2 + (p.rain.train[[i]][j] - o.rain.train[j])^2 + (p.ip.train[[i]][j] - o.ip.train[j])^2 + (p.fzra.train[[i]][j] - o.fzra.train[j])^2 
  }
  BS.ref.train[i] = sum(BS.ref.train.indv, na.rm = T)
  BS.train[i] = sum(BS.train.indv, na.rm=T)
  BSS.train[i] = 1 - (BS.train[i]/BS.ref.train[i])
  return(-BSS.train[i])
}


#Storing the testing and training sets
train.temp = matrix()
train.type = matrix()
test.temp = array()
test.type = array()
n.train = array()
n.test = array()

#16 x 13 matrix where each column is a testing set and each row corresponds to a vertical temperature
xbar.snow = matrix(nrow = 16, ncol = 12)
xbar.rain = matrix(nrow = 16, ncol = 12)
xbar.ip = matrix(nrow = 16, ncol = 12)
xbar.fzra = matrix(nrow = 16, ncol = 12)

#12 16x16 covariance matrices set up as lists
sigma.snow = list()
sigma.rain = list()
sigma.ip = list()
sigma.fzra = list()

sigma.reg.snow = list()
sigma.reg.rain = list()
sigma.reg.ip = list()
sigma.reg.fzra = list()

#probabilities for each observation
p.snow = list()
p.rain = list()
p.ip = list()
p.fzra = list()
p.snow.train = list()
p.rain.train = list()
p.ip.train = list()
p.fzra.train = list()


#For calculating the BS
BS.test = array() 
BS.ref.test = array()
BS.train = array()
BS.ref.train = array()
BSS.train = array()
correct.class = list()
correct.class.tot = array()

short.date = substr(dates, 1, 6)

#For Regularizing Sigma
ab = c(1,10)
ab.all = matrix(0,nrow = 12, ncol = 2)
scale = c(0.1, 30)

# ab.all <- read.table("ab_vals2.txt")
# a <- ab.all$V1
# b <- ab.all$V2

powers.snow=matrix(nrow=16,ncol=12)
powers.rain=matrix(nrow=16,ncol=12)
powers.ip=matrix(nrow=16,ncol=12)
powers.fzra=matrix(nrow=16,ncol=12)

for(i in 1:12){
  #start with year 2001 
  
  
  ####################
  #training data sets:
  date.indicies.train = which(short.date%in%(199609:200105 + i*100 - 100))
  obs.indicies.train=which(date.ind%in%date.indicies.train)
  n.train[i] = length(obs.indicies.train)
  
  train.temp = Twb.prof[obs.indicies.train,1:16]
  train.type = ptype[obs.indicies.train]
  
  #############################
  #means of training data sets
  ##row = vertical distance, col = training set
  ##############################
  snow.train = train.temp[train.type == "SN",]
  pow.snow=apply(snow.train,2,powerTransform)
  for(j in 1:16){
    powers.snow[j,i]=pow.snow[[j]]$lambda
  }
  powers.snow[,i]
  snow.train.norm = bcPower(snow.train,lambda=powers.snow[,i])
  xbar.snow[,i] = apply(snow.train.norm, 2, mean)
  snow.train = apply(train.temp,2,function(x) x^powers.snow[,i])
  
  rain.train = train.temp[train.type == "RA",]
  xbar.rain[,i] = apply(rain.train, 2, mean)
  
  ip.train = train.temp[train.type == "IP",]
  pow.ip=apply(ip.train,2,powerTransform)
  for(j in 1:16){
    powers.ip[j,i]=pow.ip[[j]]$lambda
  }
  powers.ip[,i]
  ip.train.norm = bcPower(ip.train,lambda=powers.ip[,i])
  xbar.ip[,i] = apply(ip.train.norm, 2, mean)
  ip.train = apply(train.temp,2,function(x) x^powers.ip[,i])
  
  fzra.train = train.temp[train.type == "FZRA",]
  xbar.fzra[,i] = apply(fzra.train, 2, mean)
  
  ###################################
  #covariances of training data sets:
  ################################### 
  sigma.snow[[i]] = cov(snow.train.norm) 
  sigma.rain[[i]] = cov(rain.train) 
  sigma.ip[[i]] = cov(ip.train.norm) 
  sigma.fzra[[i]] = cov(fzra.train) 
  
  ##################################
  #Plotting the means and variances:
  ##################################
  # pdf(file = paste("figures/Mean_plot", i,".pdf", sep=""))
  # plot(xbar.snow[,i], 1:16, xlim = c(260, 285), main = paste("Means, Training ", i, sep = ""), xlab = "temperature")
  # points(xbar.rain[,i], 1:16, col=rgb(1,0,0))
  # points(xbar.ip[,i], 1:16, col=rgb(0,1,0))
  # points(xbar.fzra[,i], 1:16, col=rgb(0,0,1))
  # dev.off()
  # 
  # pdf(file = paste("figures/Cov_plot", i,".pdf", sep=""))
  # par(mfrow = c(2,2))
  # image.plot(1:16,1:16,sigma.snow[[i]][,16:1],xaxt = "n", yaxt = "n", main = paste("Snow Reg-Cov, Training ", i, sep = ""))
  # image.plot(1:16,1:16,sigma.rain[[i]][,16:1],xaxt = "n", yaxt = "n", main = paste("Rain Reg-Cov, Training ", i, sep = ""))
  # image.plot(1:16,1:16,sigma.ip[[i]][,16:1],xaxt = "n", yaxt = "n", main = paste("Ice Pellets Reg-Cov, Training ", i, sep = ""))
  # image.plot(1:16,1:16,sigma.fzra[[i]][,16:1],xaxt = "n", yaxt = "n", main = paste("Freezing Rain Reg-Cov, Training ", i, sep = ""))
  # dev.off()
  
  ################### 
  #testing data sets:
  ###################
  date.indicies.test = which(short.date%in%(200109:200205+i*100 -100))
  obs.indicies.test = which(date.ind%in%date.indicies.test)
  n.test[i] = length(obs.indicies.test)
  
  test.temp = Twb.prof[obs.indicies.test,1:16]
  test.type = ptype[obs.indicies.test]
  
  test.snow = apply(test.temp,2,function(x) x^powers.snow[,i])
  test.ip = apply(test.temp,2,function(x) x^powers.ip[,i])
  
  ab = optim(ab, ab.BSS, control = list(parscale = scale))$par
  ab.all[i,] = ab
  
  print(i)
  #############################################
  #Probabilities of each obs being in each type
  #for each testing set
  #############################################
  o.snow = array(dim = n.test[i])
  o.rain = array(dim = n.test[i])
  o.ip = array(dim = n.test[i])
  o.fzra = array(dim = n.test[i])
  
  BS.test.indv = array(dim = n.test[i])
  BS.ref.test.indv = array(dim = n.test[i])
  correct.class[[i]] = array()
  
  p.snow[[i]] = array()
  p.rain[[i]] = array()
  p.ip[[i]] = array()
  p.fzra[[i]] = array()
  
  a = ab.all[i, 1]
  b = ab.all[i, 2]
  sigma.reg.snow[[i]] = a*sigma.snow[[i]] + b*diag(1,16,16)
  sigma.reg.rain[[i]] = a*sigma.rain[[i]] + b*diag(1,16,16)
  sigma.reg.ip[[i]] = a*sigma.ip[[i]] + b*diag(1,16,16)
  sigma.reg.fzra[[i]] = a*sigma.fzra[[i]] + b*diag(1,16,16)
  
  for(j in 1:n.test[i]){
    #j is the indicie of observation in testing set
    
    phi.snow = dmvnorm(test.snow[j,], mean = xbar.snow[,i], sigma = sigma.reg.snow[[i]])
    phi.rain = dmvnorm(test.temp[j,], mean = xbar.rain[,i], sigma = sigma.reg.rain[[i]])
    phi.ip = dmvnorm(test.ip[j,], mean = xbar.ip[,i], sigma = sigma.reg.ip[[i]])
    phi.fzra = dmvnorm(test.temp[j,], mean = xbar.fzra[,i], sigma = sigma.reg.fzra[[i]])
    
    #get station and month for this observation
    obs.ind = obs.indicies.test[j]
    obs.month = months[date.ind[obs.ind]]
    obs.station = stations[station.ind[obs.ind]]
    obs.month.col =  which(sort(unique(months)) == obs.month)
    obs.station.row = which(stations == obs.station)
    
    snow.num = (pi.snow[obs.station.row, obs.month.col]*phi.snow)
    rain.num = (pi.rain[obs.station.row, obs.month.col]*phi.rain)
    ip.num   = (pi.ip[obs.station.row, obs.month.col]*phi.ip)
    fzra.num = (pi.fzra[obs.station.row, obs.month.col]*phi.fzra)
    
    denom = ifelse(snow.num + rain.num + ip.num + fzra.num == 0, 1, snow.num + rain.num + ip.num + fzra.num) 
    
    p.snow[[i]][j] = snow.num/denom #pi.snow as list. each element is an array of probs for each obs in testing set i
    p.rain[[i]][j] = rain.num/denom
    p.ip[[i]][j] = ip.num/denom
    p.fzra[[i]][j] = fzra.num/denom
    
    #calculate BS for each testing set, combine at end
    o.snow[j] = ifelse(test.type[j] == "SN", 1, 0)
    o.rain[j] = ifelse(test.type[j] == "RA", 1, 0)
    o.ip[j]   = ifelse(test.type[j] == "IP", 1, 0)
    o.fzra[j] = ifelse(test.type[j] == "FZRA", 1, 0)
    
    pi.denom = pi.snow[obs.station.row,obs.month.col] + pi.rain[obs.station.row,obs.month.col] + pi.ip[obs.station.row,obs.month.col] + pi.fzra[obs.station.row,obs.month.col]
    
    BS.ref.test.indv[j] = ((pi.snow[obs.station.row,obs.month.col]/pi.denom) - o.snow[j])^2 + ((pi.rain[obs.station.row,obs.month.col]/pi.denom) - o.rain[j])^2 + ((pi.ip[obs.station.row,obs.month.col]/pi.denom) - o.ip[j])^2 + ((pi.fzra[obs.station.row,obs.month.col]/pi.denom) - o.fzra[j])^2
    BS.test.indv[j] =  (p.snow[[i]][j] - o.snow[j])^2 + (p.rain[[i]][j] - o.rain[j])^2 + (p.ip[[i]][j] - o.ip[j])^2 + (p.fzra[[i]][j] - o.fzra[j])^2 
    
    p.vector = c(p.snow[[i]][j], p.rain[[i]][j], p.ip[[i]][j], p.fzra[[i]][j])
    
    # class.freq.snow[j] = ifelse(max(p.type) == "SN", 1, 0) 
    # class.freq.rain[j] = ifelse(max(p.type) == "RA", 1, 0)
    # class.freq.ip[j] = ifelse(max(p.type) == "IP", 1, 0)
    # class.freq.fzra[j] = ifelse(max(p.type) == "FZRA", 1, 0)
    
    if(test.type[j] == "SN"){type = 1}
    if(test.type[j] == "RA"){type = 2}
    if(test.type[j] == "IP"){type = 3}
    if(test.type[j] == "FZRA"){type = 4}
    correct.class[[i]][j] = ifelse(which.max(p.vector) == type, 1, 0)  
    
  }
  BS.ref.test[i] = sum(BS.ref.test.indv, na.rm = T)
  BS.test[i] = sum(BS.test.indv, na.rm=T)
  
  ###############
  #Classification
  ###############
  correct.class.tot[i] = sum(correct.class[[i]][which(correct.class[[i]] == 1)], na.rm = T) 
  # class.freq.tot.snow[i] = sum(class.freq.snow, na.rm = T)
  # class.freq.tot.rain[i] = sum(class.freq.rain, na.rm = T)
  # class.freq.tot.ip[i] = sum(class.freq.ip, na.rm = T)
  # class.freq.tot.fzra[i] = sum(class.freq.fzra, na.rm = T)
}
BS.ref = sum(BS.ref.test, na.rm = T)/sum(n.test) #0.218555
BS = sum(BS.test, na.rm = T)/sum(n.test)#0.0854115
BSS = 1 - (BS/BS.ref) #0.6091995

accuracy = sum(correct.class.tot, na.rm = T)/sum(n.test) #0.942133

# freq.snow = sum(class.freq.tot.snow)/sum(n.test)
# freq.rain = sum(class.freq.tot.rain)/sum(n.test)
# freq.ip = sum(class.freq.tot.ip)/sum(n.test)
# freq.fzra = sum(class.freq.tot.fzra)/sum(n.test)

write.table(ab.all, file = "ab_vals5.txt")

################
#Investigate MVN 
################
#set an i
# pdf(file = paste("figures/Histogram", i,".pdf", sep=""))
# par(mfrow = c(2,2))
# x <- seq(min(snow.train[,1]),max(snow.train[,1]),length = length(snow.train[,1]))
# hist(snow.train[,1], xlab = "temperature", main = paste("Snow at Level 1, Training ", i, sep=""), freq=F, ylim = c(0,0.1))
# curve(dnorm(x, mean = xbar.snow[1,i], sd= sqrt(sigma.snow[[i]][1,1])), col="dark blue",lwd = 3,add=T)
# 
# x <- seq(min(rain.train[,1]),max(rain.train[,1]),length = length(rain.train[,1]))
# hist(rain.train[,1], xlab = "temperature", main = paste("Rain at Level 1, Training ", i, sep=""), freq=F, ylim = c(0,0.1))
# curve(dnorm(x, mean = xbar.rain[1,i], sd= sqrt(sigma.rain[[i]][1,1])), col="dark blue", lwd = 3, add=T)
# 
# x <- seq(min(ip.train[,1]),max(ip.train[,1]),length = length(ip.train[,1]))
# hist(ip.train, xlab = "temperature", main = paste("Ice Pellets at Level 1, Training ", i, sep=""), freq=F, ylim = c(0,0.2))
# curve(dnorm(x, mean = xbar.ip[1,i], sd= sqrt(sigma.ip[[i]][1,1])), col="dark blue", lwd = 3, add=T)
# 
# x <- seq(min(fzra.train[,1]),max(fzra.train[,1]),length = length(fzra.train[,1]))
# hist(fzra.train, xlab = "temperature", main = paste("Freezing Rain at Level 1, Training ", i, sep=""), freq=F, ylim = c(0,0.2))
# curve(dnorm(x, mean = xbar.fzra[1,i], sd= sqrt(sigma.fzra[[i]][1,1])), col="dark blue", lwd = 3, add=T)
# dev.off()
# 
# par(mfrow = c(2,2))
# qqnorm(snow.train[,1], main = "Snow at Level 1, Training 1 Q-Q Plot");qqline(snow.train[,1])
# qqnorm(rain.train[,1], main = "Rain at Level 1, Training 1 Q-Q Plot");qqline(rain.train[,1])
# qqnorm(ip.train[,1], main = "Ice Pellets at Level 1, Training 1 Q-Q Plot");qqline(ip.train[,1])
# qqnorm(fzra.train[,1], main = "Freezing Rain at Level 1, Training 1 Q-Q Plot");qqline(fzra.train[,1])
# 
# par(mfrow=c(2,4))
# plot(snow.train[,1], snow.train[,2], xlab = "first level", ylab = "second level", main = "Snow at Two Levels, Training 1")
# plot(snow.train[,1], snow.train[,16], xlab = "first level", ylab = "sixteenth level", main = "Snow at Two Levels, Training 1")
# 
# plot(rain.train[,1], rain.train[,2], xlab = "first level", ylab = "second level", main = "Rain at Two Levels, Training 1")
# plot(rain.train[,1], rain.train[,16], xlab = "first level", ylab = "sixteenth level", main = "Rain at Two Levels, Training 1")
# 
# plot(ip.train[,1], ip.train[,2], xlab = "first level", ylab = "second level", main = "Ice Pellets at Two Levels, Training 1")
# plot(ip.train[,1], ip.train[,16], xlab = "first level", ylab = "sixteenth level", main = "Ice Pellets at Two Levels, Training 1")
# 
# plot(fzra.train[,1], fzra.train[,2], xlab = "first level", ylab = "second level", main = "Freezing Rain at Two Levels, Training 1")
# plot(fzra.train[,1], fzra.train[,16], xlab = "first level", ylab = "sixteenth level", main = "Freezing Rain at Two Levels, Training 1")

print(proc.time() - time)
