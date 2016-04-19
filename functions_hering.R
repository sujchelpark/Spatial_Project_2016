##############################################################
############### Spatial/Practicum Function File ##############
##############################################################

#----Fitting with No Nugget-----#
lin=function(t1,t2,h){
  t2*h
}

pow=function(t1,t2, h){
  t1*h^(t2)
}

lin.bound=function(t1,t2, h){
  less.than=which(h<=t2)
  more.than=which(h>t1)
  fitted1=(t1/t2)*h[less.than]
  fitted2=(t1)*rep(1,length(more.than))
  return(c(fitted1, fitted2))
}

circular=function(t1,t2,h){
  less.than=which(h<=t2)
  more.than=which(h>t1)
  fitted1=t1*(1-(2/pi)/cos(h[less.than]/t2)+(2*h[less.than]/(pi*t2))*sqrt(1-(h[less.than]/t1)^2))
  fitted2=(t1)*rep(1,length(more.than))
  return(c(fitted1, fitted2))
}

spherical=function(t2,t3,h){
  less.than=which(h<=t3)
  more.than=which(h>t3)
  fitted1=t2*(1.5*h[less.than]/t3-0.5*(h[less.than]/t3)^3)
  fitted2=(t2)*rep(1,length(more.than))
  return(c(fitted1,fitted2))
}	

rational.quadratic=function(t1,t2,h){
  t1*(h^2/(1+(h^2/t2)))
}

exponential=function(t1,t2,h){
  t1*(1-exp(-h/t2))
}

gauss=function(t1,t2,h){
  t1*(1-exp(-(h/t2)^2))
}

wave=function(t1,t2,h){
  t1*(1-(t2/h)*sin(h/t2))
}

# Calculates Weighted least squares for optim function
WRSS=function(thetas,func,vg){
  if(func=="linear"){
    gam.theta=lin(thetas[1],thetas[2],vg[,2])
  }
  if(func=="power"){
    gam.theta=pow(thetas[1],thetas[2],vg[,2])
  }
  if(func=="linear bound"){
    gam.theta=lin.bound(thetas[1],thetas[2],vg[,2])
  }
  if(func=="circular"){
    gam.theta=circular(thetas[1],thetas[2],vg[,2])
  }
  if(func=="spherical"){
    gam.theta=spherical(thetas[1],thetas[2],vg[,2])
  }
  if(func=="rational quadratic"){
    gam.theta=rational.quadratic(thetas[1],thetas[2],vg[,2])
  }
  if(func=="exponential"){
    gam.theta=exponential(thetas[1],thetas[2],vg[,2])
  }
  if(func=="gaussian"){
    gam.theta=gauss(thetas[1],thetas[2],vg[,2])
  }
  if(func=="wave"){
    gam.theta=wave(thetas[1],thetas[2],vg[,2])
  }
  #thetas = parameter vector
  return(sum((vg[,1]/(gam.theta^2))*((vg[,3]-gam.theta)^2)))
}




#------Estimate the parameters-----#
param=optim(intial.values,WRSS,func="wave",vg=smvg)$par
# params.smvg.model=function(smvg,initial.param,model){
#   param=optim(initial.param,WRSS,func=model,vg=smvg)$par
#   return(param)
# }

#-------BEST MODEL from BIC-----#

#can loop through length(model.list) to choose type
model.list=c("linear", "power", "linear bound", "circular", "spherical", "rational quadratic", "exponential", "gaussian", "wave")

BIC=array()
for(i in 1:length(model.list)){
  
  param=optim(intial.values,WRSS,func=model.list[i],vg=smvg, method="L-BFGS-B", lower=c(0,0), upper=c(Inf,Inf))$par
  #compute BIC here
  BIC[i]=
    
  #then choose the lowest val
}

choose.model= function(smvg,initial.param){
  bic.calc = array()
  K = length(smvg[,1])
  model=""
  param=c(0,0)
  
  param.lin = params.smvg.model(smvg,initial.param,"linear")
  bic.calc[1] = K*log(WRSS(param.lin,"linear",smvg)/K)+2*2
  
  param.pow=params.smvg.model(smvg,initial.param,"power")
  bic.calc[2] = K*log(WRSS(param.pow,"power",smvg)/K)+2*2
  
  param.lin.bou= params.smvg.model(smvg,initial.param,"linear bound")
  bic.calc[3] = K*log(WRSS(param.lin.bou,"linear bound",smvg)/K)+2*2
  
  param.circ=params.smvg.model(smvg,initial.param,"circular")
  bic.calc[4] = K*log(WRSS(param.circ,"circular",smvg)/K)+2*2
  
  param.spher=params.smvg.model(smvg,initial.param,"spherical")
  bic.calc[5] = K*log(WRSS(param.spher,"spherical",smvg)/K)+2*2
  
  param.rat.quad=params.smvg.model(smvg,initial.param,"rational quadratic")
  bic.calc[6] = K*log(WRSS(param.rat.quad,"rational quadratic",smvg)/K)+2*2
  
  param.exp=params.smvg.model(smvg,initial.param,"exponential")
  bic.calc[7] = K*log(WRSS(param.exp,"exponential",smvg)/K)+2*2
  
  param.gaus=params.smvg.model(smvg,initial.param,"gaussian")
  bic.calc[8] = K*log(WRSS(param.gaus,"gaussian",smvg)/K)+2*2
  
  param.wa=params.smvg.model(smvg,initial.param,"wave")
  bic.calc[9] = K*log(WRSS(param.wa,"wave",smvg)/K)+2*2
  
  BIC=which.min(bic.calc)
  
  if(BIC==1){
    model="linear"
    param=param.lin
  }
  if(BIC==2){
    model="power"
    param=param.pow
  }
  if(BIC==3){
    model="linear bound"
    param=param.lin.bou
  }
  if(BIC==4){
    model="circular"
    param=param.circ
  }
  if(BIC==5){
    model="spherical"
    param=param.spher
  }
  if(BIC==6){
    model="rational quadratic"
    param=param.rat.quad
  }
  if(BIC==7){
    model="exponential"
    param=param.exp
  }
  if(BIC==8){
    model="gaussian"
    param=param.gaus
  }
  if(BIC==9){
    model="wave"
    param=param.wa
  }
  best.model=cbind(model,param,BIC)
  return(best.model)
}

#------Plot the Semivariogram------#
plot.smvg.model=function(vg,param,model,new.h,xlimit,ylimit){
  new=array()
  if(model=="linear"){
    title = "Linear Model"
    new=lin(param[1],param[2],new.h)
  }
  if(model=="power"){
    title = "Power Model"
    new=pow(param[1],param[2],param[3],new.h)
  }
  if(model=="linear bound"){
    title = "Linear Bound Model"
    new=lin.bound(param[1],param[2],param[3],new.h)
  }
  if(model=="circular"){
    title = "Circular Model"
    new=circular(param[1],param[2],param[3],new.h)
  }
  if(model=="spherical"){
    title = "Spherical Model"
    new=spherical(param[1],param[2],param[3],new.h)
  }
  if(model=="rational quadratic"){
    title = "Rational Quadratic Model"
    new=rational.quadratic(param[1],param[2],param[3],new.h)
  }
  if(model=="exponential"){
    title = "Exponential Model"
    new=exponential(param[1],param[2],param[3],new.h)
  }
  if(model=="gaussian"){
    title = "Gaussian Model"
    new=gauss(param[1],param[2],param[3],new.h)
  }
  if(model=="wave"){
    title = "Wave Model"
    new=wave(param[1],param[2],param[3],new.h)
  }
  #else{
  #print("You typed an invalid model value")
  #title=NULL
  #}
  new
  plot(vg[,2],vg[,3],pch=19,col=1,xlab="h",ylab=expression(paste("Estimated ",gamma(h))), main=title,xlim=xlimit,ylim=ylimit)
  lines(new.h,new,lwd=3)
}