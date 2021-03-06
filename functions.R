##############################################################
############### Spatial/Practicum Function File ##############
##############################################################

lin=function(t2,h){
  t2*h
}
lin.nug=function(t2,t1,h){
  t1+t2*h
}

pow=function(t1,t2, h){
  t1*h^(t2)
}
pow.nug=function(t1,t2,t3,h){
  t3+t1*h^(t2)
}

white=function(t1,h){
  t1
}

lin.bound=function(t1,t2, h){
  less.than=which(h<=t2)
  more.than=which(h>t1)
  fitted1=(t1/t2)*h[less.than]
  fitted2=(t1)*rep(1,length(more.than))
  return(c(fitted1, fitted2))
}
lin.bound.nug=function(t1,t2,t3,h){
  less.than=which(h<=t2)
  more.than=which(h>t1)
  fitted1=(t3+(t1/t2))*h[less.than]
  fitted2=(t3+t1)*rep(1,length(more.than))
  return(c(fitted1, fitted2))
}

circular=function(t1,t2,h){
  less.than=which(h<=t2)
  more.than=which(h>t1)
  fitted1=t1*(1-(2/pi)/cos(h[less.than]/t2)+(2*h[less.than]/(pi*t2))*sqrt(1-(h[less.than]/t1)^2))
  fitted2=(t1)*rep(1,length(more.than))
  return(c(fitted1, fitted2))
}
circular.nug=function(t1,t2,t3,h){
  less.than=which(h<=t2)
  more.than=which(h>t1)
  fitted1=t3+t1*(1-(2/pi)/cos(h[less.than]/t2)+(2*h[less.than]/(pi*t2))*sqrt(1-(h[less.than]/t1)^2))
  fitted2=(t3+t1)*rep(1,length(more.than))
  return(c(fitted1, fitted2))
}

spherical=function(t2,t3,h){
  less.than=which(h<=t3)
  more.than=which(h>t3)
  fitted1=t2*(1.5*h[less.than]/t3-0.5*(h[less.than]/t3)^3)
  fitted2=(t2)*rep(1,length(more.than))
  return(c(fitted1,fitted2))
}	
spherical.nug=function(t2,t3,t1,h){
  less.than=which(h<=t3)
  more.than=which(h>t3)
  fitted1=t1+t2*(1.5*h[less.than]/t3-0.5*(h[less.than]/t3)^3)
  fitted2=(t1+t2)*rep(1,length(more.than))
  return(c(fitted1,fitted2))
}

rational.quadratic=function(t1,t2,h){
  t1*(h^2/(1+(h^2/t2)))
}
rational.quadratic.nug=function(t1,t2,t3,h){
  t3+t1*(h^2/(1+(h^2/t2)))
}

exponential=function(t1,t2,h){
  t1*(1-exp(-h/t2))
}
exponential.nug=function(t1,t2,t3,h){
  t3+t1*(1-exp(-h/t2))
}

gauss=function(t1,t2,h){
  t1*(1-exp(-(h/t2)^2))
}
gauss.nug=function(t1,t2,t3,h){
  t3+t1*(1-exp(-(h/t2)^2))
}

wave=function(t1,t2,h){
  t1*(1-(t2/h)*sin(h/t2))
}

wave.nug=function(t1,t2,t3,h){
  t3+t1*(1-(t2/h)*sin(h/t2))
}


# Calculates Weighted least squares for optim function
WRSS=function(thetas,func,vg){
  if(func=="linear"){
    #gam.theta=lin(thetas[1],thetas[2],vg[,2])
    gam.theta=lin(thetas[1],vg[,2])
  } else if(func=="power"){
    #gam.theta=pow(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=pow(thetas[1],thetas[2],vg[,2])
  } else if(func=="white"){
    gam.theta=white(thetas[3],vg[,2])
  } else if(func=="linear bound"){
    #gam.theta=lin.bound(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=lin.bound(thetas[1],thetas[2],vg[,2])
  } else if(func=="circular"){
    #gam.theta=circular(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=circular(thetas[1],thetas[2],vg[,2])
  } else if(func=="spherical"){
    #gam.theta=spherical(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=spherical(thetas[1],thetas[2],vg[,2])
  } else if(func=="rational quadratic"){
    #gam.theta=rational.quadratic(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=rational.quadratic(thetas[1],thetas[2],vg[,2])
  } else if(func=="exponential"){
    #gam.theta=exponential(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=exponential(thetas[1],thetas[2],vg[,2])
  } else if(func=="gaussian"){
    #gam.theta=gauss(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=gauss(thetas[1],thetas[2],vg[,2])
  } else if(func=="wave"){
    #gam.theta=wave(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=wave(thetas[1],thetas[2],vg[,2])
  } else {
    print("You mispelt the model name.")
    gam.theta = NULL
  }
  #thetas = parameter vector
  #gam.theta=func(thetas[1],thetas[2],thetas[3],vg[,2])
  return(sum((vg[,1]/(gam.theta^2))*((vg[,3]-gam.theta)^2)))
}

WRSS.nug=function(thetas,func,vg){
  if(func=="linear"){
    #gam.theta=lin(thetas[1],thetas[2],vg[,2])
    gam.theta=lin.nug(thetas[1],thetas[3],vg[,2])
  } else if(func=="power"){
    #gam.theta=pow(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=pow.nug(thetas[1],thetas[2],thetas[3],vg[,2])
  } else if(func=="white"){
    gam.theta=white(thetas[3],vg[,2])
  } else if(func=="linear bound"){
    #gam.theta=lin.bound(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=lin.bound.nug(thetas[1],thetas[2],thetas[3],vg[,2])
  } else if(func=="circular"){
    #gam.theta=circular(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=circular.nug(thetas[1],thetas[2],thetas[3],vg[,2])
  } else if(func=="spherical"){
    #gam.theta=spherical(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=spherical.nug(thetas[1],thetas[2],thetas[3],vg[,2])
  } else if(func=="rational quadratic"){
    #gam.theta=rational.quadratic(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=rational.quadratic.nug(thetas[1],thetas[2],thetas[3],vg[,2])
  } else if(func=="exponential"){
    #gam.theta=exponential(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=exponential.nug(thetas[1],thetas[2],thetas[3],vg[,2])
  } else if(func=="gaussian"){
    #gam.theta=gauss(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=gauss.nug(thetas[1],thetas[2],thetas[3],vg[,2])
  } else if(func=="wave"){
    #gam.theta=wave(thetas[1],thetas[2],thetas[3],vg[,2])
    gam.theta=wave.nug(thetas[1],thetas[2],thetas[3],vg[,2])
  } else {
    print("You mispelt the model name.")
    gam.theta = NULL
  }
  #thetas = parameter vector
  #gam.theta=func(thetas[1],thetas[2],thetas[3],vg[,2])
  return(sum((vg[,1]/(gam.theta^2))*((vg[,3]-gam.theta)^2)))
}


#------Estimate the parameters-----#
#param=optim(intial.values,WRSS,nug=nugget,func="wave",vg=smvg)$par
# params.smvg.model=function(smvg,initial.param,model){
#   param=optim(initial.param,WRSS,func=model,vg=smvg)$par
#   return(param)
# }

#-------BEST MODEL from BIC-----#

#can loop through length(model.list) to choose type
choose.model= function(smvg,initial.values, nug){
  initial.vals=c(initial.values, nug)
  #model.list=c("linear", "power", "white", "linear bound", "circular", "spherical", "rational quadratic", "exponential", "gaussian", "wave")
  model.list=c("linear", "power", "white", "rational quadratic", "exponential", "gaussian", "wave")
  
  no.nug=TRUE
  BIC=array()
  BIC.nug=array()
  params = matrix(0,length(model.list),3)
  params.nug = matrix(0,length(model.list),3)
  K = length(smvg[,1])

for(i in 1:length(model.list)){
  params[i,]=optim(initial.vals,WRSS,func=model.list[i],vg=smvg)$par#, method="L-BFGS-B", lower=c(0,0), upper=c(Inf,Inf))$par
  params.nug[i,]=optim(initial.vals,WRSS.nug,func=model.list[i],vg=smvg)$par#, method="L-BFGS-B", lower=c(0,0), upper=c(Inf,Inf))$par
  #compute BIC here
  BIC[i]= K*log(WRSS(params[i,],model.list[i],smvg)/K)+2*2
  BIC.nug[i]=K*log(WRSS.nug(params[i,],model.list[i],smvg)/K)+2*2
    #then choose the lowest val 
}
#Choose Model with nugget or without?
  min.bic = min(BIC)
  min.bic.nug=min(BIC.nug)
  wit.witout=which.min(c(min.bic,min.bic.nug))
  if(wit.witout==1){
    BIC.choose=BIC
    no.nug=TRUE
  } else {
    BIC.choose=BIC.nug
    no.nug=FALSE
    }
  Best.BIC = which.min(BIC.choose)
  
  model=model.list[Best.BIC]
  param=params[Best.BIC,]
  b.i.c=BIC.choose[Best.BIC]
  
  best.model=as.matrix(c(model,param,b.i.c,no.nug))
  row.names(best.model)= c("model","sill","range","nug","BIC","no nug")
  return(best.model)
}
#choose.model(vg.snow,initial.values, nugget)

# Plot the Semivariogram
plot.smvg.model=function(vg,param,model,new.h,main){
  new=array()
  if(model=="linear"){
    title = paste(main,"Linear Model",sep = " ")
    #new=lin(param[1],param[2],new.h)
    new=lin(param[1],new.h)
  } else if(model=="power"){
    title = paste(main,"Power Model",sep=" ")
    #new=pow(param[1],param[2],param[3],new.h)
    new=pow(param[1],param[2],new.h)
  } else if(model=="white"){
    title = paste(main,"White Noise Model",sep = " ")
    new=white(param[3], new.h)*rep(1,length(new.h))
  } else if(model=="linear bound"){
    title = paste(main,"Linear Bound Model",sep = " ")
    #new=lin.bound(param[1],param[2],param[3],new.h)
    new=lin.bound(param[1],param[2],new.h)
  }else if(model=="circular"){
    title = paste(main,"Circular Model",sep = " ")
    #new=circular(param[1],param[2],param[3],new.h)
    new=circular(param[1],param[2],new.h)
  } else if(model=="spherical"){
    title = paste(main,"Spherical Model",sep = " ")
    #new=spherical(param[1],param[2],param[3],new.h)
    new=spherical(param[1],param[2],new.h)
  } else if(model=="rational quadratic"){
    title = paste(main,"Rational Quadratic Model",sep = " ")
    #new=rational.quadratic(param[1],param[2],param[3],new.h)
    new=rational.quadratic(param[1],param[2],new.h)
  } else if(model=="exponential"){
    title = paste(main,"Exponential Model",sep = " ")
    #new=exponential(param[1],param[2],param[3],new.h)
    new=exponential(param[1],param[2],new.h)
  } else if(model=="gaussian"){
    title = paste(main,"Gaussian Model",sep = " ")
    #new=gauss(param[1],param[2],param[3],new.h)
    new=gauss(param[1],param[2],new.h)
  } else if(model=="wave"){
    title = paste(main,"Wave Model",sep = " ")
   #new=wave(param[1],param[2],param[3],new.h)
    new=wave(param[1],param[2],new.h)
  } else {
    print("You typed an invalid model value")
    title=NULL
  }
  new
  plot(vg[,2],vg[,3],pch=19,col=1,xlab="h",ylab=expression(paste("Estimated ",gamma(h))), main=title)
  lines(new.h,new,lwd=3)
}


# Plot the Semivariogram
plot.smvg.nug.model=function(vg,param,model,new.h,main){
  new=array()
  if(model=="linear"){
    title = paste(main,"Linear Nugget Model",sep = " ")
    #new=lin(param[1],param[2],new.h)
    new=lin.nug(param[1],param[2],new.h)
  } else if(model=="power"){
    title = paste(main,"Power Nugget Model",sep = " ")
    #new=pow(param[1],param[2],param[3],new.h)
    new=pow.nug(param[1],param[2],param[3],new.h)
  } else if(model=="white"){
    title = paste(main,"White Noise Nugget Model",sep = " ")
    new=white(param[3], new.h)*rep(1,length(new.h))
  } else if(model=="linear bound"){
    title = paste(main,"Linear Bound Nugget Model",sep = " ")
    #new=lin.bound(param[1],param[2],param[3],new.h)
    new=lin.bound.nug(param[1],param[2],param[3],new.h)
  }else if(model=="circular"){
    title = paste(main,"Circular Nugget Model",sep = " ")
    #new=circular(param[1],param[2],param[3],new.h)
    new=circular.nug(param[1],param[2],param[3],new.h)
  } else if(model=="spherical"){
    title = paste(main,"Spherical Nugget Model",sep = " ")
    #new=spherical(param[1],param[2],param[3],new.h)
    new=spherical.nug(param[1],param[2],param[3],new.h)
  } else if(model=="rational quadratic"){
    title = paste(main,"Rational Quadratic Nugget Model",sep = " ")
    #new=rational.quadratic(param[1],param[2],param[3],new.h)
    new=rational.quadratic.nug(param[1],param[2],param[3],new.h)
  } else if(model=="exponential"){
    title = paste(main,"Exponential Nugget Model",sep = " ")
    #new=exponential(param[1],param[2],param[3],new.h)
    new=exponential.nug(param[1],param[2],param[3],new.h)
  } else if(model=="gaussian"){
    title = paste(main,"Gaussian Nugget Model",sep = " ")
    #new=gauss(param[1],param[2],param[3],new.h)
    new=gauss.nug(param[1],param[2],param[3],new.h)
  } else if(model=="wave"){
    title = paste(main,"Wave Nugget Model",sep = " ")
    #new=wave(param[1],param[2],param[3],new.h)
    new=wave.nug(param[1],param[2],param[3],new.h)
  } else {
    print("You typed an invalid model value")
    title=NULL
  }
  new
  plot(vg[,2],vg[,3],pch=19,col=1,xlab="h",ylab=expression(paste("Estimated ",gamma(h))), main=title)
  lines(new.h,new,lwd=3)
}
