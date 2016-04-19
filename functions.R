##############################################################
############### Spatial/Practicum Function File ##############
##############################################################

lin=function(t1,t2,h,no.nug){
  if(no.nug==T){t2*h}
  if(no.nug==F){t1+t2*h}
}

pow=function(t1,t2,t3,h){
  if(t1==0){t2*h^(t3)}
  if(t1!=0){t1+t2*h^(t3)}
}

lin.bound=function(t1,t2,t3,h){
  if(h==0){0}
  if(h<t3 & h>0){
    if(t1==0){(t2/t3)*h}
    if(t1!=0){t1+(t2/t3)*h}
    }
  if(h>t3){
    if(t1==0){t2}
    if(t1!=0){t1+t2}
    }
}

circular=function(t1,t2,t3,h){
  if(h==0){0}
  if(h<t3 & h>0){
    ifelse(t1==0,t2*(1-(2/pi)/cos(h/t3)+(2*h/(pi*t3))*sqrt(1-(h/t3)^2)),t1+t2*(1-(2/pi)/cos(h/t3)+(2*h/(pi*t3))*sqrt(1-(h/t3)^2)))
  }
  if(h>t3){
    ifelse(t1==0,t2,t1+t2)
    }
  
}

spherical=function(t1,t2,t3,h){
  if(h==0){0}
  if(h<t3 & h>0){
    ifelse(t1==0,t2*((1.5*h/t3)-(0.5*(h/t3)^3)),t1+t2*((1.5*h/t3)-(0.5*(h/t3)^3)))
    }
  if(h>t3){
    ifelse(t1==0,t2,t1+t2)
    }
  
}

rational.quadratic=function(t1,t2,t3,h){
  if(t1==0){t2*(h^2/(1+(h^2/t3)))}
  if(t1!=0){t1+t2*(h^2/(1+(h^2/t3)))}
}

exponential=function(t1,t2,t3,h){
  if(t1==0){t2*(1-exp(-h/t3))}
  if(t1!=0){t1+t2*(1-exp(-h/t3))}
}

gauss=function(t1,t2,t3,h){
  if(t1==0){t2*(1-exp(-(h/t3)^2))}
  if(t1!=0){t1+t2*(1-exp(-(h/t3)^2))}
}

wave=function(t1,t2,t3,h){
  if(t1==0){t2*(1-(t3/h)*sin(h/t3))}
  if(t1!=0){t1+t2*(1-(t3/h)*sin(h/t3))}
}

# Calculates Weighted least squares for optim function
WRSS=function(thetas,func,vg,no.nug){
  if(func=="linear"){
    gam.theta=lin(thetas[1],thetas[2],vg[,2], no.nug==no.nug)
  }
  if(func=="power"){
    gam.theta=pow(thetas[1],thetas[2],thetas[3],vg[,2])
  }
  if(func=="linear bound"){
    gam.theta=lin.bound(thetas[1],thetas[2],thetas[3],vg[,2])
  }
  if(func=="circular"){
    gam.theta=circular(thetas[1],thetas[2],thetas[3],vg[,2])
  }
  if(func=="spherical"){
    gam.theta=spherical(thetas[1],thetas[2],thetas[3],vg[,2])
  }
  if(func=="rational quadratic"){
    gam.theta=rational.quadratic(thetas[1],thetas[2],thetas[3],vg[,2])
  }
  if(func=="exponential"){
    gam.theta=exponential(thetas[1],thetas[2],thetas[3],vg[,2])
  }
  if(func=="gaussian"){
    gam.theta=gauss(thetas[1],thetas[2],thetas[3],vg[,2])
  }
  if(func=="wave"){
    gam.theta=wave(thetas[1],thetas[2],thetas[3],vg[,2])
  }
  #thetas = parameter vector
  #gam.theta=func(thetas[1],thetas[2],thetas[3],vg[,2])
  return(sum((vg[,1]/(gam.theta^2))*((vg[,3]-gam.theta)^2)))
}



#Estimate the parameters
params.smvg.model=function(smvg,initial.param,model,no.nug){
  param=optim (initial.param,WRSS,func=model,vg=smvg, no.nug=no.nug)$par
  return(param)
}

# Plot the Semivariogram
plot.smvg.model=function(vg,param,model,new.h,no.nug,xlimit,ylimit){
  new=array()
  if(model=="linear"){
    title = "Linear Model"
    new=lin(param[1],param[2],new.h, no.nug)
  } else if(model=="power"){
    title = "Power Model"
    new=pow(param[1],param[2],param[3],new.h)
  } else if(model=="linear bound"){
    title = "Linear Bound Model"
    new=lin.bound(param[1],param[2],param[3],new.h)
  }else if(model=="circular"){
    title = "Circular Model"
    new=circular(param[1],param[2],param[3],new.h)
  } else if(model=="spherical"){
    title = "Spherical Model"
    new=spherical(param[1],param[2],param[3],new.h)
  } else if(model=="rational quadratic"){
    title = "Rational Quadratic Model"
    new=rational.quadratic(param[1],param[2],param[3],new.h)
  } else if(model=="exponential"){
    title = "Exponential Model"
    new=exponential(param[1],param[2],param[3],new.h)
  } else if(model=="gaussian"){
    title = "Gaussian Model"
    new=gauss(param[1],param[2],param[3],new.h)
  } else if(model=="wave"){
    title = "Wave Model"
   new=wave(param[1],param[2],param[3],new.h)
  } else {
    print("You typed an invalid model value")
    title=NULL
  }
  new
  plot(vg[,2],vg[,3],pch=19,col=1,xlab="h",ylab=expression(paste("Estimated ",gamma(h))), main=title,xlim=xlimit,ylim=ylimit)
  lines(new.h,new,lwd=3)
}

choose.model= function(smvg,initial.param,no.nug){
  bic.calc = array()
  K = length(smvg[,1])
  model=""
  param=c(0,0,0)
  
  param.lin = params.smvg.model(smvg,initial.param,"linear",no.nug)
  bic.calc[1] = K*log(WRSS(param.lin,"linear",smvg,no.nug)/K)+2*2
  
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
  
  BIC.which=which.min(bic.calc)
  BIC=bic.calc[BIC.which]
  
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
