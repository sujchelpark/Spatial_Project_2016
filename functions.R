##############################################################
############### Spatial/Practicum Function File ##############
##############################################################

lin=function(t1,t2,h){
  if(t1==0){t2*h}
  else{t1+t2*h}
}

pow=function(t1,t2,t3,h){
  if(t1==0){t2*h^(t3)}
  else{t1+t2*h^(t3)}
}

lin.bound=function(t1,t2,t3,h){
  if(h==0){0}
  if(h<t3 & h>0){if(t1==0){(t2/t3)*h}else{t1+(t2/t3)*h}}
  else{if(t1==0){t2}else{t1+t2}}
}

circular=function(t1,t2,t3,h){
  if(h==0){0}
  if(h<t3 & h>0){if(t1==0){t2*(1-(2/pi)*arccos(h/t3)+(2*h/(pi*t3))*sqrt(1-(h/t3)^2))}else{t1+t2*(1-(2/pi)*arccos(h/t3)+(2*h/(pi*t3))*sqrt(1-(h/t3)^2))}}
  else{if(t1==0){t2}else{t1+t2}}
  
}

spherical=function(t1,t2,t3,h){
  if(h==0){0}
  if(h<t3 & h>0){if(t1==0){t2*((1.5*h/t3)-(0.5*(h/t3)^3))}else{t1+t2*((1.5*h/t3)-(0.5*(h/t3)^3))}}
  else{if(t1==0){t2}else{t1+t2}}
  
}

rational.quadratic=function(t1,t2,t3,h){
  if(t1==0){t2*(h^2/(1+(h^2/t3)))}
  else{t1+t2*(h^2/(1+(h^2/t3)))}
}

exponential=function(t1,t2,t3,h){
  if(t1==0){t2*(1-exp(-h/t3))}
  else{t1+t2*(1-exp(-h/t3))}
}

gauss=function(t1,t2,t3,h){
  if(t1==0){t2*(1-exp(-(h/t3)^2))}
  else{t1+t2*(1-exp(-(h/t3)^2))}
}

wave=function(t1,t2,t3,h){
  if(t1==0){t2*(1-(t3/h)*sin(h/t3))}
  else{t1+t2*(1-(t3/h)*sin(h/t3))}
}

# Calculates Weighted least squares for optim function
WRSS=function(thetas){
  #thetas = parameter vector
  gam.theta=func(thetas[1],thetas[2],thetas[3],vg[,2])
  sum((vg[,1]/(gam.theta^2))*((vg[,3]-gam.theta)^2))
}



#Estimate the parameters
params.smvg.model=function(smvg,inital.param,model){
  vg=smvg
  func=model
  param=optim(initial.param,WRSS)$par
  return(param)
}

# Plot the Semivariogram
plot.smvg.model=function(vg,param,model,new.h,title,xlimit,ylimit){
  func=model
  
  plot(vg[,2],vg[,3],pch=19,col=1,xlab="h",ylab=expression(paste("Estimated ",gamma(h))), main=title,xlim=xlimit,ylim=ylimit)
  new.exp=func(param.exp[1],param.exp[2],new.h)
  lines(new.h,new.exp,lwd=3)
}
