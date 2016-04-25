krige.loocv = function(smvg, model,formula,locations,data){
  
  for(i in 1:length(data)){
    new.data=data[-i,]
    smvg.krig=krige(formula,locations=coordinates(probs.snow),data=probs.snow,newdata=probs.snow.new,model=fit.vg.snow)
    
    
    
  }
  
  
}