
fix.probs = function(probs){
  new = probs
  if(max(probs)>=1){new=probs/max(probs)}
  if(min(probs)<=0){new=probs-min(probs)}
  new = new/sum(new) 
  return(new)
}

sum(fix.probs(probs))

