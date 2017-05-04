#  James Rekow

##  Contains several simple functions used in the DOC analysis.

overlapSet = function(x, y){
  xPresent = x > 0
  yPresent = y > 0
  return(xPresent & yPresent)
} #  end overlapSet

normalize = function(x){
  return(x / sum(x))
} #  end normalize

#  This function does not assume the data is normalized
overlap = function(x, y){
  
  x = normalize(x)
  y = normalize(y)
  
  bothPresent = overlapSet(x, y)
  return(0.5 * sum(x[bothPresent] + y[bothPresent]))
  
} #  end overlap

#  This function does not assume the data is normalized
DrJSD = function(x, y){
  
  bothPresent = overlapSet(x, y)
  x = normalize(x[bothPresent])
  y = normalize(y[bothPresent])
  m = 0.5 * (x + y)
  
  if(sum(m) == 0){
    return(0)
  }
  
  return(sqrt(0.5 * (sum(x * log(x / m) + y * log(y / m)))))
  
} #  end DrJSD
