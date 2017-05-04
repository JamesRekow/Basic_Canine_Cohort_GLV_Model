#  James Rekow

changePointComputer = function(x){

  #  ARGS: x - numerical vector
  #
  #  RETURNS: the index of the last nonnegative value of x
  
  n = length(x)
  changePoint = n
  
  #  check if x has a positive length
  posLength = n > 0
  
  #  return NULL if input does not have positive length
  if(!identical(posLength, TRUE)){
    return(NULL)
  } #  end if
  
  #note that this is different than in the reference paper for DOC analysis
  #in the paper they use <0 (strict) but this leads to problems b/c the last points
  #often have first deriv 0 even in downward sloping smoothed lowess curves
  
  #  compute the greatest index (changePoint) such that all elements of x with index greater
  #  than the changePoint are negative
  while(x[changePoint] <= 0){
    
    changePoint = changePoint - 1
    
    #  break out of loop if changePoint reaches 0
    #  this occurs if the entire DOC curve has a negative slope
    if(changePoint == 0){
      break
    } #  end if
    
  } #  end while
  
  #  add 1 to changePoint if it's less than n, because the current value of changePoint
  #  is the maximal index of elements in x that are less than changePoint
  if(changePoint < n){
    changePoint = changePoint + 1
  } #  end if
  
  return(changePoint)
  
} #  end changePointComputer function
