#  James Rekow

eulerIntegrate = function(smpl, threshold = 10 ^ (-6), maxSteps = 10 ^ 4, tStep = 10 ^ (-2)){
  
  #  ARGS: smpl - a sample of the form list(abd, gr, imat)
  #        maxSteps - maximum number of time steps allowed before termination
  #        tStep - length of time covered in each iteration (e.g. delta t)
  #        threshold - cutoff below which a change in the abd vector is considered
  #                    a small step
  #
  #  RETURNS: abd - numeric vector containing the integrated abd vector corresponding to smpl
  
  #  extract list elements
  abd = smpl[[1]]
  gr = smpl[[2]]
  imat = smpl[[3]]
  
  #  stores the number of consecutive small steps
  smallSteps = 0
  
  #  integrate the system of ODEs from the GLV model with parameters from smpl using Euler's method
  for(j in 1:maxSteps){
    
    #  compute delta x
    abdChange = abd * {gr + c(imat %*% abd)}
    deltaabd = tStep * abdChange
    abd = abd + deltaabd
    
    #  ensure nonnegativity of abundances
    abd = 0.5 * {abs(abd) + abd}
    
    #  check if delta x is less than the threshold value (euclidean distance with threshold pre-squared)
    smallStep = sum(abdChange * abdChange) < threshold
    
    #  if delta x is less than the threshold value five times in a row, we assume x has 
    #  reached a steady state and terminate the loop
    if(smallStep){
      smallSteps = smallSteps + 1
      
      if(smallSteps == 5){
        break
      } #  end if
    } #  end if
    
    #  reset small step counter if delta x is greater than the threshold value
    else{
      smallSteps = 0
    } #  end else
    
  } #  end integration loop
  
  #  round (in part to set very small abundances to zero)
  abd = round(abd, 4)
  
  return(abd)
  
} #  end eulerIntegrate function
