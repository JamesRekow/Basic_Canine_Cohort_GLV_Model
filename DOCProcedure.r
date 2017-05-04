#  James Rekow

DOCProcedure = function(abdList){
  
  #  ARGS:  abdList - list of integrated abundances from a cohort
  #
  #  RETURNS:  doc - a list of the form list(overlap, dissimilarity) containing
  #            two numerical vectors
  
  source("universality_functions.r")
  
  M = length(abdList)
  
  #  if abdList is empty, return NULL
  if(M == 0){
    return(NULL)
  } #  end if
  
  #  list the indices of all pairs of samples and count the number of such pairs
  pairs = combn(1:M, 2)
  numPairs = choose(M, 2)
  
  #  preallocate list
  doc = list(x = numeric(numPairs), y = numeric(numPairs))
  
  #compute the overlap and dissimilarity of the samples in the cohort
  for(i in 1:numPairs){
    
    x = abdList[[pairs[, i][1]]]
    y = abdList[[pairs[, i][2]]]

    doc[[1]][i] = overlap(x, y)
    doc[[2]][i] = DrJSD(x, y)
    
  } #  end for
  
  return(doc)
  
} #  end DOCProcedure function
