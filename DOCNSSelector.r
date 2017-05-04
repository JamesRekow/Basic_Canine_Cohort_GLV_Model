#  James Rekow

DOCNSSelector = function(doc){
  
  #  ARGS: doc - input doc
  #
  #  RETURNS: docNS - doc with points outside the NS region (region of OD plane with
  #                   overlap < xns) removed
  
  source("trimmer.r")
  source("firstDeriv.r")
  source("changePointComputer.r")
  
  #  pass NULL inputs through the function (in case a doc with value NULL was produced)
  if(is.null(doc)){
    return(NULL)
  } #  end if
  
  #  extract overlap and dissimilarity vectors
  over = doc[[1]]
  diss = doc[[2]]
  
  #  grid of overlap values on which derivatives and means will be computed
  x = seq(0, 1, 10 ^ (-4))
  
  #  create a DOC curve fitting the points in the DOC plane
  smoothDOC = lowess(over, diss, f = 0.2)
  
  #  trim redundant values in the x and y elements of the lowess object (e.g. if an (x, y) point
  #  is repeated multiple times, remove duplicates)
  #  needed so we can use approx in the next step
  smoothDOC = trimmer(smoothDOC)
  
  #  approximate the values of the DOC curve on the grid of overlap values x
  smoothDOC = approx(smoothDOC$x, smoothDOC$y, xout = x, method = "linear", rule = 2)
  
  #  compute centered first derivative of the DOC curve
  firstDeriv = firstDeriv(smoothDOC$y, 10 ^ (-4))
  
  #  TEST (remove this after NULL check?)
  if(length(firstDeriv) == 0) print(firstDeriv)
  
  #  compute change point
  changePoint = changePointComputer(firstDeriv)
  
  #  if an error was generated in changePointComputer and NULL was output, carry NULL through this
  #  function
  if(is.null(changePoint)){
    return(NULL)
  } #  end if
  
  #  compute xns, the maximal overlap value such that the DOC curve has a negative slope for all overlap
  #  values greater than or equal to xns
  xns = smoothDOC$x[changePoint]
  
  #  select the indices of the overlap vector whose corresponding elements are in the NS region (greater
  #  than or equal to xns)
  NS = over >= xns
  
  #  select the overlap and dissimilarity values in the NS region
  overNS = over[NS]
  dissNS = diss[NS]
  
  #  create a new data frame with the NS values
  docNS = data.frame(x = overNS, y = dissNS)
  
  return(docNS)
  
} #  end DOCNSSelector function
