#  James Rekow

trimmer = function(docLowess){
  
  #  ARGS: docLowess - input DOC curve
  #
  #  RETURNS: docLowess with duplicate (x, y) pairs removed

  tmpLowess = docLowess
  copy = diff(tmpLowess$x) == 0
  copy = c(FALSE, copy)
  tmpLowess$x = tmpLowess$x[!copy]
  tmpLowess$y = tmpLowess$y[!copy]
  
  return(tmpLowess)
  
} #  end trimmer function
