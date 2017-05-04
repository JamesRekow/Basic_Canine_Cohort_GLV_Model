#  James Rekow

conGraphCreator = function(graphType = c("smallworld", "scaleFree", "gnp", "gnm"), meanDegree, v){
  
  #  ARGS: graphType - type of graph to create
  #        meanDegree - mean degree of graph. Some graph types have a limited number of possible
  #                     mean degrees, at least in the igraph implementations.
  #        v - number of vertices in conGraph
  #
  #  RETURNS: conGraph - connectivity graph describing which samples can interact with each other
  
  ##  Mean degree notes: 
  #  sample_pa = 2 * m, m must be an integer
  #  sample_gnp = p * (n - 1)
  #  sample_gnm = 2 * m / n
  #  sample_smallworld = 2 * nei (if dim = 1, only nei and dim affect mean degree)
  
  #  NOTE: the actual mean degree of the produced conGraph is only approximate, but almost always
  #  within 0.4 of the given meanDegree. There is better agreement for larger graphs.
  
  library(igraph)
  
  #  vector of graph types which currenltly have a method to create them
  supportedGraphTypes = c("smallWorld", "scaleFree", "gnp", "gnm")
  
  #  check that the given graph type is supported
  if(!(graphType %in% supportedGraphTypes)){
    stop("Invalid graphType")
  } #  end if
  
  #  check that the given mean degree is even
  evenMeanDegree = meanDegree %% 2 == 0
  
  #  issue a warning if the given mean degree is not even
  if(!evenMeanDegree){
    warning("In the current implementation smallworld, scalefree, and gnm methods may not work as
            intended if meanDegree is odd.")
  } #  end if
  
  ##  Implement the method corresponding to graphType to create a connectivity graph of the given
  ##  type with the given mean degree and number of vertices
  
  if(graphType == "smallWorld"){
    nei = 0.5 * meanDegree
    conGraph = sample_smallworld(dim = 1, size = v, nei = nei, p = 0.2, loops = FALSE, multiple = FALSE)
  } #  end if
  
  if(graphType == "scaleFree"){
    m = 0.5 * meanDegree
    conGraph = sample_pa(n = v, power = 1, m = m, directed = FALSE)
  } #  end if
  
  if(graphType == "gnp"){
    p = meanDegree / (v - 1)
    conGraph = sample_gnp(n = v, p = p, directed = FALSE, loops = FALSE)
  } #  end if
  
  if(graphType == "gnm"){
    m = meanDegree * v / 2
    conGraph = sample_gnm(n = v, m = m, directed = FALSE, loops = FALSE)
  } #  end if
  
  return(conGraph)
  
} #  end conGraphCreator function
