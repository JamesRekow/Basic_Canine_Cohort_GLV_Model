#  James Rekow

conSynthChrtCreator = function(synthM, conGraph, abdList, int = TRUE, numSample = 100, maxIters = 100,
                               numSynthChrt = 10){
  
  #  ARGS: synthM - number of samples in each synthetic cohort
  #        int - if TRUE, create an interacting synthetic cohort. If FALSE, create a non-interacting one
  #        numSample - argument to be passed to disconnectedSubgraphs function if int == FALSE
  #        maxIters - argument to be passed to disconnectedSubgraphs function if int == FALSE
  #        numSynthChrt - number of synthetic cohorts to create
  #
  #  RETURNS: synthChrt - a list of numSynthChrt sublists. Each sublist is the abundance list for a
  #                       synthetic cohort. A synthetic cohort is a subset of the abundance list of a
  #                       cohort that satisfies a given property (e.g. interacting - samples form
  #                       a connected subgraph, or non-interacting - samples form a subgraph with no edges)
  
  source("connectedSubgraphs.r")
  source("disconnectedSubgraphs.r")
  
  #  if creating an interacting cohort, create a list of vectors of indices that form connected subgraphs
  #  of conGraph
  if(int){
    synthChrtIx = connectedSubgraphs(g = conGraph, v = synthM, numSubgraphs = numSynthChrt)
  } #  end if interacting
  
  #  if creating a non-interacting cohort, create a list of vectors of indices that form subgraphs of
  #  conGraph with no edges (e.g. completely disconnected subgraphs)
  if(!int){
    synthChrtIx = disconnectedSubgraphs(g = conGraph, v = synthM, numSubgraphs = numSynthChrt,
                                        numSample = numSample, maxIters = maxIters)
  } #  end if not interacting
  
  #  function that takes a vector of indices and returns a list of abundance vectors corresponding to
  #  those indices
  abdSelector = function(ixVec) abdList[ixVec]
  
  #  for each vector of indices, select the appropriate list of abundance vectors in order to create the
  #  synthetic cohorts
  synthChrt = lapply(synthChrtIx, abdSelector)
  
  return(synthChrt)
  
} #  end conSynthChrtCreator function
