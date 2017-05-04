#  James Rekow

intraCohortInteraction = function(chrt, M, N, lambda = 0, threshold = 10 ^ (-6),
                                  maxSteps = 10 ^ 4, tStep = 10 ^ (-2), intTime = 100,
                                  interSmplMult= 0.01, conGraph = NULL){
  
  #  ARGS:  chrt - input cohort in the form list(list(abd, gr, imat)), where abundances have
  #                already been integrated
  #         lambda - rate of exponential random variable describing time between interactions
  #         intTime - time over which interactions are to be simulated
  #         interSmplMult - the fraction of a sample's abundances that get transmitted
  #                         upon contact with another sample (transmission is bi-directional
  #                         and does not deplete the abundance of the dog from which it is 
  #                         transmitted)
  #         conGraph - if a graph describing the connectivity of samples in the cohort is
  #                    supplied, then samples can interact if and only if their corresponding
  #                    vertices are connected by an edge
  #
  #  RETURNS:  abdList - list of interacted and re-integrated abundances from chrt
  #
  #  NOTE: exactly one interaction occurs during each iteration of the loop. If cumTime is less
  #        than intTime but the next wait time makes cumTime greater than intTime then the interaction
  #        corresponding to that wait time does not occur.
  
  library(igraph)
  source("eulerIntegrate.r")

  ##  store the indices of all interacting pairs of samples and count how many such pairs there are
  
  #  if conGraph is not specified, all samples interact with all other samples
  #  (e.g. they form a complete graph)
  if(is.null(conGraph)){
    
    smplPairs = combn(M, 2)
    numPairs = choose(M, 2)
    
  } #  end if
  
  #  if conGraph is specified, use edges in conGraph for smplPairs
  if(!is.null(conGraph)){
    
    smplPairs = t(ends(conGraph, E(conGraph)))
    numPairs = ncol(smplPairs)
    
  } #  end if
  
  #  precompute coefficient for use in interaction step
  tempSmplMult = 1 - interSmplMult
  
  #  track cumulative time and terminate interaction step once cumulative time exceeds intTime
  cumTime = rexp(n = 1, rate = lambda)
  
  #  interact samples while the cumulative time is less than the interaction time
  while(cumTime < intTime){
    
    iPairIx = sample(numPairs, 1)
    iPair = smplPairs[ , iPairIx]
    
    #  select the abundance vectors of the interacting samples
    abd1 = chrt[[iPair[1]]][[1]]
    abd2 = chrt[[iPair[2]]][[1]]
    
    #  compute the weighted sum
    tempSum = interSmplMult * {abd1 + abd2}
    
    #  update sample abundances in the cohort
    chrt[[iPair[1]]][[1]] = tempSmplMult * abd1 + tempSum
    chrt[[iPair[2]]][[1]] = tempSmplMult * abd2 + tempSum
    
    #  re-integrate samples that have interacted
    for(ii in iPair){
      chrt[[ii]][[1]] = eulerIntegrate(chrt[[ii]], threshold = threshold, maxSteps = maxSteps,
                                      tStep = tStep)
    } #  end for
    
    #  sample the exponential wait time again and add this new wait time to the cumulative time
    cumTime = cumTime + rexp(n = 1, rate = lambda)
    
  } #  end while loop
  
  #  select the abundance vectors from the cohort
  abdList = lapply(chrt, "[[", 1)
  
  return(abdList)
  
} #  end intraCohortInteraction function
