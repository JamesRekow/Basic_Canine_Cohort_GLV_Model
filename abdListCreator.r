#  James Rekow

abdListCreator = function(M = 40, N = 20, iStrength = 1, univ = 1, sigmaMax = 0.1,
                          thresholdMult = 10 ^ (-1), maxSteps = 10 ^ 4, tStep = 10 ^ (-2),
                          interSmplMult = 0.01, intTime = 100, lambda = 0, returnParams = FALSE,
                          conGraph = NULL){
  
  #  ARGS: M - number of samples
  #        N - number of bacterial species in each sample
  #        iStrength - interaction strength. In range [0, 1]. Describes how strongly the bacterial species
  #                    within each individual sample interact with each other. Higher value indicates a
  #                    greater degree of interaction
  #        univ - universality. In range [0, 1]. Describes how similar the underlying dynamics governing
  #               the time evolution of the bacterial species are in different samples. At univ = 1 the
  #               dynamics are identical in each sample. This determines the variance in the growth
  #               rate vector and interaction matrix of each sample
  #        sigmaMax - maximum standard deviation of elements in the interaction matrix describing
  #                   the bacterial dynamics of a sample. Every sample in a cohort has the same
  #                   value of sigmaMax (e.g. it is a property of the whole cohort). Universality
  #                   acts like a scaling factor being multiplied by sigmaMax when determing sd of
  #                   elements in the interaction matrix. Although universality also effects the
  #                   variance of the growth rate vectors, while sigmaMax does not
  #        thresholdMult - determines the cutoff at which changes are considered small steps. When
  #                        integrating samples if eucidean length of the change in the abundance
  #                        vector during one times step is less than the threshold value then that
  #                        change is considered a small step. If there are 5 consecutive small steps in
  #                        a row the system is assumed to be at equilibrium, and the integration step is
  #                        terminated
  #        maxSteps - maximum number of steps allowed in integration step. If the integration has not
  #                   been completed after this many steps the current abundance vector is returned
  #        tStep - time step in the integration step
  #        interSmplMult - the fraction of a sample's abundances that get transmitted
  #                        upon contact with another sample. Transmission is bi-directional
  #                        and does not deplete the abundance of the dog from which it is 
  #                        transmitted
  #        intTime - total time during the interaction step. Unrelated to tStep
  #        lambda - rate parameter for the exponentially distributed wait times between interactions
  #        returnParams - if TRUE grList and imatList will also be returned
  #        conGraph - "connectivity graph". Vertices represent samples and edges connect
  #                   samples that can interact
  #
  #  RETURNS: abdList - a list of length M containing the integrated and interacted (if lambda > 0)
  #                     numerical abundance vectors of length N from a single cohort
  #
  #  NOTE: typical non-zero value for lambda is 0.5 - 0.15. Expected value of waiting time is 1 / rate
  
  library(igraph)
  source("cohortCreator.r")
  source("eulerIntegrate.r")
  source("intraCohortInteraction.r")
  
  #  give a warning if conGraph is supplied but there are no interactions
  if(lambda == 0 && (!is.null(conGraph))){
    warning("Changing the cohort connectivity (conGraph) will have no effect if 
            there are no interactions.")
  } #  end if
  
  #  used to determine if the system is near equilibrium
  threshold = {thresholdMult * tStep} ^ 2
  
  #  define integrator function to apply eulerIntegrate with desired parameter vals
  integrator = function(smpl){
    return(eulerIntegrate(smpl, threshold = threshold, maxSteps = maxSteps, tStep = tStep))
  } #  end integrator function
  
  #  create cohort and abundance list
  chrt = cohortCreator(M = M, N = N, iStrength = iStrength, univ = univ, sigmaMax = sigmaMax)
  abdList = lapply(chrt, integrator)
  
  #  simulate intra-cohort interaction if lambda > 0
  if(lambda > 0){
    
    #  update the abundances of each sample in the cohort
    for(i in 1:M){
      chrt[[i]][[1]] = abdList[[i]]
    } #  end for
    
    #  simulate intra-cohort interactions
    abdList = intraCohortInteraction(chrt, M = M, N = N, lambda = lambda, threshold = threshold,
                                     maxSteps = maxSteps, tStep = tStep, intTime = intTime,
                                     interSmplMult= interSmplMult, conGraph = conGraph)
    
  } #  end if
  
  #  if returnParams is TRUE return grList and imatList along with abdList
  if(returnParams){
    
    grList = lapply(chrt, "[[", 2)
    imatList = lapply(chrt, "[[", 3)
    
    return(list(abdList, grList, imatList))
    
  } #  end if returnParams
  
  #  default is to only return abdList
  else{
    return(abdList)
  } #  end else
  
} #  end abdListCreator function
