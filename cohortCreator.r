#  James Rekow

cohortCreator = function(M, N, iStrength = 1, univ = 1, sigmaMax = 0.1){
  
  #  ARGS: usual suspects
  #
  #  RETURNS: chrt - A list of M samples (which are lists of length 3). Each sample is of 
  #                  the form list(abundance, growth rate, interaction matrix)
  
  #  generate the attributes of the base GLV model
  #  GR = growth rate, IME = interaction matrix elements
  baseGR = runif(N, min = 0, max = 1)
  baseIME = rnorm(N ^ 2, mean = 0, sd = iStrength * sigmaMax)
  
  sampleCreator = function(x = NULL){
    
    #  ARGS: none
    #
    #  RETURNS: a sample (a list of the form list(abd, gr, imat))
    
    #  create abundance vector
    abd = numeric(N)
    numKept = round(N * runif(1, min = 0.2, max = 1))
    initAbd = runif(numKept, min = 0, max = 1)
    initPresent = sample(1:N, numKept)
    abd[initPresent] = initAbd
    
    #  create growth rate vector
    grMultVec = runif(N, min = univ, max = (2 - univ))
    gr = grMultVec * baseGR
    
    #  create interaction matrix
    imatMultVec = runif(N ^ 2, min = univ, max = (2 - univ))
    imat = matrix(imatMultVec * baseIME, nrow = N, ncol = N)
    diag(imat) = -1
    
    return(list(abd, gr, imat))
  } #  end sampleCreator function
  
  #  create cohort
  chrt = lapply(as.list(1:M), sampleCreator)
  
  return(chrt)
  
} #  end cohortCreator function
