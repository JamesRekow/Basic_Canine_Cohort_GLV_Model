#  James Rekow

connectedSubgraphs = function(g, v, numSubgraphs){
  
  #  ARGS: g - graph
  #        v - number of vertices in the connected subgraphs to be returned
  #        numSubgraphs - number of connected subgraphs to return
  #
  #  RETURNS: vList - a list of numSubgraphs number of numeric vectors of vertices forming connected
  #                   subgraphs of g, all of which have v vertices, or NULL if no such subgraphs are 
  #                   found. May contain repeated subgraphs.
  
  library(igraph)
  
  #  get number of vertices in g
  numVertices = vcount(g)
  
  #  create adjacency matrix for connectivity graph
  adjMat = as.matrix(get.adjacency(g, type = "both"))
    
  #  make a list to store samples that are not part of any connected subgraphs with
  #  v or more vertices
  unfitVertices = c()
  
  constructConSubgraph = function(v, numVertices, adjMat, unfitVertices){
    
    #  ARGS: numVertices - number of vertices in created subgraph
    #        adjMat - adjacency matrix of parent graph
    #        unfitVertices - vector of vertices which are not in any connected subgraph with
    #                        at least numVertices number of vertices (e.g. these might be
    #                        a disconnected subgraph with only a few vertices)
    #
    #  RETURNS: vVec - if a viable connected subgraph was found, this is a numeric vector of
    #                  length numVertices containing the indices of vertices which form a
    #                  connected subgraph of the graph corresponding to the adjacency matrix
    #                  adjMat. If the chosen starting vertex is not part of any connected subgraph
    #                  of sufficient size (at least numVertices number of vertices), then multiply
    #                  the first index of vVec by -1 and return all vertices in the connected subgraph
    #                  of which the starting vertex is a part
    
    #  samples in the subgraph
    vVec = 0 * (1:v)
    
    #  potential vertices at which to begin constructing the connected subgraph
    validStartVertices = 1:numVertices
    
    #  if any unfit vertices have been identified, remove those from the vector of valid
    #  staring samples
    if(length(unfitVertices > 0)){
      validStartVertices = validStartVertices[-which(validStartVertices %in% unfitVertices)]
    }
    
    #  select the first vertex
    newVertex = sample(validStartVertices, 1)
    
    #  add the first vertex to the connected subgraph
    vVec[1] = newVertex
    
    #  remove extraneous edges connecting to newVertex from the adjacency matrix
    adjMat[newVertex, ] = 0
    
    #  create a vector of vertices which are connected to any of the samples in the existing
    #  connected subgraph
    potentialVertices = which(adjMat[ , vVec[1]] == 1)
    
    for(ii in 2:v){
      
      #  if there are no potential samples, but the synthetic cohort is not large enough, mark all
      #  of the synthSamples produced on this iteration as unfit, as they form a connected subgraph
      #  that is disconnected from all of the other samples, but insufficiently large to create a
      #  synthetic cohort with synthM samples
      if(length(potentialVertices) == 0){
        
        #  return the negative indices of samples if they are not part of any connected synthetic
        #  cohort with at least synthM vertices
        vVec[1] = -vVec[1]
        break
        
      } #  end if
      
      ##  select a new sample to add to the synthetic cohort
      if(length(potentialVertices) == 1){
        newVertex = potentialVertices[1]
      } #  end if
      
      if(length(potentialVertices) > 1){
        newVertex = sample(potentialVertices, 1) 
      } #  end if
      
      #  remove extraneous edges connecting to the newly added vertex from the adjacency matrix
      adjMat[newVertex, ] = 0
      
      #  remove the newly added vertex from the vector of potential vertices
      potentialVertices = potentialVertices[potentialVertices != newVertex]
      
      #  identify new potential vertices from non-reduntant edges incident to the newly added vertex
      newPotentialVertices = which(adjMat[ , newVertex] == 1)
      
      #  add new potential vertices to potential vertices
      potentialVertices = append(potentialVertices, newPotentialVertices)
      
      #  add the new vertex to vVec
      vVec[ii] = newVertex
      
    } #  end while loop
    
    return(vVec)
    
  } #  end constructConSubgraph function
  
  #  preallocate list of subgraphs
  vList = vector("list", numSubgraphs)
  
  #  create variable to store the number of successfully created connected subgraphs
  countSubgraphs = 0
  
  #  create connected subgraphs until the desired number of connected subgraphs have been obtained,
  #  or it has been determined that there are no connected subgraphs with the desired number of
  #  vertices
  while(length(unfitVertices) < (numVertices - v)){
    
    #  obtain a connected subgraph, or a vector containing unfit vertices
    vVec = constructConSubgraph(v = v, numVertices = numVertices, adjMat = adjMat,
                                unfitVertices = unfitVertices)
    
    #  if vVec contains unfit vertices, then add those vertices to unfitVertices
    if(vVec[1] < 0){
      
      vVec[1] = -vVec[1]
      unfitVertices = append(unfitVertices, vVec)
      
    } #  end if
    
    #  if vVec contains the indices of vertices that comprise a connected subgraph, add it to vList
    else{
      
      vList[[countSubgraphs + 1]] = vVec
      countSubgraphs = countSubgraphs + 1
      
      #  if the desired number of connected subgraphs have been created, break out of the loop
      if(countSubgraphs == numSubgraphs){
        break
      } #  end if
      
    } #  end else
    
  } #  end while loop
  
  #  return NULL and issue a warning if there are no connected subgraphs with the desired number of
  #  vertices
  if(all(unlist(lapply(vList, is.null)))){
    warning("There are no connected subgraphs with the desired number of vertices.")
    vList = NULL
  } #  end if
  
  return(vList)
  
} #  end connectedSubgraphs function
