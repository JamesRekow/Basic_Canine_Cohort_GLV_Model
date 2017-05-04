#  James Rekow

##  A basic example to show how to chain some common functions together.

library(igraph)
source("abdListCreator.r")
source("DOCProcedure.r")
source("DOCNSSelector.r")


#  create a scale free connectivity graph with 40 vertices
conGraph = sample_pa(n = 40, power = 1, directed = FALSE)

#  create a list of integrated and interacted abundance vectors
#  there are 40 samples (M) and 20 bacterial species in each sample (N)
#  the rate parameter for the exponential wait time is 0.1 (lambda)
#  uses the connectivity graph created above when determining which samples can interact with each other
abdList = abdListCreator(M = 40, N = 20, lambda = 0.1, conGraph = conGraph)

#  create a list of two vectors: overlap values (x) and dissimilarity values (y) by applying DOCProcedure
#  to the created list of abundance vectors
doc = DOCProcedure(abdList = abdList)

#  remove points from doc that are not in the NS region
docNS = DOCNSSelector(doc)
