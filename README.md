## GLV_Model

James Rekow
4/14/2017

## Introduction ##

The software contained with this README simulates generalized Lotka-Volterra (GLV) models describing
bacterial populations inside host organisms following the methodology described in "Universality of
human microbial dynamics" by Amir Bashan, Travis E. Gibson, Jonathan Friedman, Vincent J. Carey, Scott
T. Weiss, Elizabeth L. Hohmann & Yang-Yu Liu. There are two differences between their software in Matlab
and the R code contained here. First, the criterion for determining if a sample is at a steady state was
based on having a sufficiently small ratio deltax / x in their model, and is based on having a sufficient
number of consecutive sufficiently small delta x values in this model. Second, this model contains a framework
for modeling interactions between samples in a single cohort that builds on their work but is an independent
creation. This work is primarily in the function intraCohortInteraction.r.

## Data Structures ##

sample: A list with three elements: abundance (abd), growth rate (gr), and interaction matrix (imat).
        For a sample with N bacterial species (for example a dog with 10 species of bacteria in its nose)
        abd and gr are each numerical vectors of length N, while imat is a numerical N x N matrix.

	abd - gr - imat

cohort: A list of samples.

	abd - gr - imat
      	      |
              |
	abd - gr - imat
              |
              |
	abd - gr - imat
              |
              |
      	     ...

doc: A list frame containing the overlap and dissimilarity values of of each pair of samples from a cohort.
     The first vector is the overlap values, the second vector is the dissimilarity values. To take advantage
     of R's plotting internals the overlap and dissimilarity vectors are named "x" and "y" respectively.

## Common Variables ##

Note: This section gives a description of the most common variables found in the functions packaged with this
      README. Variables described here are not listed in the ARGS section of functions in which they appear.
      Instead, they are referred to as the "usual suspects".

M             - number of samples in a cohort
N             - number of bacterial species in a sample
iStrength     - interaction strength. A measure of the interaction strength between bacterial species in a sample,
                in the range [0, 1]. Assumed equal for each sample in a given cohort
univ          - universality, in the range [0, 1]. A measure of how similar the interactions between bacterial species
                within a sample are between different samples in the same cohort. If univ = 1 then all samples have
                identical interaction matrices describing the inter-species interactions of bacteria within the sample
                and all samples have identical growth rates
sigmaMax      - parameter controlling the spread of values in the interaction matrices describing inter-bacterial
                interactions within a sample
thresholdMult - controls the sensitivity of the integration step. Lowering this makes the integration
                more accurate
maxSteps      - the maximum number of integration steps allowed by the eulerIntegrate function before terminating
tStep         - time step. Delta t in the integration performed by eulerIntegrate
intTime       - total time during the interaction step. Unrelated to tStep
interSmplMult - the fraction of a sample's abundances that get transmitted upon contact with another sample 
                (transmission is bi-directional and does not deplete the abundances of the sample from which it
                is transmitted)
lambda        - rate parameter for the exponentially distributed wait times between interactions
conGraph      - connectivity graph. Vertices represent samples and edges connect samples that can interact. This
                is an igraph object

## Standard Names of Intermediate Variables ##

smpl     - sample
chrt     - cohort
abd      - abundance
gr       - growth rate
imat     - interaction matrix
abdList  - list of abundance vectors for each sample in a cohort
grList   - list of growth rate vectors for each sample in a cohort
imatList - list of interaction matrices for each sample in a cohort
conGraph - connectivity graph
doc/DOC  - doc object, as described above
DOCDF    - data frame doc instead of a list
DOCNS    - doc that only contains points in the negative slope (NS) region of the overlap-dissimilarity (OD) plane. That
           is, points with overlap greater than or equal to the largest overlap such that the DOC curve has a negative
           slope for all overlap values greater than or equal to that one
