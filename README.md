# variational Bayes expectation maximization of multivariate gaussians for single particle tracking diffusion analysis

This is a mish-mash collection of exploratory MATLAB scripts to employ vbem analysis.  The main script is main_VBEM.m. This script simulates a set of particle tracks with properties specified by the user and then performs variational inference to determine the number of diffusion states and their properties using the expectation-maximization algorithm.  The algorithm will inherently prune excessive number of states and determine the number of diffusion states supported by the data.  

This is a work in progress... Unfortunately, the simulation files is not well organized/coded, albeit it works.  It may be easier to generate the simulations from pEMv2 (https://github.com/p-koo/pEMv2) and modify the main_VBEM.m script to just load the tracks.  

Also, in some cases when there are many diffusive states that have similar properties, vbem will yield the wrong number of diffusive states.  This may be due to too strong priors.  In the end, I opted for pEMv2 instead  (http://arxiv.org/abs/1608.01419), but I still feel like vbem is an elegant approach and can yield good results with better prior tuning.
