####################################################################
#	Exercise the machinery for estimating a starting guess
#
# mrm & Rajenki Das: Man Uni, 21 Oct 2021
####################################################################

rm( list=ls() ) # Wipe the slate clean.
source( "StartingGuess.R" )

# Read the test data
load( "TestData.RData" )

# Estimate params and look at the estimates
kmean.guess <- subject.data.to.mcmc.start( 3, traj.list, 2 )
kmean.guess

###########################################################
# Reorder the components so that the grand totals of the
# shape parameters in the transition matrices are
# increasing
###########################################################

# Get the totals of the alphas in the transition matrices
alphaGrandTotal = rep( 0, 3 )
for( j in 1:3)  {
  crnt.mc.params <- kmean.guess$component.params[[j]]
  crnt.trans.mat.alphas <- crnt.mc.params$trans.mat.Dirichlet.params
  alphaGrandTotal[j] <- sum( crnt.trans.mat.alphas )
}

# Find a permutation that puts alphaGrandTotal[] into 
# increasing order
alphaGrandTotal
perm <- order( alphaGrandTotal )
perm

# Make a new starting guess whose components are ordered by 
# increasing grand total of shape parameters
ordered.guess <- kmean.guess
ordered.guess$membership.probs <- ordered.guess$membership.probs[perm]
for( j in 1:3) {
  ordered.guess$component.params[[j]] <- kmean.guess$component.params[[perm[j]]]
}

# Look at the ordered result
ordered.guess
