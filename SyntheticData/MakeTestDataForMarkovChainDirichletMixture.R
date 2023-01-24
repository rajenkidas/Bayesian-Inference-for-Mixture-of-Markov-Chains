####################################################################
#	Test a circle of functions that fit a mixture model to a set
# of state trajectories. The components of the mixture are *distributions*
# over Markov chains specified by Dirichlet distributions: one per row
# of the transition matrix and one for the initial state.
#
# mrm & Rajenki Das: Whalley Range & Delhi, 15 June 2021
####################################################################

rm( list=ls() ) # Wipe the slate clean.
source( "GenerateMarkovChainDirichletMixtureData.R" )

# Set up a test problem similar to one proposed by Thomas House. It's a mixture 
# with three components, each of which is a distribution over two-state Markov chains.
# The conventions used to specify the mixture are the same as those
# used in GenerateMarkovChainDirichletMixtureData.R

mc.distrib.A <- list(
	transition.matrix.alphas = matrix( c(95, 10, 14, 96), byrow=TRUE, nrow=2 ),
	initial.state.prob.alpha = c( 50, 1 )
)

mc.distrib.B <- list(
  transition.matrix.alphas = matrix( c(55, 15, 45, 45), byrow=TRUE, nrow=2 ),
  initial.state.prob.alpha = c( 1, 50 )
)

mc.distrib.C <- list(
  transition.matrix.alphas = matrix( c(10, 40, 40, 10), byrow=TRUE, nrow=2 ),
  initial.state.prob.alpha = c( 20, 60 )
)

true.DirichletMC.mix <- list( 
	n.components	= 3,
	n.states = 2,
	pi 				= c( 0.2, 0.5, 0.3 ), 
	dirichlet.params	= list( mc.distrib.A, mc.distrib.B, mc.distrib.C ) 	
)

true.DirichletMC.mix # display the test problem

# Sample from the two components separately, to see that 
# the draws look different.
sample.one.chain( mc.distrib.A )
sample.one.chain( mc.distrib.B )

# Generate some random trajectories from the mixture specified above
n.subjects <- 500
target.length <- 500 
synthetic.data <- DirichletMC.mixture.data( n.subjects, true.DirichletMC.mix, target.length )
traj.list <- synthetic.data$state.seqs 

# Spit the trajectories, as well as the params, out to a file
save( true.DirichletMC.mix, traj.list, synthetic.data, file="TestData.RData")


