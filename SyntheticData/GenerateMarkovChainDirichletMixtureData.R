######################################################
#	Generate data from a mixture of distributions 
# over Markov chains.
#
# mrm & Rajenki Das: Whalley Range & Delhi, 15 June 2020
#
# N.B. The convention here is that, in a transition
# matrix P, P_{ij} is the probability of a 
# transition i -> j.
######################################################

######################################################
#	Construct, at random, a mixture of Markov
# chains. Most of the randomly-generated parameters
# are drawn from suitable Dirichlet distributions.
######################################################
library( gtools ) 

# Check that the params all seem to refer to a chain with the
# same number of states.
validate.Dirichlet.MC <- function( dp ) # dp is short for "Dirichlet params"
{
	verdict <- TRUE ;
	n.states <- length( dp$initial.state.prob.alpha ) ;
	if( (nrow( dp$transition.matrix.alphas ) != n.states ) ||
		(ncol( dp$transition.matrix.alphas ) != n.states )
	) {
		verdict <- FALSE ;
	}

	return( verdict ) ;
}

# Draw a single chain from 
sample.one.chain <- function( dp ) # dp is short for "Dirichlet params"
{
	# Examine the inputs and check for sanity
	n.states <- length( dp$initial.state.prob.alpha ) ;
	if( !validate.Dirichlet.MC(dp) ) {
		stop( "sample.one.Markov.chain: inputs are inconsistently-shaped." ) ;
	}
	
	# Generate transition probabilities
	trans.mat <- matrix( rep(0, n.states*n.states), nrow=n.states )
	for( j in 1:n.states ) {
		trans.mat[j,] <- rdirichlet(1, dp$transition.matrix.alphas[j,] )  ;
	}
	
	# Draw the probs for the initial state
	init.probs <- as.vector( rdirichlet(1, dp$initial.state.prob.alpha) ) ;
	
	# Package up the results and return
	my.markov.params <- list(
		transition.probs = trans.mat,
		initial.state.probs = init.probs
	)
	
	return( my.markov.params )
}

######################################################
#	Generate simulated data.
######################################################

sample.one.state.seq <- function( n.transitions, start.probs, trans.probs ) {	
	# Initialise the result
	state.seq <- rep( 0, n.transitions+1 ) ;
	
	# Choose the initial state
	n.states <- length( start.probs )
	state.seq[1] <- sample( n.states, 1, prob=start.probs ) ;
	
	# Now sample the rest in a Markovian way
	pos <- 1 ;
	while( pos < length(state.seq) ) {
		crnt.state <- state.seq[pos] ;
		next.state <- sample( n.states, 1, prob=trans.probs[crnt.state,] )
		state.seq[pos + 1] <- next.state ;
		pos <- pos + 1 ;
	}
	
	return( state.seq )
}

DirichletMC.mixture.data <- function( n.samples, DMC.mix, target.length=NULL ) 
{
	# Unpack some useful info about the mixture
	n.states <- DMC.mix$n.states ;
	n.comps <- DMC.mix$n.components ;
	
	# If need be, set target.length, then get params for the negative binomial
	# distribution from which we'll draw the lengths.
	if( is.null(target.length) ) {
		target.length <- 200 ; # do something better here
	}
	
	# Choose the desired numbers of observed transitions for the sequences
	# by drawing them from a negative binomial distrib. This is an 
	# arbitrary choice and one could use many other approaches
	nb.size <- 2 ; # Target number of successes for neg. binomial
	n.transitions <- rnbinom( n.samples, nb.size, mu=target.length )
	
	# Choose which components of the mixture to sample from
	my.comp.nums <- sample( n.comps, size=n.samples, replace=TRUE, prob=DMC.mix$pi )
	
	# Sample the Markov chains from the mixture components
	my.chains <- lapply( my.comp.nums,
		function(c) { sample.one.chain( DMC.mix$dirichlet.params[[c]] ) ; } 
	)
	
	# Finally, generate sequences of observed states 
	my.data <- list( 1:n.samples ) # Initialise the result
	for( j in 1:n.samples ) {
		# Get the parameters of the chain from which this sample
		# will be drawn.
		crnt.comp.num <- my.comp.nums[j] ;
		crnt.chain <- my.chains[[j]] ;
		crnt.start.probs <- crnt.chain$initial.state.probs ;
		crnt.trans.probs <- crnt.chain$transition.probs ;
		crnt.n.trans <- n.transitions[j] ;
		
		# Invoke a separate function to do the actual sampling
		my.data[[j]] <- sample.one.state.seq(
			crnt.n.trans, crnt.start.probs, crnt.trans.probs
		) ;
	}
	
	# Add some convenient data, including a list of the
	# components from which the sequences were sampled, to the
	# result.
	result = list( 
		n.comps = n.comps,
		n.states = n.states,
		true.component = my.comp.nums,
		true.chains = my.chains,
		state.seqs = my.data
	) ;
	
	return( result )
}
