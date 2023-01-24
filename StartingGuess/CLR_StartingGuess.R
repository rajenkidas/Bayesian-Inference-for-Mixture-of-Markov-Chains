######################################################################
#	Here is a circle of routines to take state trajectories,
# one per subject, and use them to get starting guesses for
# an MCMC simulation that underpins a Bayesian estimation of
# the parameters for a mixture of Markov chains.
#
# In outline, we:
#
#	(0) Start with a collection of state trajectories, one per subject.
#		Say there are n states and S subjects.
#	(1) Convert each state trajectory into a matrix of counts of
#		observed transitions
#	(2) Send each matrix of counts to a Euclidean space of
#		dimension (n*n) by:
#		(2a) Perfomring the centred log-ratio transform on each row,
#			 yielding n rows, each of which lies in R^n
#		(2b) Concatenating the rows to get a vector in R^(n*n)
#	(3) Doing k-means on the resulting set of S vectors.
#	(4) For each cluster found by k-means:
#		(4a) Estimate the parameters of a Dirichlet distribution
#			 for the per-subject distribution over initial states.
#		(4b) Estimating the parameters of a Dirichlet distribution
#			 for each row of the transition matrix.
#
# mrm & Rajenki Das:	Whalley Range & Delhi, 15 July 2021
######################################################################

library( MCMCprecision ) # fit_dirichlet()

######################################################################
# Take in count data from a subject and send it to CLR space
######################################################################

count.vec.to.CLR <- function( cv ) {
	# Do the centred log-ratio transform
	log.cv <- log( cv )
	log.geom.mean <- mean( log.cv )
	return( log.cv - log.geom.mean )
}

######################################################################
# Take count data from one subject and use CLR transforms to send
# it to a Euclidean space with dimension (n.states * n.states)
######################################################################

count.mat.to.CLR <- function( transition.mat ) {
	# Do some paranoid error checking
	t.dim <- dim( transition.mat )
	n.states <- t.dim[1]
	if( t.dim[1] != t.dim[2] ) {
		warning( "Transition matrix is not square?!?" )
		return( NULL )
	}

	# Do all the CLR transforms
	cm.clr <- matrix( rep(0.0, n.states*n.states), nrow=n.states)
	for( j in 1:n.states ) {
		cm.clr[j,] <- count.vec.to.CLR( transition.mat[j,] )
	}

	# Munge everything into one giant vector
	return( as.vector(cm.clr) )
}

######################################################################
# Given a set of state trajectories for subjects who are all said
# to belong to the same cluster, get starting guesses for parameters
# for one component of our Markov-chain mixture
######################################################################

estimate.component.params <- function( subj.traj.list, n.states ) {
	# Estimate a posterior distribution for the initial state,
	# starting from a uniform prior.
	initial.state.alphas <- rep( 1.0, n.states )
	for( j in 1:length(subj.traj.list) ) {
		crnt.traj <- subj.traj.list[[j]]
		crnt.init.state <- crnt.traj[1]
		initial.state.alphas[crnt.init.state] = initial.state.alphas[crnt.init.state] + 1
	}

	# Now estimate parameters for distributions over the
	# rows of the transition matrices.
	n.subjects <- length( subj.traj.list )
	all.trans.probs <- array(
		rep(0.0, n.subjects*n.states*n.states),
		dim=c( n.subjects, n.states, n.states)
	)

	for( j in 1:n.subjects ) {
		# Build a matrix of counts
		crnt.count.mat <- state.traj.to.count.mat( subj.traj.list[[j]], n.states )

		# Convert rows of counts into estimated probabilities and store
		# them in our giant array of probabilities
		for( k in 1:n.states ) {
			all.trans.probs[j,k,] <- crnt.count.mat[k,] / sum(crnt.count.mat[k,])
		}
	}

	# Estimate Dirichlet distributions for each row of the transition matrix
	trans.mat.alphas <- matrix( rep( 0.0, n.states*n.states ), nrow=n.states )
	for( j in 1:n.states ) {
		# Get an n.subjects-by-n.states matrix that tabulates the
		# j-th rows of all the subjects' estimated transition matrices
		all.subj.row.mat <- all.trans.probs[,j,]

		# Estimate a Dirichlet distribution from which these
		# rows might have been drawn.
		fitted.params <- fit_dirichlet( all.subj.row.mat )
		trans.mat.alphas[j,] <- fitted.params$alpha
	}

	# Bundle up the results and return
	result <- list(
		"initial.state.Dirichlet.params" = initial.state.alphas,
		"trans.mat.Dirichlet.params" = trans.mat.alphas
	)

	return( result )
}

######################################################################
# Given a list of trajectories and a number of components,
# use k-means on CLR-transformed data to get initial cluster
# assignments
######################################################################

subject.data.to.mcmc.start <- function( k.clusters, subj.traj.list, n.states ) {
	# Examine the inputs
	n.subjects <- length( subj.traj.list )

	# Initialise the result as a list with no values
	count.mats <- vector(mode = "list", length = n.subjects)
	
	# Build all the count matrices
	count.mats <- lapply( subj.traj.list, 
		function(traj) { state.traj.to.count.mat( traj, n.states ) }
	)
	
	# Use CLR/k.means to get class assignments
	km.result <- kmeans.clusters( k.clusters, count.mats ) 
	
	# Estimate the cluster membership probs
	membership.probs <- km.result$size / n.subjects

	# For each cluster, estimate the parameters of (n.states + 1)
	# Dirichlet distributions: one for each row of a transition matrix
	# and one for the initial state.
	component.params <- list( NULL, k.clusters ) # will receive params
	for( j in 1:k.clusters ) {
		# Get the data for subjects in the j-th cluster
		crnt.subj.traj.list <- subj.traj.list[km.result$cluster == j]

		# Estimate params for the j-th cluster and record them
		crnt.params <- estimate.component.params( crnt.subj.traj.list, n.states )
		component.params[[j]] <- crnt.params
	}

	# Package up the results and return
	result <- list(
		"cluster" = km.result$cluster,
		"membership.probs" = membership.probs,
		"component.params" = component.params
	)

	return( result )
}

kmeans.CLR.clusters <- function( k.clusters, count.mats ) {
	# Examine the inputs
	n.subjects <- length( count.mats )
	n.states <- nrow( count.mats[[1]] )

	# Transform the rows of the count matrices so as to get
	# a matrix whose rows are points representing transition matrices.
	subj.points <- t( sapply( count.mats, count.mat.to.CLR ) )

	# Do k-means
	km.result <- kmeans( subj.points, k.clusters )
	return( km.result )
}

