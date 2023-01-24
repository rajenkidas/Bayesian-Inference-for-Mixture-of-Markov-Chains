######################################################################
#	Here are some utilities to support the fitting of mixtures
# of Markov chains to lists of state trajectories. 
#
# mrm:	Triumphstadt, Aalen, 11 Aug 2022
######################################################################

######################################################################
#	A utility to get the data into a useful form
######################################################################

trajectories.to.count.mats <- function( traj.list, n.states, pseudocount ) {
	# Examine the inputs
	n.subjects <- length( traj.list )

	# Initialise the result as a list with no values
	count.mats <- vector(mode = "list", length = n.subjects)
	
	# Build all the count matrices
	count.mats <- lapply( traj.list, 
		function(traj) { state.traj.to.count.mat( traj, n.states, pseudocount ) }
	)
	
	return( count.mats ) 
}

######################################################################
# Take a state trajectory and return a matrix whose i, j entry
# is the number of observed transitions from state to to state j.
######################################################################

state.traj.to.count.mat <- function( traj, n.states, pseudocount=0.5 ) {
	count.mat <- matrix( rep(pseudocount, n.states*n.states), nrow=n.states )
	for( j in 2:length(traj) ) {
		from <- traj[j-1]
		to <- traj[j]
		count.mat[from, to] <- count.mat[from, to] + 1
	}

	return( count.mat )
}


