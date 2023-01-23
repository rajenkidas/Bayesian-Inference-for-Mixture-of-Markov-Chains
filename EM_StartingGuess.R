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
#	(2) Fit a mixture of Markov chains to the counts, using
#		a Dirichlet-Multinomial model for the rows of the 
#		transition matrix. This is a generalisation of the
#		methods in:
#
#			I. Holmes, K. Harris and C. Quince (2012), 
#			Dirichlet multinomial mixtures: Generative models 
#			for microbial metagenomics, PLoS ONE, 7:1–15.
#			DOI: 10.1371/journal.pone.0030126
#
# mrm:	Whalley Range, 26 & 27 July 2022
######################################################################

library( MCMCprecision ) # fit_dirichlet()

# Set switches that determine whether we do some  (potentially
# computationally expensive) tests. These should be FALSE in
# production code.
test.grad.loglike <- FALSE

######################################################################
# 	Take a state trajectory and return a matrix whose i, j entry
# is the number of observed transitions from state to to state j.
######################################################################

state.traj.to.count.mat <- function( traj, n.states ) {
	count.mat <- matrix( rep(0.5, n.states*n.states), nrow=n.states )
	for( j in 2:length(traj) ) {
		from <- traj[j-1]
		to <- traj[j]
		count.mat[from, to] <- count.mat[from, to] + 1
	}

	return( count.mat )
}

######################################################################
# 	Given a vector of logs {v_1, v_2, ..., v_n}, compute
# log( \sum_{j=1}^n exp( v_j ) )
# See http://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
######################################################################

log.sum.exp <- function( v ) {
	v.max <- max( v )
	u <- v - v.max 
	sum.exp.u <- sum( exp(u) )
	return( log(sum.exp.u) + v.max )
}

######################################################################
#	Given a vector of shape parameters (alpha_1, ..., alpha_n),
# compute B(alpha) such that the normalisation factor for the 
# Dirichlet distribution is 1 / B(alpha).
######################################################################

log.Dirichlet.B <- function( alpha ) {
	# Use cool, fast, vectorized calculations
	logB <- sum( lgamma(alpha) ) ;
	logB <- logB - lgamma( sum(alpha) ) ;
	
	return( logB ) ;
}

######################################################################
#	Given a vector of shape parameters (alpha_1, ..., alpha_n),
# compute the gradient of B(alpha). The key fact is that 
# digamma() is the derivative of the log of the gamma function
######################################################################

grad.log.Dirichlet.B <- function( alpha ) {
	sum.of.alphas <- sum(alpha)
	grad.log.B <- digamma( alpha ) 
	grad.log.B <- grad.log.B - digamma( sum.of.alphas ) 
	return( grad.log.B )
}

######################################################################
#	Given a vector of counts (y_1, ..., y_n) and a vector of
# shape parameters (alpha_1, ..., alpha_n), compute the probability
# mass function for observing the given pattern of counts.
######################################################################

DirMult.lpmf <- function( y, alpha, normalised=FALSE ) {
	# Compute those terms that depend on alpha
	result <- log.Dirichlet.B( alpha + y ) ;
	result <- result - log.Dirichlet.B( alpha ) ;

	# If desired, include those that depend on the
	# counts alone. These don't affect inference for the alphas.
	if( normalised ) {
		result <- result - sum( lgamma( y + 1 ) ) ;
		result <- result + lgamma( sum(y) + 1 ) ;
	}
	
	return( result ) ;
}

######################################################################
#	Given a vector of counts (y_1, ..., y_n) and a vector of
# shape parameters (alpha_1, ..., alpha_n), compute the gradient
# (with respect to the shape params) of the probability
# mass function for the given pattern of counts.
######################################################################

grad.DirMult.lpmf <- function( y, alpha ) {
	grad.lpmf <- grad.log.Dirichlet.B( alpha + y ) 
	grad.lpmf <- grad.lpmf - grad.log.Dirichlet.B( alpha ) 
	
	return( grad.lpmf )
}

######################################################################
#	Given a matrix of transition counts and a matrix of shape
# parameters characterising a distribution over transition matrices,
# compute the probability of getting the counts.
######################################################################

count.mat.lpmf <- function( counts, alphas, normalised=FALSE ) {
	n.states <- nrow( counts )
	row.lpmf.vals <- sapply( 1:n.states, 
		function(j) {
			DirMult.lpmf( counts[j,], alphas[j,], normalised )
		}
	)
		
	return( sum(row.lpmf.vals) ) ;
}

######################################################################
#	Given a set of count matrices and the parameters of a 
# Dirichlet-Multinomial mixture of Markov chains, compute a matrix
# L whose entries are
#
#	L_{j,k} = log( P( count.mat[j] | alpha.mat[k] ) )
#
# If normalised is FALSE, then we leave out certain combinatorial
# factors that depend only on the counts.
######################################################################

log.likelihood.matrix <- function( count.mats, mixture.params, normalised=FALSE ) {
	# Get various useful sizes
	n.states <- nrow( count.mats[[1]] )
	n.subjects <- length( count.mats )
	n.components <- length( mixture.params$pi )
	
	# Initialise the result
	mat.entries <- rep( 0.0, n.subjects*n.components )
	log.like.mat <- matrix( mat.entries, ncol=n.components )
	
	# Now fill it up
	log.pi <- log( mixture.params$pi )
	for( j in 1:n.subjects ) {
		log.like.vals <- sapply( mixture.params$alpha.matrices, 
			function( alpha.mat ) {
				count.mat.lpmf( count.mats[[j]], alpha.mat, normalised )
			}
		)
		
		log.like.mat[j,] <- log.like.vals + log.pi
	}
		
	return( log.like.mat )
}

######################################################################
#	Given a set of count matrices and the parameters of a 
# Dirichlet-Multinomial mixture of Markov chains, compute a matrix
# gamma of membership probabilities where 
#
#	gamma_{i,j} = P( countMat_i comes from chain_j )
######################################################################

log.membership.probs <- function( count.mats, mixture.params ) {
	# Get various useful sizes
	n.components <- length( mixture.params$pi )
	
	# Get the terms that depend on the alphas and counts
	log.raw.gamma.mat <- log.likelihood.matrix( count.mats, mixture.params, FALSE )
	
	# Normalise the rows
	log.normalisation <- apply( log.raw.gamma.mat, 1, log.sum.exp ) 
	mat.entries <- rep( log.normalisation, times=n.components )
	normalisation.mat <- matrix( mat.entries, ncol=n.components )
	log.gamma.mat <- log.raw.gamma.mat - normalisation.mat
	
	return( log.gamma.mat )
}

######################################################################
#	Given a set of count matrices and the parameters of a 
# Dirichlet-Multinomial mixture of Markov chains compute the
# log-likelihood.
######################################################################

log.likelihood <- function( count.mats, mixture.params, normalised=FALSE ) {
	log.like.mat <- log.likelihood.matrix( count.mats, mixture.params, normalised )
	loglike.per.subject <- apply( log.like.mat, 1, log.sum.exp ) 
	return( sum(loglike.per.subject) )
}

######################################################################
#	Given a set of count matrices and the parameters of a 
# Dirichlet-Multinomial mixture of Markov chains compute the
# gradient of the log-likelihood. It gets stored in a list of 
# arrays that are the same shape as alpha.mats
######################################################################

grad.log.likelihood <- function( count.mats, mixture.params ) {
	# Get various useful sizes
	n.states <- nrow( count.mats[[1]] )
	n.subjects <- length( count.mats )
	n.components <- length( mixture.params$pi )
	
	# Get the membership probs
	log.gamma.mat <- log.membership.probs( count.mats, mixture.params )
	gamma.mat <- exp( log.gamma.mat )
	
	# Get the gradients 
	n.states.sq <- n.states * n.states
	zero.mat <- matrix( rep(0.0, n.states.sq), nrow=n.states )
	grad.mats <- lapply( 1:n.components, function(j) { zero.mat } )
	for( i in 1:n.subjects ) {
		count.mat <- count.mats[[i]]
		for( j in 1:n.components ) {
			alpha.mat <- mixture.params$alpha.matrices[[j]]
			for( k in 1:n.states ) {
				crnt.grad <- grad.DirMult.lpmf( count.mat[k,], alpha.mat[k,] )
				weighted.grad <- gamma.mat[i,j] * crnt.grad 
				grad.mats[[j]][k,] <- grad.mats[[j]][k,] + weighted.grad
			}
		}
	}
	
	return( grad.mats )
}

######################################################################
#	A pair of utilities to convert back and forth between a 
# list of alpha matrices and a vector of the logs of their entries.
# We use logs so as to avoid the constraint that alpha > 0.
######################################################################

alpha.matrices.to.log.vec <- function( alpha.matrices ) {
	vec.list <- lapply( alpha.matrices, as.vector )
	vec <- do.call( c, vec.list )
	return( log(vec) )
}

alpha.log.vec.to.matrices <- function( log.vec, n.rows, n.matrices ) {
	n.entries <- n.rows * n.rows # number of entries in one matrix
	stopifnot( length(log.vec) == n.matrices * n.entries ) # paranoid error checking
	mat.list <- lapply( 1:n.matrices,
		function(j) {
			first <- 1 + (j-1)*n.entries
			last <- j*n.entries
			matrix( exp( log.vec[first:last]) , nrow=n.rows )
		}
	)
	
	return( mat.list )
}

######################################################################
#	Given count data and membership probs, estimate matrices of
# shape parameters
######################################################################

estimate.alphas <- function( count.mats, log.gamma.mat ) {
	# Assign each subject to their most probable component
	comp.num <- apply( log.gamma.mat, 1, which.max )
	
	# Get sizes of things
	n.components <- ncol( log.gamma.mat )
	n.states <- nrow( count.mats[[1]] )
	
	# Initialise the result with all 1's corresponding to a flat 
	# Dirichlet distrib.
	flat.mat <- matrix( rep(1, n.states*n.states), nrow=n.states )
	alpha.mats <- lapply( 1:n.components, function(j) { flat.mat } )
	
	# Within each group, estimate shape parameters
	for( i in 1:n.components ) {
		# Get the count matrices for  
		# subjects assigned to component i.
		my.counts <- count.mats[comp.num == i]
		
		# Fit a dirichlet to each row in turn
		crnt.alpha.mat <- flat.mat
		if( length(my.counts) > 0 ) {
			for( j in 1:n.states ) {
				# Build an n.subjects-by-n.states matrix whose
				# i-th row is the j-th row of subject i's count matrix
				tmp <- t(sapply( my.counts, function(mat) {  mat[j,] } ))
			
				# Drop any rows that are all zeroes
				row.sums <- rowSums( tmp )
				idx <- (row.sums != 0)
				row.mat <- tmp[idx,]
			
				# Reduce the counts to proportions
				prob.mat <- t(apply( row.mat, 1, function(v) { v / sum(v) } ))
			
				# Fit a dirichlet to the proportions
				fitted.params <- fit_dirichlet( prob.mat )
				crnt.alpha.mat[j,] <- fitted.params$alpha
			}
		}
		
		alpha.mats[[i]] <- crnt.alpha.mat
	}
	
	return( alpha.mats )
}

######################################################################
#	Given membership probs, make a maximum-likelihood estimate the 
# parameters of the mixture. We can't do this analytically, so 
# need to use an optimiser.
######################################################################

maxlike.mixture.params <- function( count.mats, 
                                    log.gamma.mat, 
                                    alpha.penalty = 0.0, 
                                    alpha.mats = NULL ) {
	# Get various useful sizes
	n.components <- ncol( log.gamma.mat )
	n.subjects <- nrow( log.gamma.mat )
	n.states <- nrow( count.mats[[1]] )
	n.states.sq <- n.states * n.states # number of entries in a transition matrix
	
	# Set up starting guesses for the parameters
	gamma.mat <- exp( log.gamma.mat )
	if( is.null(alpha.mats) ) {
		alpha.mats <- estimate.alphas( count.mats, gamma.mat )
	} 
	
	# Assemble the starting guesses into a list
	crnt.params <- list(
		pi = colSums( gamma.mat ) / n.subjects,
		alpha.matrices = alpha.mats
	)
	
	# Specify the optimisation target
	neg.loglike <- function( log.alpha.vec ) {
		# Repackage the alphas as a list of matrices
		alpha.mats <- alpha.log.vec.to.matrices( log.alpha.vec, n.states, n.components )
		
		# Make a copy of the parameter object
		params <- crnt.params
		params$alpha.matrices <- alpha.mats
		
		# Compute the log-likelihood
		loglike <- log.likelihood( count.mats, params )
		
		# Add a penalty for the sum of all the alphas
		loglike <- loglike - alpha.penalty * sum( exp(log.alpha.vec) )
		return( -loglike )
	}
	
	# Also its gradient
	grad.neg.loglike <- function( log.alpha.vec ) {
		# Repackage the alphas as a list of matrices
		alpha.mats <- alpha.log.vec.to.matrices( log.alpha.vec, n.states, n.components )
		
		# Make a copy of the parameter object
		params <- crnt.params
		params$alpha.matrices <- alpha.mats
		
		# Build a list of matrices whose entries are components
		# of the gradient w.r.t. the alphas, then repackage it as a vector
		grad.mats <- grad.log.likelihood( count.mats, params )
		grad.list <- lapply( grad.mats, as.vector )
		grad <- do.call( c, grad.list )
		
		# Include the penalty term in the gradient
		grad <- grad - alpha.penalty
		
		# Take care of the fact that we're optimising with log(alpha) 
		grad <- grad * exp(log.alpha.vec)  # elementwise product
		return( -grad )
	}
	
	# If we're doing testing, check that gradient works as it should.
	# These tests are computationally expensive, so should be turned
	# off in production code.
	if( test.grad.loglike ) { # See note near the top of the file
		# Get the parameters we'll actually send to optim(): they're 
		# logs of shape parameters.
		log.alpha.vec <- alpha.matrices.to.log.vec( alpha.mats )
		n.alphas <- length( log.alpha.vec )
		
		# Estimate the gradient with finite differences.
		deriv.eps <- 1.0e-4
		fd.grad <- rep( 0.0, n.alphas )
		for( j in 1:n.alphas) {
			# Use the approximation (f(x+h) - f(x - h))/2h
			log.alpha.vec.plus <- log.alpha.vec
			log.alpha.vec.plus[j] <- log.alpha.vec.plus[j] + deriv.eps
			neg.loglike.plus <- neg.loglike( log.alpha.vec.plus )
			
			log.alpha.vec.minus <- log.alpha.vec
			log.alpha.vec.minus[j] <- log.alpha.vec.minus[j] - deriv.eps
			neg.loglike.minus <- neg.loglike( log.alpha.vec.minus )
			
			fd.grad[j] <- (neg.loglike.plus - neg.loglike.minus) / (2.0 * deriv.eps)
		}
		
		# Evaluate the gradient using our code
		direct.grad <- grad.neg.loglike( log.alpha.vec  )
		
		# Compute the ratio of difference over mean absolute value.
		mean.abs.grad <- 0.5 * sum( abs(fd.grad) + abs(direct.grad) )
		grad.diff <- fd.grad - direct.grad
		ratio <- grad.diff / mean.abs.grad
		print( direct.grad )
		print( fd.grad )
		print( ratio )
		
		# Force the run to stop after this step
		max.cycles <- 1 
	}

	# Invoke optim
	alpha.guess <- alpha.matrices.to.log.vec( crnt.params$alpha.matrices )
	result <- optim( alpha.guess, fn=neg.loglike, gr=grad.neg.loglike, method="BFGS", hessian=TRUE )
	if( result$convergence != 0 ) {
		warning( "optim() failed to converge" )
	}
	
	crnt.params$alpha.matrices <- alpha.log.vec.to.matrices( result$par, n.states, n.components )
	return( crnt.params ) 
}

######################################################################
#	A utility: given a list of class assignments, generate a 
# starting guess for the EM algorithm.
######################################################################

class.assignments.to.EM.start <- function( comp.nums, count.mats ) {
	# Work out the sizes of things
	n.subjects <- length( comp.nums )
	n.components <- max( comp.nums )
	
	# Set up a class membership matrix using the class assignments,
	# but make sure all entries are non-zero.
	gamma.mat.eps = 0.00001 / n.states
	small.entry <- log( gamma.mat.eps / (n.components - 1) )
	large.entry <- log( 1.0 - gamma.mat.eps )
	mat.entries <- rep( small.entry, n.subjects*n.components )
	log.gamma.mat <- matrix( mat.entries, nrow=n.subjects )
	for( j in 1:n.subjects ) {
		k <- comp.nums[j]
		log.gamma.mat[j,k] <- large.entry
	}
	
	# Use the class assignments to get starting guesses for the parameters
	gamma.mat <- exp( log.gamma.mat )
	initial.params <- list(
		pi = colSums( gamma.mat ) / n.subjects,
		alpha.matrices = estimate.alphas( count.mats, log.gamma.mat )
	)
	
	return( initial.params )
}

######################################################################
#	The main event: given a collection of count matrices, use
# the EM algorithm to fit a mixture model to them.
######################################################################

fit.Dirichlet.Markov.mixture <- function( 
	n.components, traj.list, n.states, 
	initial.comp.nums=NULL, alpha.penalty=0.0,
	max.cycles=100, tol=1.0e-4 ) {
	# Get sizes
	n.subjects <- length( traj.list )
	
	# Reduce the state trajectories to count matrices
	count.mats <- trajectories.to.count.mats( traj.list, n.states, 0.0 ) # final arg is a pseudocount
	
	# If need be, make an initial guess for the class assignments
	if( is.null(initial.comp.nums) ) {
		# Use our old kmeans machinery to get a starting guess.
		pseudo.count <- 0.01 # ensures entries in the count matrices are all non-zero (needed for CLR)
		km.count.mats <- trajectories.to.count.mats( traj.list, n.states, pseudo.count ) 
		km.result <- kmeans.CLR.clusters( n.components, km.count.mats )
		
		# Set comp.nums from the kmeans result
		initial.comp.nums <- km.result$cluster
	}
		
	# Use the class assignments to estimate the params
	initial.params <- class.assignments.to.EM.start( initial.comp.nums, count.mats )
		
	# Plunge into the EM algorithm
	n.cycles <- 0
	converged <- FALSE
	crnt.params <- initial.params
	while( (n.cycles < max.cycles) && !converged ) {
		# Recompute class membership probs
		log.gamma.mat <- log.membership.probs( count.mats, crnt.params )
		
		# Re-estimate the parameters
		next.params <- maxlike.mixture.params( count.mats, log.gamma.mat, alpha.penalty )
		
		# Check for convergence 
		converged <- convergence.test( crnt.params, next.params, tol )
		
		# Get ready for the next cycle
		n.cycles <- n.cycles + 1 
		crnt.params <- next.params
	}
	
	if( !converged ) {
		warning( "fit.Dirichlet.Markov.mixture() didn't converge." ) ;
	}
	
	crnt.params$converged <- converged
	return( crnt.params )
}

######################################################################
#	A convergence test
######################################################################

max.frac.abs.diff <- function( vec.a, vec.b ) {
	avg.abs <- 0.5 * (abs(vec.a) + abs(vec.b)) 
	abs.diff <- abs(vec.a - vec.b) 
	idx <- avg.abs  > 0.0
	return( max( abs.diff[idx] / avg.abs[idx] ) )
}

# Check for small fractional changes in pi and the alphas
convergence.test <- function( crnt.params, next.params, tol ) {
	# Get sizes
	n.components <- length( crnt.params$pi )
	n.states <- nrow( crnt.params$alpha.matrices[[1]] )
	
	# Do the test
	verdict <- (max.frac.abs.diff(crnt.params$pi, next.params$pi) < tol ) 
	for( i in 1:n.components ) {
		mat.vec.A <- as.vector( crnt.params$alpha.matrices[[i]] )
		mat.vec.B <- as.vector( next.params$alpha.matrices[[i]] )
		verdict <- verdict && (max.frac.abs.diff(mat.vec.A, mat.vec.B) < tol)
	}
	
	return( verdict )
}
