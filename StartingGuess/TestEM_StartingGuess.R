####################################################################
#	Exercise the machinery for estimating a starting guess
#
# mrm: Whalley Range, 27 July 2022
####################################################################

rm( list=ls() ) # Wipe the slate clean.
source( "EM_StartingGuess.R" )
source( "MarkovChainMixtureUtils.R" )
source( "CLR_StartingGuess.R" )

###################################################################
#	Begin by testing some of the routines in EM_StartingGuess.R
##################################################################

# Check the lpmf
y <- c( 1, 0, 0 )
alpha <- c(1, 1, 1)
log.pmf <- DirMult.lpmf( y, alpha, normalised=TRUE )
pmf <- exp( log.pmf )
abs.diff.ratio <- abs( pmf - 1.0/3.0) / .Machine$double.eps
if( abs.diff.ratio < 1 ) {
	print( "DirMult.lpmf: passed test." ) 
} else { 
  print( paste( "DirMult.lpmf: failed test. abs.dif.ratio =", abs.diff.ratio ) )
}

# Check the gradient of the lpmf
deriv.h <- 1.0e-4
log.pmf <- DirMult.lpmf( y, alpha, normalised=TRUE )
computed.grad <- grad.DirMult.lpmf( y, alpha )
approx.grad <- rep( 0.0, length(alpha) )
for( j in 1:length(alpha) ) {
	delta.vec <- rep( 0.0, length(alpha) )
	delta.vec[j] <- deriv.h
	log.pmf.plus <- DirMult.lpmf( y, alpha + delta.vec, normalised=TRUE )
	log.pmf.minus <- DirMult.lpmf( y, alpha - delta.vec, normalised=TRUE )
	approx.grad[j] <- (log.pmf.plus - log.pmf.minus) / (2.0 * deriv.h)
}

max.diff <- max(abs(approx.grad - computed.grad))
if( max.diff < deriv.h*deriv.h ) {
	print( "grad.DirMult.lpmf: passed test." ) 
} else { 
  print( "grad.DirMult.lpmf: failed test." )
}

# Check the probability for a count matrix
y.mat <- diag( nrow=3 ) # A 3-by-3 identity matrix
alpha.mat <- matrix( rep( 1.0, 9 ), nrow=3 )
log.pmf <- count.mat.lpmf( y.mat, alpha.mat, normalised=TRUE )
pmf <- exp( log.pmf )
abs.diff.ratio <- abs( pmf - 1.0/27.0) / .Machine$double.eps
if( abs.diff.ratio < 1 ) {
	print( "count.mat.lpmf: passed test." ) 
} else { 
  print( paste( "count.mat.lpmf: failed test. abs.dif.ratio =", abs.diff.ratio ) )
}

# Check the parameter-packing utilities
y1 <- matrix( 1:9, nrow=3 )
y2 <- 2 * y1
orig.mat.list <- list( y1, y2 )
vec <- alpha.matrices.to.log.vec( orig.mat.list )
new.mat.list <- alpha.log.vec.to.matrices( vec, 3, 2 )
verdicts <- sapply( 1:2, 
	function(j){ 
		orig.vec <- as.vector( orig.mat.list[[j]] )
		new.vec <- as.vector( new.mat.list[[j]] )
		avg.abs <- 0.5 * (abs(orig.vec)  +  abs(new.vec))
		abs.diff <- abs(orig.vec - new.vec)
		ratios <- abs.diff / avg.abs
		return( (max(ratios) / .Machine$double.eps) < 10.0 ) 
		
	} 
) 
if( all( verdicts ) ) {
	print( "param packing: passed test." ) 
} else { 
  print( "param packing: failed test." )
}

# Check part of the convergence test
v1 <- c (1, 1.9, 3 )
v2 <- c( -1, 2.1, 3.1 )
result <- max.frac.abs.diff( v1, v2 )
abs.diff.ratio <- abs( result - 0.1) / .Machine$double.eps
if( abs.diff.ratio < 1 ) {
	print( "max.frac.abs.diff: passed test." ) 
} else { 
  print( "max.frac.abs.diff: failed test." )
}

###################################################################
#	Read the test data, then use the true params to test various
# components of the machinery.
##################################################################

load( "../TestData.RData" )

# Get the correct dimensions of things
n.states <- true.DirichletMC.mix$n.states
n.components <- true.DirichletMC.mix$n.components

# Reduce the state trajectories to count matrices
count.mats <- trajectories.to.count.mats( traj.list, n.states, 0.0 )

# Repackage the true params in the form expected by EM_StartingGuess
true.params <- list(
	n.states = true.DirichletMC.mix$n.states,
	n.components = true.DirichletMC.mix$n.components,
	pi = true.DirichletMC.mix$pi,
	alpha.matrices = lapply( true.DirichletMC.mix$dirichlet.params,
			function( params ) { params$transition.matrix.alphas }
	)
)

# Evaluate class membership probs based on the true parameters
log.gamma.mat <- log.membership.probs( count.mats, true.params )
gamma.mat <- exp( log.gamma.mat )
fitted.comp <-  apply( gamma.mat, 1, which.max )
true.comp <- synthetic.data$true.component
print( cbind( gamma.mat, true.comp, fitted.comp ) )
all( fitted.comp == true.comp )

# Initialise the EM algorithm with the true class assignments
class.assignments.to.EM.start( true.comp, count.mats )

# Run the EM algorithm
em.guess <- fit.Dirichlet.Markov.mixture( n.components, traj.list, n.states, true.comp, alpha.penalty=0.005 ) 
em.guess

###################################################################
#	Finally, test the code by trying to recover the original
# params and class assignments without using any of the true
# data beyond n.states and n.components
##################################################################

em.guess <- fit.Dirichlet.Markov.mixture( n.components, traj.list, n.states, alpha.penalty=0.005 ) 
em.guess


