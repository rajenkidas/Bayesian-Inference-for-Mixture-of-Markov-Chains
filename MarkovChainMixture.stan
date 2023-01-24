functions {
	real logDirichletB( vector alpha ) {
	//**************************************************
	//	Given a vector alpha = (alpha_1, ..., alpha_n),
	// compute B(alpha) such that the normalisation
	// factor for the Dirichlet distribution
	// is 1 / B(alpha) 
	//
		real logB ;
	
		// Do cool, fast, vectorized calculations
		logB = sum( lgamma(alpha) ) ;
		logB -= lgamma( sum(alpha) ) ;
		
		return( logB ) ;
	}

	real DirMultRaw_lpmf( int[] y, vector alpha ) {
	//**************************************************
	//	Compute the alpha-dependent part of the log
	// of the Dirichlet-Multinomial distribution.
	//
		real result ;
		
		result = logDirichletB( alpha + to_vector(y) ) ;
		result -= logDirichletB( alpha ) ;
		
		return( result ) ;
	}
	
	real DirMultNormalised_lpmf( int[] y, vector alpha ) {
	//**************************************************
	//	Compute the log of the Dirichlet-Multinomial 
	// distribution, including normalising factors that
	// depend on the counts y.
	//
		real result = DirMultRaw_lpmf( y | alpha ) ;
		
		result -= sum( lgamma( to_vector(y) + 1.0 ) ) ;
		result += lgamma( sum(y) + 1.0 ) ;
		return( result ) ;
	}
}

data {
	// Fix the sizes of things
	int<lower = 1>		nSubjects ;
	int<lower = 2>		nStates ;
	int<lower = 1>		nComponents ;
	
	// The observed counts: for each subject, an
	// nStates-by-nStates matrix of counts.
	int<lower = 0>		counts[nSubjects, nStates, nStates] ;
	
	// Fix the parameters of the priors
	real<lower = 0>		thetaPriorAlpha ;
	real<lower = 0>		alphaPriorAlpha ;
	real<lower = 0>		alphaPriorBeta ;
}

transformed data {
	// Assemble parameters for priors
	vector[nComponents] weightPriorAlphaVec ;	// Dirichlet
	vector[nStates]	thetaPriorAlphaVec[nStates]  ;	// Dirichlet
	real alphaGrandTotalPriorAlpha = nStates * nStates * alphaPriorAlpha ; // Gamma
	
	// The Dirichlet distrib takes a vector of alphas as
	// its parameters. For the matrices of transition probs,
	// we want them all to be equal to thetaAlpha.
	for( j in 1:nStates ) { 
		for( k in 1:nStates ) { 
			thetaPriorAlphaVec[j, k] = thetaPriorAlpha ; 
		}
	} 
	
	// For the mixture weights, we want a uniform prior over
	// the relevant simplex.
	for( j in 1:nComponents ) {
		weightPriorAlphaVec[j] = 1.0 ;
	}
}

parameters {
	simplex[nComponents]	weightOmega ;
	ordered[nComponents]	alphaGrandTotal ;
	simplex[nStates]		phiSmplx[nComponents] ;
	simplex[nStates]		theta[nComponents, nStates] ;		
}

transformed parameters {
	// Here we construct the parameters of the 
	// Dirichlet distribs from the alphaGrandTotals,
	// the phis and the thetas
	vector[nStates]		alphaHat ;
	vector[nStates]		alpha[nComponents, nStates] ;
	
	for( j in 1:nComponents ) {
		alphaHat = alphaGrandTotal[j] * phiSmplx[j] ;
		for( k in 1:nStates ) {
			alpha[j,k] = alphaHat[k] * theta[j,k] ;
		}
	}
}

model {
	// Declare variables to hold intermediate results
	vector[nComponents]	logWeights ;
	vector[nComponents]	logProbsPerComponent ;
	
	// Add contributions from the priors
	weightOmega ~ dirichlet( weightPriorAlphaVec ) ;
	for( j in 1:nComponents ) {
		alphaGrandTotal[j] ~ gamma( alphaGrandTotalPriorAlpha, alphaPriorBeta ) ;
		for( k in 1:nStates ) {
			theta[j,k] ~ dirichlet( thetaPriorAlphaVec[k] ) ;
		}
	}
	
	// Add the contributions from the likelihood
	logWeights = log( weightOmega ) ;
	for( i in 1:nSubjects ) {
		logProbsPerComponent = logWeights ;
		for( j in 1:nComponents ) {
			for( k in 1:nStates ) {
		  		logProbsPerComponent[j] += DirMultRaw_lpmf( counts[i,k] | alpha[j,k] ) ;
		  	}
		}
	  
		target += log_sum_exp( logProbsPerComponent ) ;
	}
}



