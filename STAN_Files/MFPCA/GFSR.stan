////////////////////////////////////////////////////////////
// February 1 2014
// Jeff Goldsmith
//
// This file contains code for multilevel bayesian 
// generalized function-on-scalar regression with MFPCA for 
// functional residuals. In this code we use a binary 
// outcome and logit link.
////////////////////////////////////////////////////////////

data {
	int<lower=0> I;                  // number of subjects
	int<lower=0> J;                  // number of visits per subject
	int<lower=0> IJ;                 // total number of subjects
	int<lower=0> D;                  // grid length
	int<lower=0> p;                  // number of fixed effects

	int Kt;                          // number of spline basis functions
	int Kp1;                         // number of PC basis functions
	int Kp2;                         // number of PC basis functions
	
	matrix[IJ, D] Y;                 // outcome matrix
	vector[p] X[IJ];                 // fixed effect design matrix
	matrix[D, Kt] BS;                 // B-spline evaluation matrix
	
  cov_matrix[Kt] PenMat;           // prior precision matrix for spline effects (penalty matrix)
}

transformed data {
  vector[Kt] mu_beta;              // prior mean for spline effects

	for (k in 1:Kt) {
		mu_beta[k] = 0;
	}
}

parameters {
	matrix[p,Kt] beta;               // matrix of fixed effect spline coefficients
	matrix[Kp1,Kt] beta_psi1;        // matrix of PC spline coefficients
	matrix[Kp2,Kt] beta_psi2;        // matrix of PC spline coefficients
	vector[Kp1] c1[I];               // matrix of PC loadings
	vector[Kp2] c2[IJ];              // matrix of PC loadings
	vector<lower=0,upper=100>[p] beta_sig;     // tuning variance
	vector<lower=0,upper=100>[Kp1] psi1_sig;   // tuning variance
	vector<lower=0,upper=100>[Kp2] psi2_sig;   // tuning variance
	real<lower=0> sigma2;            // Error in observation
}

transformed parameters {
	vector<lower=0.01>[p] beta_tau2;    // tuning parameter
	vector<lower=0.01>[Kp1] psi1_tau2;  // tuning parameter
	vector<lower=0.01>[Kp2] psi2_tau2;  // tuning parameter
  
	for(pcur in 1:p) {
		beta_tau2[pcur] = pow(beta_sig[pcur], -1);
	}

	for(kcur in 1:Kp1) {
		psi1_tau2[kcur] = pow(psi1_sig[kcur], -1);
	}	

	for(kcur in 1:Kp2) {
		psi2_tau2[kcur] = pow(psi2_sig[kcur], -1);
	}	
}


model {
	/////////////////////////////////////
	// Prior distributions
	/////////////////////////////////////

	// Prior for variance components controlling smoothness in beta
	for (pcur in 1:p) {
		beta_sig[pcur] ~ inv_gamma(.001,.001); 
//		beta_sig[pcur] ~ uniform(0,25); 
//		beta_sig[pcur] ~ cauchy(0, 5); 
	}
	
	// Prior for variance components controlling smoothness in psi1
	for (kcur in 1:Kp1) {
		psi1_sig[kcur] ~ inv_gamma(.001,.001); 
//		psi1_sig[kcur] ~ uniform(0,25); 
//		psi1_sig[kcur] ~ cauchy(0, 5);
	}
	
	// Prior for variance components controlling smoothness in psi2
	for (kcur in 1:Kp2) {
		psi2_sig[kcur] ~ inv_gamma(.001,.001); 
//		psi2_sig[kcur] ~ uniform(0,25); 
//		psi2_sig[kcur] ~ cauchy(0, 5);
	}
			
	// Prior for spline coefficients for beta
	for (pcur in 1:p) {
		(beta[pcur])' ~ multi_normal_prec(mu_beta, beta_tau2[pcur] * PenMat);
	}

	// Prior for spline coefficients for psi1
	for (kcur in 1:Kp1) {
		(beta_psi1[kcur])' ~ multi_normal_prec(mu_beta, psi1_tau2[kcur] * PenMat);
	}

	// Prior for spline coefficients for psi2
	for (kcur in 1:Kp2) {
		(beta_psi2[kcur])' ~ multi_normal_prec(mu_beta, psi2_tau2[kcur] * PenMat);
	}
	
	// Prior for subject-level PC scores
	for (i in 1:I) {
		c1[i] ~ normal(0.0, 1.0);
	}

	// Prior for subject/visit-level PC scores
	for (i in 1:IJ) {
		c2[i] ~ normal(0.0, 1.0);
	}
	
	// Prior on error
  sigma2 ~ inv_gamma(1, 0.01);

	/////////////////////////////////////
	// Outcome likelihood
	/////////////////////////////////////
	for (i in 1:I) {
		for (j in 1:J) {
			Y[((i-1)*J+j)] ~ normal((BS * beta') * X[((i-1)*J+j)] +      // fixed effects
			                        (BS * beta_psi1') * c1[i] +          // subject-level PC effects
			                        (BS * beta_psi2') * c2[((i-1)*J+j)], // subject/visit-level PC effects
			                        sqrt(sigma2));
		}
	}
}
