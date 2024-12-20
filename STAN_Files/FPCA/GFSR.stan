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
	int<lower=0> N;                  // number of time series
	int<lower=0> D;                  // grid length
	int<lower=0> p;                  // number of fixed effects

	int Kt;                          // number of spline basis functions
	int Kp;                          // number of PC basis functions
	
	matrix[N, D] Y;                  // outcome matrix
	vector[p] X[N];                  // fixed effect design matrix
	matrix[D, Kt] BS;                // B-spline evaluation matrix
	
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
	matrix[Kp,Kt] beta_psi;          // matrix of PC spline coefficients
	vector[Kp] c[N];                 // matrix of PC loadings
	vector<lower=0,upper=100>[p] beta_sig;     // tuning variance
	vector<lower=0,upper=100>[Kp] psi_sig;     // tuning variance
	real<lower=0> sigma2;            // Error in observation
}

transformed parameters {
  vector[D] UDV[N];
	
	vector<lower=0.01>[p] beta_tau2;    // tuning parameter
	vector<lower=0.01>[Kp] psi_tau2;    // tuning parameter
  
  for(i in 1:N){
    UDV[i] = (BS * beta_psi') * c[i];
  }
  
	for(pcur in 1:p) {
		beta_tau2[pcur] = pow(beta_sig[pcur], -1);
	}

	for(kcur in 1:Kp) {
		psi_tau2[kcur] = pow(psi_sig[kcur], -1);
	}	

}


model {
	/////////////////////////////////////
	// Prior distributions
	/////////////////////////////////////

	// Prior for variance components controlling smoothness in beta
	for (pcur in 1:p) {
		beta_sig[pcur] ~ inv_gamma(.001,.001); 
	}
	
	// Prior for variance components controlling smoothness
	for (kcur in 1:Kp) {
		psi_sig[kcur] ~ inv_gamma(.001,.001); 
	}
			
	// Prior for spline coefficients for beta
	for (pcur in 1:p) {
		(beta[pcur])' ~ multi_normal_prec(mu_beta, beta_tau2[pcur] * PenMat);
	}

	// Prior for spline coefficients for psi1
	for (kcur in 1:Kp) {
		(beta_psi[kcur])' ~ multi_normal_prec(mu_beta, psi_tau2[kcur] * PenMat);
	}
	
	// Prior for PC scores
	for (i in 1:N) {
		c[i] ~ normal(0.0, 1.0);
	}
	
	// Prior on error
  sigma2 ~ inv_gamma(1, 0.01);

	/////////////////////////////////////
	// Outcome likelihood
	/////////////////////////////////////
	for (i in 1:N) {
	  Y[i] ~ normal((BS * beta') * X[i] +       // fixed effects
			            UDV[i], sqrt(sigma2));      // PC effects
	}
}
