data {
  int N;    // Number of time series
  int M;    // Number of observations
  int Q;    // Number of spline bases
  int K;    // Number of Eigenfunctions
  
  matrix[N, M] Y;   // Original data
  matrix[M, Q] B;   // Orthogonalized basis
  matrix[Q, Q] P;  // Penalty matrix for splines
}

parameters {
  real<lower=0> sigma2; // Error in observation
  
  // Fixed-effect components
  vector[Q] w_mu;               // Population mean parameters
  real<lower=0> h_mu;           // Population mean smoothing parameter
  
  // Components/weights
  positive_ordered[K] lambda;        // Eigenvalues
  vector<lower=0>[K] H;              // EF Smoothing parameters
  matrix[Q, K] X;                    // Unconstrained EF weights (X)
  matrix[N, K] Scores;               // EF scores (xi)
}

transformed parameters{
  // Fixed effects - population mean
  vector[M] mu = B * w_mu;
  
  // Eigenvalues
  vector[K] eval = reverse(lambda);
  
  // Orthogonal basis weights
  matrix[Q, K] Psi;
  
  // Polar decomposition
  {
    matrix[K,K] evec_XtX = eigenvectors_sym(crossprod(X)); 
    vector[K] eval_XtX = eigenvalues_sym(crossprod(X));
    Psi = X*evec_XtX*diag_matrix(1/sqrt(eval_XtX))*evec_XtX';
  }
}

model {
  // Smoothing weight priors 
  H ~ gamma(1, 0.005);
  h_mu ~ gamma(1, 0.005); 
  
   // Eigenvalue priors
  lambda ~ inv_gamma(0.1, 0.001);
  
  // Error
  sigma2 ~ inv_gamma(1, 0.001);
  
  // Smoothing additions to the target density
  target += -(1 * h_mu) / (2 * sigma2) * w_mu' * P * w_mu;
  for(i in 1:K){
    target += -(1 * H[i]) / (2 * sigma2) * Psi[,i]' * P * Psi[,i]; 
  }
  
  // Uniform priors through matrix normals
  to_vector(X) ~ normal(0, 1);
  
  // Score priors 
  for(i in 1:K){
    to_vector(Scores[,i]) ~ normal(0, sqrt(eval[i]));
  }
  
  // Likelihood
  {
    // Smooth component estimates - not saved as they can be calculated
    matrix[N, M] Theta = Scores * (B * Psi)';
    
    for(i in 1:N){
      Y[i] ~ normal(mu + to_vector(Theta[i,]), sqrt(sigma2));
    }
  }
}
