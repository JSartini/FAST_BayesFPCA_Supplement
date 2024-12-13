data {
  int N;         // Total number of functional observations
  int I;         // Greatest observed ID
  int J;         // Greatest observed visit
  int ID[N];     // Index-based ID per function
  int Visit[N];  // Index-based visit per function
  int M;         // Number of observations in functional domain

  int Q;    // Dimension of spline bases
  int K1;   // Number of Eigenfunctions at group level
  int K2;   // Number of Eigenfunctions at visit level
  
  matrix[N, M] Y;     // Original data
  matrix[M, Q] B;     // Orthogonalized basis
  matrix[Q, Q] P;     // Penalty matrix for splines
}

parameters {
  real<lower=0> sigma2; // Error in observation
  
  // Fixed-effect components
  vector[Q] w_mu;       // B-spline weights for overall mean
  real<lower=0> H_mu;   // Smoothness parameter for overall mean
  
  // Person-specific components
  positive_ordered[K1] lambda_1;  // Variances of each eigen-function
  vector<lower=0>[K1] H_1;        // Smoothing multipliers
  matrix[Q, K1] X_1;              // Full-rank representation decomposed into group EF
  matrix[I, K1] Scores_1;         // Group scores
  
  // Person-visit components
  positive_ordered[K2] lambda_2;  // Variances of each eigen-function
  vector<lower=0>[K2] H_2;        // Smoothing multiplier
  matrix[Q, K2] X_2;              // Full-rank representation decomposed into visit EF
  matrix[N, K2] Scores_2;         // Instance scores
}

transformed parameters{
  // Fixed effects
  vector[M] mu = B * w_mu;
  
  // Eigenvalues
  vector[K1] eval_1 = reverse(lambda_1);
  vector[K2] eval_2 = reverse(lambda_2);
  
  // Orthogonal bases
  matrix[Q, K1] Psi_1;
  matrix[Q, K2] Psi_2;
  
  // Polar decomposition
  {
    vector[K1] trans_eval_V1;   // Inverse sqrt of eigenvalues of X_1'*X_1
    matrix[K1,K1] evec_V1;      // Eigenvectors of X_1'*X_1
    
    trans_eval_V1 = 1/sqrt(eigenvalues_sym(X_1'*X_1));
    evec_V1 = eigenvectors_sym(X_1'*X_1); 
    Psi_1 = X_1*evec_V1*diag_matrix(trans_eval_V1)*evec_V1';
    
    vector[K2] trans_eval_V2;   // Inverse sqrt of eigenvalues of X_2'*X_2
    matrix[K2,K2] evec_V2;      // Eigenvectors of X_2'*X_2
    
    trans_eval_V2 = 1/sqrt(eigenvalues_sym(X_2'*X_2));
    evec_V2 = eigenvectors_sym(X_2'*X_2); 
    Psi_2 = X_2*evec_V2*diag_matrix(trans_eval_V2)*evec_V2';
  }
}

model {
  // Smoothing weight priors 
  H_mu ~ gamma(1, 0.005); 
  H_1 ~ gamma(1, 0.005); 
  H_2 ~ gamma(1, 0.005); 
  
  // Smoothing additions to the target density
  target += -(1 * H_mu) / (2 * sigma2) * w_mu' * P * w_mu;
  
  for(i in 1:K1){
    target += -(1 * H_1[i]) / (2 * sigma2) * Psi_1[,i]' * P * Psi_1[,i];
  }
  
  for(i in 1:K2){
    target += -(1 * H_2[i]) / (2 * sigma2) * Psi_2[,i]' * P * Psi_2[,i];
  } 
  
  // Uniform priors through matrix normals
  to_vector(X_1) ~ normal(0, 1);
  to_vector(X_2) ~ normal(0, 1);
  
  // Score priors 
  for(i in 1:K1){
    to_vector(Scores_1[,i]) ~ normal(0, sqrt(eval_1[i]));
  }
  for(i in 1:K2){
    to_vector(Scores_2[,i]) ~ normal(0, sqrt(eval_2[i]));
  }
  
  // Eigenvalue priors
  lambda_1 ~ inv_gamma(0.1, 0.001);
  lambda_2 ~ inv_gamma(0.1, 0.001);
  
  // Error
  sigma2 ~ inv_gamma(1, 0.001);
  
  // Likelihood 
  {
    // Smooth estimates - not saved as they can be calculated
    matrix[I, M] UDV_1 = Scores_1 * (B * Psi_1)';
    matrix[N, M] UDV_2 = Scores_2 * (B * Psi_2)';
    
    for(i in 1:N){
      Y[i,] ~ normal(mu + to_vector(UDV_1[ID[i],] + UDV_2[i,]), sqrt(sigma2));
    }
  }
}
