data{
  int<lower=1> n;         // Number of time series
  int<lower=1> p;         // Dimension of observation grid
  matrix[n,p] Y;          // Data matrix
  int<lower=1> k;         // Number of principal components
  real t[p];              // Locations of data in domain
  real<lower=0> nu;       // Hyperparameter of inverse gamma prior for sig2 
  real<lower=0> s2;       // Hyperparameter of inverse gamma prior for sig2
  real<lower=0> tau;      // SD parameter of prior for D
  real<lower=0> alpha;    // Hyperparameter of inverse gamma prior for smoothing parameters
  real<lower=0> beta;     // Hyperparameter of inverse gamma prior for smoothing parameters
}

parameters{
  vector[n*k] z_U; 
  vector[p*k] z_V; 
  vector<lower=0>[k] D; 
  vector[p] mu_func;       // fixed effect mu function
  real<lower=0> sig2;      // Marginal variance of autoregressive process
  real<lower=0> rho;       // Length scale parameter - eigenfunctions
  real<lower=0> lambda;    // Length scale parameter - mean function
}

transformed parameters{
  matrix[n,k] U; 
  matrix[p,k] V;
  vector[p] expect_mu;     // Expected prior values for mean function mu
  matrix[p,p] Sigma;       // Smoothing prior for fixed effects population mean
  matrix[n, k] Scores;
  {
    matrix[p,p] K;           // Parameter of MACG prior for V
    matrix[p,p] L_K;         // Cholesky decomposition of K
    matrix[p,k] X_V;         // V is the orth. component of the polar decomp. of X_V
    vector[k] eval_V;        // Eigenvalues of X_V'*X_V
    vector[k] eval_trans_V;  // Transformation of eigenvalues for polar decomp.
    matrix[k,k] evec_V;      // Eigenvectors of X_V'*X_V
    matrix[n,k] X_U;         // U is the orth. component of the polar decomp. of X_U
    vector[k] eval_U;        // Eigenvalues of X_U'*X_U
    vector[k] eval_trans_U;  // Transformation of eigenvalues for polar decomp.
    matrix[k,k] evec_U;      // Eigenvectors of X_U'*X_U
    
    // Constructing U as orthogonal component of the polar decomposition of X_U 
    // U = to_matrix(z_U, n, k); 
    X_U = to_matrix(z_U, n, k); 
    eval_U = eigenvalues_sym(X_U'*X_U);
    for(i in 1:k){
      eval_trans_U[i] = 1.0/sqrt(eval_U[i]);
    }
    evec_U = eigenvectors_sym(X_U'*X_U); 
    U = X_U*evec_U*diag_matrix(eval_trans_U)*evec_U'; 
  
    // Constructing V as orthogonal component of the polar decomposition of X_V 
    K = cov_exp_quad(t, 1, rho) + diag_matrix(rep_vector(1e-7, p));
    L_K = cholesky_decompose(K);
    X_V = L_K*to_matrix(z_V, p, k);
    eval_V = eigenvalues_sym(X_V'*X_V);
    for(i in 1:k){
      eval_trans_V[i] = 1.0/sqrt(eval_V[i]);
    }
    evec_V = eigenvectors_sym(X_V'*X_V); 
    V = X_V*evec_V*diag_matrix(eval_trans_V)*evec_V';
  }
  
  // Create score matrix
  Scores = U * diag_matrix(D);
  
  // Construct the expectation vector for GP mean function
  for(i in 1:p){
      expect_mu[i] = 0;
  }
  
  // Construct the covariance matrix for GP mean function
  Sigma = cov_exp_quad(t, 1, lambda) + diag_matrix(rep_vector(1e-7, p));
}

model{ 
  matrix[n,p] M;
  M = Scores*V'; 
  
  // Prior specifications
  z_V ~ normal(0, 1.0); 
  z_U ~ normal(0, 1.0); 
  D ~ normal(0, tau);
  sig2 ~ inv_gamma(nu/2.0, nu*s2/2.0);
  rho ~ inv_gamma(alpha, beta); 
  lambda ~ inv_gamma(alpha, beta);
  mu_func ~ multi_normal(expect_mu, Sigma);
  
  // Likelihood
  for(i in 1:n){
    Y[i] ~ normal(mu_func + to_vector(M[i]), sqrt(sig2));
  }
  
}
