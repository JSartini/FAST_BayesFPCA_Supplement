
# Legendre polynomials of orders 0-5 on [0,1]
lp0 <- function(x){
  v = rep(1, length(x))
  return(v)
}

lp1 <- function(x){
  v = 2*x-1
  return(v)
}

lp2 <- function(x){
  v = 2*x-1
  return((3*v^2 - 1)/2)
}

lp3 <- function(x){
  v = 2*x-1
  return((5*v^3 - 3*v)/2)
}

lp4 <- function(x){
  v = 2*x-1
  return((35*v^4 - 30*v^2 + 3)/8)
}

lp5 <- function(x){
  v = 2*x-1
  return((63*v^5 - 70*v^3 + 15*v)/8)
}

# Generate data adhering to FPCA model (mean function provided)
gen_FPCA <- function(V, Lambdas, domain, N, sigma_eps, mu_func){
  nobs = length(domain)
  nf = ncol(V)
  
  # Store outputs
  outcome_data = matrix(0, nrow = N, ncol = nobs)
  
  Scores = matrix(0, nrow = N, ncol = nf)
  
  # Build eigen-bases and sample scores
  for(i in 1:nf){
    Scores[,i] = rnorm(N, sd = sqrt(Lambdas[i]))
  }
  
  # Build outputs
  smooth_data = Scores %*% t(V) + matrix(rep(mu_func, N), nrow = N, byrow = T)
  
  # Add noise
  outcome_data = smooth_data + matrix(rnorm(N*nobs, sd = sigma_eps), nrow = N, ncol = nobs)
  
  return(list(Y = outcome_data, True_Y = smooth_data, Scores = Scores, EF = V))
}

# Prior for Jauch, Dunson, and Hoff
JDH_priors <- function(crossings, ratio, K, Y, N, p, domain){
  expectation = 1/(2*pi*crossings)
  stddev = expectation*ratio
  
  solve_fn <- function(x){
    E <- x[2]/(x[1] - 1) - expectation
    V <- x[2]^2/((x[1] - 1)^2 * (x[1] - 2)) - stddev^2
    return(c(E, V))
  }
  
  smooth_priors = nleqslv(c(5, 5), solve_fn)$x
  alpha = smooth_priors[1]
  beta = smooth_priors[2]
  
  approx_comps = svd(Y, nv = K)
  best_approx = approx_comps$u[,1:K] %*% diag(approx_comps$d[1:K]) %*% t(approx_comps$v)
  nu = 1
  s2 = 3*var(as.numeric(Y - best_approx))
  
  tau2 = sum(diag(t(best_approx) %*% best_approx))/K
  return(list(n = N, p = p, Y = Y, k = K, t = domain, nu = nu,
              s2 = s2, tau = sqrt(tau2), alpha = alpha, beta = beta))
}

# Two realistic simulation scenarios based on real-data

# Based upon postprandial continuous glucose monitoring from DASH4D
gen_dataset_CGM <- function(n_ts, n_obs, fit_EF = 3, spline_q = 15){
  # Basic constants
  domain = seq(0, 1, length.out = n_obs)
  n_EF = 3
  
  # True$Y underlying basis functions
  V = cbind(lp0(domain) * sqrt(1/n_obs),
            (lp1(domain) - 0.5*lp3(domain)) * sqrt(84/(31 * n_obs)), 
            -1*(lp2(domain) * sqrt(5/n_obs)))
 
  
  # Eigenvalues and measurement noise based upon exploration
  lambdas = c(45000, 7000, 2000)
  sim_sigma = 2
  
  # Simulate hierarchical functional data
  mu_func = 145 - 10*lp2(domain)
  gen_data = gen_FPCA(V, lambdas, domain, n_ts, 
                      sigma = sqrt(sim_sigma), mu_func)
  
  # Generate spline basis with derivatives
  splines = bSpline(domain, df = spline_q, intercept = T)
  spline_derivs = deriv(splines, 2)
  
  # Orthogonalize with eigen-decomposition
  orthog_splines = eigen(splines %*% t(splines))$vectors[,1:spline_q]
  basis_svd = svd(splines)
  pi_splines = basis_svd$v %*% diag(basis_svd$d^-1) %*% t(basis_svd$u)
  transform = pi_splines %*% orthog_splines
  orthog_derivs = spline_derivs %*% transform
  
  # FE matrix - just a population intercept
  X_mat = matrix(0, nrow = n_ts, ncol = 1)
  X_mat[,1] = 1
  
  # Penalty Matrix
  p0 = (1/n_obs) * t(orthog_splines) %*% orthog_splines
  p2 = (1/n_obs) * t(orthog_derivs) %*% orthog_derivs
  PM = 0.1*p0 + 0.9*p2
  
  # Data for direct sampling model
  direct_dat = list(N = n_ts, M = n_obs, Q = spline_q, 
                    K = fit_EF, Y = gen_data$Y, 
                    B = orthog_splines, P = PM)
  
  # Data for regression model
  regress_dat = list(N = n_ts, D = n_obs, p = 1, 
                     Kt = spline_q, Kp = fit_EF, Y = gen_data$Y,
                     X = X_mat, BS = orthog_splines, PenMat = PM)
  
  # Data for conditional scores model
  cond_fit = fpca.face(gen_data$Y, knots = spline_q, npc = fit_EF)
  scores_dat = list(N = n_ts, p = n_obs, K = fit_EF)
  
  # Data for JDH model
  jdh_dat = JDH_priors(2, 1/6, fit_EF, gen_data$Y, n_ts, n_obs, domain)
  
  return(list(direct = direct_dat, regress = regress_dat,
              scores = scores_dat, jdh = jdh_dat, 
              Domain = domain, True = gen_data, Mu = mu_func, 
              fpca_mu = cond_fit$mu))
}

# Based upon actigraphy from NHANES
gen_dataset_PA <- function(n_ts, n_obs, fit_EF = 2, spline_q = 15){
  # Basic constants
  domain = seq(0, 1, length.out = n_obs)
  
  # True$Y underlying basis functions
  n_EF = 2
  V = matrix(0, ncol = n_EF, nrow = n_obs)
  V[,1] = (1 - 2*cos(2 * pi * domain)) * sqrt(1/(3 * n_obs))
  V[,2] = (5 * sin(2 * pi * domain) - sin(4 * pi * domain) - 2 * cos(4 * pi * domain)) * sqrt(1/(15 * n_obs))
  
  # Eigenvalues and measurement noise based upon exploration
  lambdas = c(25000, 8500)
  sim_sigma = 11
  
  # Simulate hierarchical functional data
  mu_func = 10 - 5*sin(2 * pi * domain) - 5*cos(2 * pi * domain)
  gen_data = gen_FPCA(V, lambdas, domain, n_ts, 
                      sigma = sqrt(sim_sigma), mu_func)
  
  # Generate spline basis with derivatives
  splines = bSpline(domain, df = spline_q, intercept = T)
  spline_derivs = deriv(splines, 2)
  
  # Orthogonalize with eigen-decomposition
  orthog_splines = eigen(splines %*% t(splines))$vectors[,1:spline_q]
  basis_svd = svd(splines)
  pi_splines = basis_svd$v %*% diag(basis_svd$d^-1) %*% t(basis_svd$u)
  transform = pi_splines %*% orthog_splines
  orthog_derivs = spline_derivs %*% transform
  
  # FE matrix - just a population intercept
  X_mat = matrix(0, nrow = n_ts, ncol = 1)
  X_mat[,1] = 1
  
  # Penalty Matrix
  p0 = (1/n_obs) * t(orthog_splines) %*% orthog_splines
  p2 = (1/n_obs) * t(orthog_derivs) %*% orthog_derivs
  PM = 0.1*p0 + 0.9*p2
  
  # Data for direct sampling model
  direct_dat = list(N = n_ts, M = n_obs, Q = spline_q, 
                    K = fit_EF, Y = gen_data$Y, 
                    B = orthog_splines, P = PM)
  
  # Data for regression model
  regress_dat = list(N = n_ts, D = n_obs, p = 1, 
                     Kt = spline_q, Kp = fit_EF, Y = gen_data$Y,
                     X = X_mat, BS = orthog_splines, PenMat = PM)
  
  # Data for conditional scores model
  cond_fit = fpca.face(gen_data$Y, knots = spline_q, npc = fit_EF)
  scores_dat = list(N = n_ts, p = n_obs, K = fit_EF)
  
  # Data for JDH model
  jdh_dat = JDH_priors(2, 1/6, fit_EF, gen_data$Y, n_ts, n_obs, domain)
  
  return(list(direct = direct_dat, regress = regress_dat,
              scores = scores_dat, jdh = jdh_dat, 
              Domain = domain, True = gen_data, Mu = mu_func, 
              fpca_mu = cond_fit$mu))
}

gen_MFPCA <- function(V1, V2, Lambda1, Lambda2, domain, N, ID, sigma_eps, mu_func){
  nobs = length(domain)
  nIf = ncol(V1)
  nPVf = ncol(V2)
  npersons = n_distinct(ID)
  
  # Store outputs
  outcome_data = matrix(0, nrow = N, ncol = nobs)
  
  # Sample scores
  Scores1 = matrix(0, nrow = n_distinct(ID), ncol = nIf)
  for(i in 1:nIf){
    Scores1[,i] = rnorm(npersons, sd = sqrt(Lambda1[i]))
  }
  
  Scores2 = matrix(0, nrow = N, ncol = nPVf)
  for(i in 1:nPVf){
    Scores2[,i] = rnorm(N, sd = sqrt(Lambda2[i]))
  }
  
  # Build outputs
  person_matrix = Scores1 %*% t(V1)
  pv_matrix = Scores2 %*% t(V2)
  
  # Combine
  outcome_data = pv_matrix
  for(i in 1:N){
    outcome_data[i, ] = outcome_data[i, ] + person_matrix[ID[i],] + mu_func
  }
  
  # Add noise
  noise_matrix = matrix(rnorm(N*nobs, sd = sigma_eps), nrow = N, ncol = nobs)
  outcome_data = outcome_data + noise_matrix
  
  return(list(Y = outcome_data, True_Y = outcome_data - noise_matrix, 
              person_matrix = person_matrix, visit_matrix = pv_matrix,
              EF_lvl1 = V1, Scores_lvl1 = Scores1, 
              EF_lvl2 = V2, Scores_lvl2 = Scores2))
}

gen_dataset_multilevel <- function(n_group, n_per_group, n_obs, Q = 10){
  # Basic constants
  n_ts = n_group * n_per_group
  IDs = rep(1:n_group, each = n_per_group)
  domain = seq(0, 1, length.out = n_obs)
  
  # True underlying basis functions
  mu_func = 5*lp1(domain) - 3*lp2(domain)
  V1 = cbind(sin(2*pi*domain) * sqrt(2/n_obs), cos(2*pi*domain) * sqrt(2/n_obs))
  V2 = cbind(sin(4*pi*domain) * sqrt(2/n_obs), cos(4*pi*domain) * sqrt(2/n_obs))
  
  # Eigenvalues and measurement noise based upon exploration
  Lambda1 = c(10000, 2000)
  Lambda2 = c(8000, 4000)
  sim_sigma = 4
  
  # Simulate hierarchical functional data
  gen_data = gen_MFPCA(V1, V2, Lambda1, Lambda2, domain, n_ts, IDs, sim_sigma, mu_func)
  
  # Generate spline basis with derivatives
  spline_q = Q
  splines = bSpline(domain, df = spline_q, intercept = T)
  spline_derivs = deriv(splines, 2)
  
  # Orthogonalize with eigen-decomposition
  orthog_splines = eigen(splines %*% t(splines))$vectors[,1:spline_q]
  basis_svd = svd(splines)
  pi_splines = basis_svd$v %*% diag(basis_svd$d^-1) %*% t(basis_svd$u)
  transform = pi_splines %*% orthog_splines
  orthog_derivs = spline_derivs %*% transform
  
  # FE matrix - just a population intercept
  X_mat = matrix(0, nrow = n_ts, ncol = 1)
  X_mat[,1] = 1
  
  # Penalty Matrix
  p0 = (1/n_obs) * t(orthog_splines) %*% orthog_splines
  p2 = (1/n_obs) * t(orthog_derivs) %*% orthog_derivs
  PM = 0.1*p0 + 0.9*p2
  
  # Data for direct sampling model
  direct_dat = list(N = n_ts, I = n_group, J = n_per_group, 
                    ID = IDs, Visit = rep(1:n_per_group, n_group), 
                    M = n_obs, Q = spline_q, K1 = ncol(V1), K2 = ncol(V2),
                    Y = gen_data$Y, B = orthog_splines, P = PM)
  
  # Data for regression model
  regress_dat = list(I = n_group, J = n_per_group, IJ = n_ts, D = n_obs, p = 1, 
                     Kt = spline_q, Kp1 = ncol(V1), Kp2 = ncol(V2),
                     Y = gen_data$Y, X = X_mat, BS = orthog_splines, 
                     PenMat = PM)
  
  # Data for conditional scores model
  cond_fit = mfpca.face(gen_data$Y, id = IDs, twoway = F,
                        knots = spline_q, npc = max(ncol(V1), ncol(V2)))
  scores_dat = list(N = n_ts, L = n_group, p = n_obs, 
                    kI = ncol(V1), kPV = ncol(V2), Y = gen_data$Y, 
                    ID = IDs, EF1 = cond_fit$efunctions$level1[,1:ncol(V1)], 
                    EF2 = cond_fit$efunctions$level2[,1:ncol(V2)])
  
  return(list(direct = direct_dat, regress = regress_dat, IDs = IDs,
              scores = scores_dat, Domain = domain, True = gen_data, 
              Mu = mu_func, init_mfpca = cond_fit))
}


