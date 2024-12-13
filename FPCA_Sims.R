source("FPCA_Processing/Libs.R")
source("FPCA_Processing/Generate.R")
source("FPCA_Processing/Extract_FE.R")
source("FPCA_Processing/Extract_EF.R")
source("FPCA_Processing/Extract_Scores.R")

results_directory = "Results"

create_dir <- function(direct){
  if (!file.exists(direct)){
    dir.create(direct, recursive = T)
  }
}

gen_funcs = list(gen_dataset_CGM = gen_dataset_CGM,
                 gen_dataset_PA = gen_dataset_PA)

args = commandArgs(trailingOnly=TRUE)

n_ts = 50
n_obs = 30
n_sim = as.numeric(args[3])
generate_dataset = gen_funcs[[as.numeric(args[4])]]

message("Output directory: ", args[1])
message("Output file: ", args[2])
message("Number of simulations: ", args[3])
message("Generation function: ", names(gen_funcs)[as.numeric(args[4])])

FE_out = data.frame(Method = c(), Mean_Cov = c(), MSE = c(), Sample = c())
EF_out = data.frame(Method = c(), EFNum = c(), Mean_Cov = c(), MSE = c(), Sample = c())
Score_out = data.frame(Method = c(), EFNum = c(), Mean_cov = c(), MSE = c(), Sample = c())
Time_out = data.frame(Method = c(), Times = c(), Sample = c())

for(x in 1:n_sim){
  print(paste0("Iteration: ", x))
  
  # Generate the data
  sim_dataset = generate_dataset(n_ts, n_obs, spline_q = 10)
  true_EF = sim_dataset$True$EF
  
  # Fast BayesFPCA
  {
    # Seed values using FPCA estimates
    start.time = Sys.time()
    pre_mod = fpca.face(sim_dataset$True$Y, knots = floor(n_obs/2), 
                        npc = sim_dataset$direct$K, var = T)
    scaling = diag(t(pre_mod$efunctions) %*% pre_mod$efunctions)[1:sim_dataset$direct$K] 
    init_list = list(w_mu = coef(lm(pre_mod$mu ~ sim_dataset$direct$B - 1)),
                     lambda = rev(pre_mod$evalues[1:sim_dataset$direct$K]*scaling), 
                     sigma2 = pre_mod$sigma2) 
    fast_mod = stan(
      file = "STAN_Files/FPCA/FAST.stan",
      data = sim_dataset$direct, 
      chains = 4, 
      cores = 4, 
      warmup = 500, 
      iter = 1000, 
      init = list(init_list, init_list, 
                  init_list, init_list), 
      verbose = F,
      refresh = 0
    )
    # Extract samples
    samples = extract(fast_mod)

    # Mu samples
    fast_mu = fast_FE(samples$w_mu, sim_dataset$direct$B, sim_dataset$Domain)
    
    # EF samples
    preproc = fast_align(samples$Scores, samples$Psi, sim_dataset$direct$B, 
                              anchor = true_EF)
    fast_EF_objs = fast_EF(samples$Psi, sim_dataset$direct$B, sim_dataset$Domain, 
                           preproc$Ord, preproc$Flip)
    
    fast_EF_estimate = fast_EF_objs$Mean
    fast_EF_sample = fast_EF_objs$Samples
    
    # Score samples
    fast_S = fast_scores(samples$Scores, preproc$Ord, preproc$Flip)
    
    # Conclude
    end.time = Sys.time()
    fast_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  
  # GFSR framework: Goldsmith, Schrack, and Zipunnikov
  {
    start.time = Sys.time()
    gfsr_mod = stan(
      file = "STAN_Files/FPCA/GFSR.stan",
      data = sim_dataset$regress,
      chains = 4, 
      cores = 4,
      warmup = 500,
      iter = 1000, 
      verbose = F,
      refresh = 0
    )
    # Extract samples
    samples = extract(gfsr_mod)

    # Mu samples
    gfsr_mu = gfsr_FE(samples$beta, sim_dataset$regress$BS, sim_dataset$Domain)
    
    # EF samples
    preproc = gfsr_align(samples$beta_psi, sim_dataset$regress$BS, samples$c,
                             sim_dataset$regress$Kp, anchor = true_EF)
    under = gfsr_latent(samples, sim_dataset$regress$BS)
    
    gfsr_EF_objs = gfsr_EF(preproc$EF_mat, sim_dataset$Domain, preproc$Flip, under)
    
    gfsr_EF_estimate = gfsr_EF_objs$Mean
    gfsr_EF_sample = gfsr_EF_objs$Samples
    
    # Score samples
    gfsr_S = gfsr_scores(sim_dataset$True$Y, preproc$EF_mat, samples$beta, sim_dataset$regress$BS,
                        preproc$Flip)
    
    # Conclude
    end.time = Sys.time()
    gfsr_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  
  # SVD-GP: Jauch, Dunson, and Hoff model
  {
    start.time = Sys.time()
    svdgp_mod = stan(
      file = "STAN_Files/FPCA/SVD-GP.stan",
      data = sim_dataset$jdh,
      chains = 4, 
      cores = 4,
      warmup = 500,
      iter = 1000,
      verbose = F,
      refresh = 0
    )
    # Extract samples
    samples = extract(svdgp_mod)

    # Mu samples
    svdgp_mu = svdgp_FE(samples$mu_func, sim_dataset$Domain)
    
    # EF samples
    preproc = svdgp_align(samples$Scores, samples$V, anchor = true_EF) 
    svdgp_EF_objs = svdgp_EF(samples$V, sim_dataset$Domain, 
                             preproc$Ord, preproc$Flip)
    
    svdgp_EF_estimate = svdgp_EF_objs$Mean
    svdgp_EF_sample = svdgp_EF_objs$Samples
    
    # Score samples
    svdgp_S = svdgp_scores(samples$Scores, preproc$Ord, preproc$Flip)
    
    # Conclude
    end.time = Sys.time()
    svdgp_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  
  # Coverage and MSE of fixed effects function mu
  {
    trueMu_df = data.frame(Mu = sim_dataset$Mu, Arg = sim_dataset$Domain)
    
    Mu_df = rbind(fast_mu %>% mutate(Method = "FAST"),
                  gfsr_mu %>% mutate(Method = "GFSR"), 
                  svdgp_mu %>% mutate(Method = "SVD-GP"))
    
    FE_df = Mu_df %>%
      group_by(Arg, Method) %>%
      summarize(Estimate = mean(Mu),
                Lower = quantile(Mu, probs = c(0.025)),
                Upper = quantile(Mu, probs = c(0.975))) %>%
      ungroup() %>%
      left_join(trueMu_df, by = c("Arg")) %>%
      mutate(Coverage = case_when(Lower <= Mu & Upper >= Mu ~ 1,
                                  TRUE ~ 0), 
             Error = (Mu - Estimate)^2) %>%
      group_by(Method) %>%
      summarize(Mean_Cov = mean(Coverage), 
                MSE = mean(Error))
    FE_df$Sample = x
  }
  
  # Coverage and MSE of Eigenfunction estimates
  {
    trueEF_df = data.frame(sim_dataset$True$EF)
    colnames(trueEF_df) = paste0("EF ", 1:ncol(trueEF_df))
    trueEF_df$Arg = sim_dataset$Domain
    trueEF_df = trueEF_df %>%
      pivot_longer(-c(Arg), names_to = "EFNum", values_to = "EigFunc")
    
    Est_df = rbind(fast_EF_estimate %>% mutate(Method = "FAST"),
                   gfsr_EF_estimate %>% mutate(Method = "GFSR"), 
                   svdgp_EF_estimate %>% mutate(Method = "SVD-GP")) %>%
      left_join(trueEF_df, by = c("Arg", "EFNum")) %>%
      mutate(Error = (EigFunc - Estimate)^2) %>%
      group_by(Method, EFNum) %>%
      summarize(MSE = mean(Error)) 
    
    Cov_df = rbind(fast_EF_sample %>% mutate(Method = "FAST"), 
                   gfsr_EF_sample %>% mutate(Method = "GFSR"), 
                   svdgp_EF_sample %>% mutate(Method = "SVD-GP")) %>%
      group_by(Arg, EFNum, Method) %>%
      summarize(Lower = quantile(EigFunc, probs = c(0.025)),
                Upper = quantile(EigFunc, probs = c(0.975))) %>%
      ungroup() %>%
      left_join(trueEF_df, by = c("Arg", "EFNum")) %>%
      mutate(Coverage = case_when(Lower <= EigFunc & Upper >= EigFunc ~ 1,
                                  TRUE ~ 0)) %>%
      group_by(Method, EFNum) %>%
      summarize(Mean_Cov = mean(Coverage)) 
    
    EF_df = left_join(Est_df, Cov_df, by = c("Method", "EFNum"))
    EF_df$Sample = x
  }
  
  # Coverage and MSE of Eigenfunction Scores
  {
    trueScore_df = data.frame(sim_dataset$True$Scores)
    colnames(trueScore_df) = paste0("EF ", 1:ncol(trueScore_df))
    trueScore_df$Curve = 1:nrow(trueScore_df)
    trueScore_df = trueScore_df %>%
      pivot_longer(-c(Curve), names_to = "EFNum", values_to = "Score")
    
    score_df = rbind(fast_S %>% mutate(Method = "FAST"), 
                     gfsr_S %>% mutate(Method = "GFSR"), 
                     svdgp_S %>% mutate(Method = "SVD-GP"))
    
    score_df = score_df %>%
      group_by(Curve, EFNum, Method) %>%
      summarize(Estimate = mean(Score),
                Lower = quantile(Score, probs = c(0.025)),
                Upper = quantile(Score, probs = c(0.975))) %>%
      ungroup() %>%
      left_join(trueScore_df, by = c("Curve", "EFNum")) %>%
      mutate(Coverage = case_when(Lower <= Score & Upper >= Score ~ 1,
                                  TRUE ~ 0), 
             Error = (Score - Estimate)^2) %>%
      group_by(EFNum, Method) %>%
      summarize(Mean_Cov = mean(Coverage), 
                MSE = mean(Error)) 
    score_df$Sample = x
  }
  
  # Collated timing information
  {
    timing_df = data.frame(Times = c(fast_time, gfsr_time, svdgp_time) %>% as.numeric(), 
                           Method = c("FAST", "GFSR", "SVD-GP"))
    timing_df$Sample = x
  }
  
  # Add most recent result
  FE_out = rbind(FE_out, FE_df)
  EF_out = rbind(EF_out, EF_df)
  Score_out = rbind(Score_out, score_df)
  Time_out = rbind(Time_out, timing_df)
  
  # Write to storage
  FE_dir = paste0(results_directory, "/FE/", args[1])
  EF_dir = paste0(results_directory, "/EF/", args[1])
  Score_dir = paste0(results_directory, "/Score/", args[1])
  Time_dir = paste0(results_directory, "/Time/", args[1])
  
  create_dir(FE_dir)
  create_dir(EF_dir)
  create_dir(Score_dir)
  create_dir(Time_dir)
    
  write.csv(FE_out, paste0(FE_dir, "/", args[2], ".csv"))
  write.csv(EF_out, paste0(EF_dir, "/", args[2], ".csv"))
  write.csv(Score_out, paste0(Score_dir,  "/", args[2], ".csv"))
  write.csv(Time_out, paste0(Time_dir, "/", args[2], ".csv"))
}


