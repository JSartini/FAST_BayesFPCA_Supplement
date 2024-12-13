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

args = commandArgs(trailingOnly=TRUE)

n_group = 50
n_per_group = 5
n_obs = 50
n_sim = as.numeric(args[3])

message("Output directory: ", args[1])
message("Output file: ", args[2])
message("Number of simulations: ", args[3])

FE_out = data.frame(Method = c(), Mean_Cov = c(), MSE = c(), Sample = c())
EF_out = data.frame(Method = c(), EFNum = c(), Mean_Cov = c(), MSE = c(), Sample = c())
Score_out = data.frame(Method = c(), EFNum = c(), Mean_cov = c(), MSE = c(), Sample = c())
Time_out = data.frame(Method = c(), Times = c(), Sample = c())

ncores = parallelly::availableCores()
message("I have ", ncores, " cores available")

nworkers = length(parallelly::availableWorkers())
message("I have ", nworkers, " workers available")

for(x in 1:n_sim){
  print(paste0("Iteration: ", x))
  
  # Generate the data
  sim_dataset = gen_dataset_multilevel(n_group, n_per_group, n_obs, Q = 15)
  sim_Y = sim_dataset$True$Y
  true_EF1 = sim_dataset$True$EF_lvl1
  true_EF2 = sim_dataset$True$EF_lvl2
  IDs = sim_dataset$IDs
  VisitIDs = rep(1:n_per_group, n_group)
  spline_q = sim_dataset$direct$Q
  orthog_splines = sim_dataset$direct$B
  n_eigen_lvl1 = sim_dataset$direct$K1
  n_eigen_lvl2 = sim_dataset$direct$K2
  
  # Direct sampling
  {
    # Seed values using MFPCA estimates
    start.time = Sys.time()
    
    pre_mod = mfpca.face(sim_Y, id = IDs, knots = spline_q, 
                         npc = max(n_eigen_lvl1, n_eigen_lvl2))
    mu_weights = coef(lm(pre_mod$mu ~ orthog_splines-1))
    scaling1 = diag(t(pre_mod$efunctions$level1) %*% pre_mod$efunctions$level1)[1:n_eigen_lvl1] 
    scaling2 = diag(t(pre_mod$efunctions$level2) %*% pre_mod$efunctions$level2)[1:n_eigen_lvl2]
    
    init_list = list(w_mu = mu_weights, 
                     sigma2 = pre_mod$sigma2,
                     lambda_1 = rev(pre_mod$evalues$level1[1:n_eigen_lvl1]*scaling1), 
                     lambda_2 = rev(pre_mod$evalues$level2[1:n_eigen_lvl2]*scaling2))
    
    fast_mod = stan(
      file = "STAN_Files/MFPCA/FAST.stan",
      data = sim_dataset$direct, 
      init = list(init_list, init_list, 
                  init_list, init_list), 
      chains = 4, 
      cores = 4, 
      warmup = 500, 
      iter = 1000, 
      control = list(max_treedepth = 12), 
      verbose = F,
      refresh = 0
    )
    
    samples = extract(fast_mod)
    
    # Extract fixed effects
    fast_mu = fast_FE(samples$w_mu, orthog_splines, sim_dataset$Domain)
    
    # Extract eigenfunctions
    preproc1 = fast_align(samples$Scores_1, samples$Psi_1, orthog_splines, 
                          anchor = true_EF1)
    preproc2 = fast_align(samples$Scores_2, samples$Psi_2, orthog_splines,
                          anchor = true_EF2)
    fast_EF1 = fast_EF(samples$Psi_1, orthog_splines, sim_dataset$Domain, 
                       preproc1$Ord, preproc1$Flip)
    fast_EF2 = fast_EF(samples$Psi_2, orthog_splines, sim_dataset$Domain, 
                       preproc2$Ord, preproc2$Flip)
    
    fast_EF_estimate = rbind(fast_EF1$Mean %>% mutate(Level = 1), 
                             fast_EF2$Mean %>% mutate(Level = 2))
    fast_EF_sample = rbind(fast_EF1$Samples %>% mutate(Level = 1), 
                           fast_EF2$Samples %>% mutate(Level = 2))
    
    # Extract scores
    fast_S1 = fast_scores(samples$Scores_1, preproc1$Ord, preproc1$Flip)
    fast_S2 = fast_scores(samples$Scores_2, preproc2$Ord, preproc2$Flip)
    fast_S = rbind(fast_S1 %>% mutate(Level = 1), 
                     fast_S2 %>% mutate(Level = 2))
    
    end.time = Sys.time()
    fast_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  
  # Regression-based model
  {
    start.time = Sys.time()
    gfsr_mod = stan(
      file = "STAN_Files/MFPCA/GFSR.stan",
      data = sim_dataset$regress,
      chains = 4, 
      cores = 4,
      warmup = 500, 
      iter = 1000,
      verbose = F,
      refresh = 0
    )
    samples = extract(gfsr_mod)
    
    # Extract fixed effects
    gfsr_mu = gfsr_FE(samples$beta, sim_dataset$regress$BS, sim_dataset$Domain)
    
    # Extract eigenfunctions
    preproc1 = gfsr_align(samples$beta_psi1, orthog_splines, samples$c1,
                          n_eigen_lvl1, anchor = true_EF1)
    preproc2 = gfsr_align(samples$beta_psi2, orthog_splines, samples$c2,
                          n_eigen_lvl2, anchor = true_EF2)
    
    under1 = gfsr_latent(samples, orthog_splines, 1)
    under2 = gfsr_latent(samples, orthog_splines, 2)
    
    gfsr_EF1 = gfsr_EF(preproc1$EF_mat, sim_dataset$Domain, preproc1$Flip, under1)
    gfsr_EF2 = gfsr_EF(preproc2$EF_mat, sim_dataset$Domain, preproc2$Flip, under2)
    
    gfsr_EF_estimate = rbind(gfsr_EF1$Mean %>% mutate(Level = 1), 
                             gfsr_EF2$Mean %>% mutate(Level = 2))
    gfsr_EF_sample = rbind(gfsr_EF1$Samples %>% mutate(Level = 1), 
                           gfsr_EF2$Samples %>% mutate(Level = 2))
    
    # Extract scores
    gfsr_S = gfsr_scores_ML(sim_dataset$True$Y, preproc1$EF_mat, preproc2$EF_mat, 
                            samples$beta, orthog_splines, n_eigen_lvl1, 
                            n_eigen_lvl2, IDs, preproc1$Flip, preproc2$Flip)
    
    end.time = Sys.time()
    gfsr_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  
  # Integrated coverage/error of fixed effects
  {
    trueMu_df = data.frame(Mu = sim_dataset$Mu, Arg = sim_dataset$Domain)
    
    Mu_df = rbind(fast_mu %>% mutate(Method = "FAST"), 
                  gfsr_mu %>% mutate(Method = "GFSR"))
    
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
  
  # Integrated coverage/error of EF
  {
    tEFDF1 = data.frame(true_EF1)
    colnames(tEFDF1) = paste0("EF ", 1:n_eigen_lvl1)
    tEFDF1$Arg = sim_dataset$Domain
    tEFDF1 = tEFDF1 %>%
      pivot_longer(-c(Arg), names_to = "EFNum", values_to = "EigFunc")
    
    tEFDF2 = data.frame(true_EF2)
    colnames(tEFDF2) = paste0("EF ", 1:n_eigen_lvl2)
    tEFDF2$Arg = sim_dataset$Domain
    tEFDF2 = tEFDF2 %>%
      pivot_longer(-c(Arg), names_to = "EFNum", values_to = "EigFunc")
    
    tEFDF = rbind(tEFDF1 %>% mutate(Level = 1), 
                  tEFDF2 %>% mutate(Level = 2))
    
    Est_df = rbind(fast_EF_estimate %>% mutate(Method = "FAST"), 
                   gfsr_EF_estimate %>% mutate(Method = "GFSR")) %>%
      left_join(tEFDF, by = c("Arg", "EFNum", "Level")) %>%
      mutate(Error = (EigFunc - Estimate)^2) %>%
      group_by(Method, Level, EFNum) %>%
      summarize(MSE = mean(Error))
    
    Cov_df = rbind(fast_EF_sample %>% mutate(Method = "FAST"), 
                   gfsr_EF_sample %>% mutate(Method = "GFSR")) %>%
      group_by(Arg, EFNum, Method, Level) %>%
      summarize(Lower = quantile(EigFunc, probs = c(0.025)),
                Upper = quantile(EigFunc, probs = c(0.975))) %>%
      ungroup() %>%
      left_join(tEFDF, by = c("Arg", "EFNum", "Level")) %>%
      mutate(Coverage = case_when(Lower <= EigFunc & Upper >= EigFunc ~ 1,
                                  TRUE ~ 0)) %>%
      group_by(Method, Level, EFNum) %>%
      summarize(Mean_Cov = mean(Coverage))
    
    EF_df = left_join(Est_df, Cov_df, by = c("Method", "EFNum", "Level"))
    EF_df$Sample = x
  }
  
  # Mean coverage/error for Scores
  {
    Score1_df = data.frame(sim_dataset$True$Scores_lvl1)
    colnames(Score1_df) = paste0("EF ", 1:n_eigen_lvl1)
    Score1_df$Curve = 1:nrow(Score1_df)
    Score1_df = Score1_df %>%
      pivot_longer(-c(Curve), names_to = "EFNum", values_to = "Score")
    
    Score2_df = data.frame(sim_dataset$True$Scores_lvl2)
    colnames(Score2_df) = paste0("EF ", 1:n_eigen_lvl2)
    Score2_df$Curve = 1:nrow(Score2_df)
    Score2_df = Score2_df %>%
      pivot_longer(-c(Curve), names_to = "EFNum", values_to = "Score")
    
    TrueScore = rbind(Score1_df %>% mutate(Level = 1), 
                      Score2_df %>% mutate(Level = 2))
    
    score_df = rbind(fast_S %>% mutate(Method = "FAST"), 
                     gfsr_S %>% mutate(Method = "GFSR"))
    
    score_df = score_df %>%
      group_by(Curve, EFNum, Method, Level) %>%
      summarize(Estimate = mean(Score),
                Lower = quantile(Score, probs = c(0.025)),
                Upper = quantile(Score, probs = c(0.975))) %>%
      ungroup() %>%
      left_join(TrueScore, by = c("Curve", "EFNum", "Level")) %>%
      mutate(Coverage = case_when(Lower <= Score & Upper >= Score ~ 1,
                                  TRUE ~ 0), 
             Error = (Score - Estimate)^2) %>%
      group_by(EFNum, Level, Method) %>%
      summarize(Mean_Cov = mean(Coverage), 
                MSE = mean(Error)) 
    score_df$Sample = x
  }
  
  # Timing information
  {
    timing_df = data.frame(Times = c(fast_time, gfsr_time) %>% as.numeric(), 
                           Method = c("Direct", "Regress"))
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


