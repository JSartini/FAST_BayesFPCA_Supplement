fast_scores <- function(Scores, Order, Flips){
  n_samples = dim(Scores)[1]
  nEF = ncol(Flips)
  
  scores_df = map(1:n_samples, function(x){
    score_sample = Scores[x,,1:nEF]
    # Perform reorder
    if(nEF > 1){
      score_sample = score_sample[,Order[x,]]
    }
    
    # Perform flips
    for(j in 1:nEF){
      if(Flips[x,j]){
        score_sample[,j] = -score_sample[,j]
      }
    }
    
    # Collate into dataframe
    out_df = data.frame(score_sample)
    colnames(out_df) = paste0("EF ", 1:nEF)
    out_df$Sample = x
    out_df$Curve = 1:nrow(out_df)
    out_df = out_df %>%
      pivot_longer(-c(Sample, Curve), names_to = "EFNum", values_to = "Score")
    return(out_df)
  }) %>% list_rbind()
  
  return(scores_df)
}

GAM_scores_multilevel <- function(fitted, n_persons, n_ts, nEF1, nEF2){
  coefficients = coef(fitted)
  
  scores_level1 = matrix(0, nrow = n_persons, ncol = nEF1)
  for(i in 1:nEF1){
    scores_level1[,i] = coefficients[((i-1)*n_persons + 1):(i*n_persons)]
  }
  
  scores_level2 = matrix(0, nrow = n_ts, ncol = nEF2)
  after_level1 = nEF1*n_persons
  for(i in 1:nEF2){
    scores_level2[,i] = coefficients[((i-1)*n_ts + after_level1+1):(i*n_ts + after_level1)]
  }
  return(list(Lvl1 = scores_level1, Lvl2 = scores_level2))
}

fit_GAM_multilevel <- function(Y, EF_1, EF_2, mu, IDs){
  n_persons = n_distinct(IDs)
  n_ts = nrow(Y)
  
  pop_matrix = rep(mu, n_ts) %>% matrix(byrow = T, nrow = n_ts)
  
  EF1 = data.frame(EF_1)
  nEF1 = ncol(EF1)
  colnames(EF1) = paste0("Phi_ID_", 1:nEF1)
  p = nrow(EF1)
  EF1$Arg = as.numeric(1:p)
  
  EF2 = data.frame(EF_2)
  nEF2 = ncol(EF2)
  colnames(EF2) = paste0("Phi_PM_", 1:nEF2)
  EF2$Arg = as.numeric(1:p)
  
  Yp = Y - pop_matrix
  colnames(Yp) = as.character(1:p)
  data_for_fit = as.data.frame(Yp)
  data_for_fit$ID = as.factor(IDs)
  data_for_fit$PM = as.factor(1:nrow(data_for_fit))
  data_for_fit = data_for_fit %>%
    pivot_longer(-c(ID, PM), names_to = "Arg", values_to = "Yp") %>%
    mutate(Arg = as.numeric(Arg)) %>%
    inner_join(EF1, by = "Arg") %>%
    inner_join(EF2, by = "Arg")
  
  form_ID = paste0("s(ID, by = Phi_ID_", 1:nEF1, ", bs = 're')", collapse = "+")
  form_PM = paste0("s(PM, by = Phi_PM_", 1:nEF2, ", bs = 're')", collapse = "+")
  form = paste0("Yp ~ ", form_ID, "+", form_PM, "-1")
  
  model = bam(formula(form), method = "fREML", data = data_for_fit, discrete = T)
  return(list(fit.df = data_for_fit, fit = model))
}

gfsr_scores_ML <- function(Original_Y, EF_matrix1, EF_matrix2, beta, B, nEF1,
                           nEF2, IDs, Flips1, Flips2){
  n_samples = dim(EF_matrix1)[1]
  
  scores_df = map(1:n_samples, function(x){
    EF1 = EF_matrix1[x,,]
    EF2 = EF_matrix2[x,,]
    
    # Extract mu and eta estimates
    FE = beta[x,,] %*% t(B)
    mu = FE[1,] %>% as.numeric()
    
    # Project onto EF estimates
    gam_obj = fit_GAM_multilevel(Original_Y, EF1, EF2, mu, IDs)
    
    # Extract scores
    score_sample = GAM_scores_multilevel(gam_obj$fit, n_distinct(IDs), length(IDs), nEF1, nEF2)
    
    # Perform flips
    for(j in 1:nEF1){
      if(Flips1[x,j]){
        score_sample$Lvl1[,j] = -score_sample$Lvl1[,j]
      }
    }
    for(j in 1:nEF2){
      if(Flips2[x,j]){
        score_sample$Lvl2[,j] -score_sample$Lvl2[,j]
      }
    }
    
    # Collate into dataframes
    out_df1 = data.frame(score_sample$Lvl1)
    colnames(out_df1) = paste0("EF ", 1:nEF1)
    out_df1$Sample = x
    out_df1$Curve = 1:nrow(out_df1)
    out_df1 = out_df1 %>%
      pivot_longer(-c(Sample, Curve), names_to = "EFNum", values_to = "Score")
    
    out_df2 = data.frame(score_sample$Lvl2)
    colnames(out_df2) = paste0("EF ", 1:nEF2)
    out_df2$Sample = x
    out_df2$Curve = 1:nrow(out_df2)
    out_df2 = out_df2 %>%
      pivot_longer(-c(Sample, Curve), names_to = "EFNum", values_to = "Score")
    return(rbind(out_df1 %>% mutate(Level = 1), out_df2 %>% mutate(Level = 2)))
  }) %>% list_rbind()
  
  return(scores_df)
}

fit_GAM <- function(Y, EF, mu){
  n_ts = nrow(Y)
  pop_matrix = rep(mu, n_ts) %>% matrix(byrow = T, nrow = n_ts)
  
  EF_df = data.frame(EF)
  nEF = ncol(EF_df)
  colnames(EF_df) = paste0("Phi_", 1:nEF)
  p = nrow(EF_df)
  EF_df$Arg = 1:p
  
  Yp = Y - pop_matrix
  colnames(Yp) = as.character(1:p)
  data_for_fit = as.data.frame(Yp)
  data_for_fit$Curve = as.factor(1:n_ts)
  data_for_fit = data_for_fit %>%
    pivot_longer(-c(Curve), names_to = "Arg", values_to = "Yp") %>%
    mutate(Arg = as.numeric(Arg)) %>%
    inner_join(EF_df, by = "Arg")
  
  form_RE = paste0("s(Curve, by = Phi_", 1:nEF, ", bs = 're')", collapse = "+")
  form = paste0("Yp ~ ", form_RE, "-1")
  
  model = bam(formula(form), method = "fREML", data = data_for_fit, discrete = T)
  return(list(fit.df = data_for_fit, fit = model))
}

GAM_scores <- function(fitted, n_ts, nEF){
  coefficients = coef(fitted)
  
  scores = matrix(0, nrow = n_ts, ncol = nEF)
  for(i in 1:nEF){
    scores[,i] = coefficients[((i-1)*n_ts + 1):(i*n_ts)]
  }
  
  return(scores)
}

gfsr_scores <- function(Original_Y, EF_matrix, beta, B, Flips){
  n_samples = dim(EF_matrix)[1]
  nEF = ncol(Flips)
  
  scores_df = map(1:n_samples, function(x){
    EF_sample = EF_matrix[x,,]
    
    # Extract mu and eta estimates
    FE = beta[x,,] %*% t(B)
    mu = FE[1,] %>% as.numeric()
    
    # Project onto EF estimates
    gam_obj = fit_GAM(Original_Y, EF_sample, mu)
    
    # Extract scores
    score_sample = GAM_scores(gam_obj$fit, nrow(Original_Y), nEF)
    
    # Perform flips
    for(j in 1:nEF){
      if(Flips[x,j]){
        score_sample[,j] = -score_sample[,j]
      }
    }
    
    # Collate into dataframe
    out_df = data.frame(score_sample)
    colnames(out_df) = paste0("EF ", 1:nEF)
    out_df$Sample = x
    out_df$Curve = 1:nrow(out_df)
    out_df = out_df %>%
      pivot_longer(-c(Sample, Curve), names_to = "EFNum", values_to = "Score")
    return(out_df)
  }) %>% list_rbind()
  
  return(scores_df)
}

svdgp_scores <- function(Scores, Order, Flips){
  n_samples = dim(Scores)[1]
  nEF = ncol(Flips)
  
  scores_df = map(1:n_samples, function(x){
    score_sample = Scores[x,,1:nEF]
    # Perform reorder
    if(nEF > 1){
      score_sample = score_sample[,Order[x,]]
    }
    
    # Perform flips
    for(j in 1:nEF){
      if(Flips[x,j]){
        score_sample[,j] = -score_sample[,j]
      }
    }
    
    # Collate into dataframe
    out_df = data.frame(score_sample)
    colnames(out_df) = paste0("EF ", 1:nEF)
    out_df$Sample = x
    out_df$Curve = 1:nrow(out_df)
    out_df = out_df %>%
      pivot_longer(-c(Sample, Curve), names_to = "EFNum", values_to = "Score")
    return(out_df)
  }) %>% list_rbind()
  
  return(scores_df)
}


