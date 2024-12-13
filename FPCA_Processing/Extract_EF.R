fast_align <- function(Scores, Psi, B, anchor = NULL){
  # Check ordering by eigenfunction
  score_vars = apply(Scores, c(1, 3), var)
  orders = do.call(rbind, map(1:nrow(score_vars), function(x){
    return(as.vector(sort(score_vars[x,], decreasing = T, index.return = T)$ix))
  }))
  
  # Extract eigenfunctions in matrix form
  n_samples = dim(Psi)[1]
  nEF = ncol(orders)
  
  EF_mat = map(1:n_samples, function(x){
    EF_sample = (B %*% Psi[x,,])[,1:nEF]
    if(nEF > 1){
      return(EF_sample[,orders[x,]])
    }
    else{
      return(EF_sample)
    }
  }) %>% abind(along = 0)
  
  # Derive negation flips by eigenfunction cross-correlation
  flips = matrix(F, nrow = n_samples, ncol = nEF)
  
  if(is.null(anchor)){
    for(i in 2:n_samples){
      for(j in 1:nEF){
        if(EF_mat[i,,j] %*% EF_mat[1,,j] < 0){
          flips[i,j] = T
        }
      }
    } 
  }
  else{
    for(i in 1:n_samples){
      for(j in 1:nEF){
        if(EF_mat[i,,j] %*% anchor[,j] < 0){
          flips[i,j] = T
        }
      }
    } 
  }
  
  return(list(Ord = orders, Flip = flips))
}

fast_EF <- function(Psi, B, domain, Order, Flips){
  n_samples = dim(Psi)[1]
  nEF = ncol(Flips)
  
  # Get all eigenfunction matrix samples in list
  EF_list = map(1:n_samples, function(x){
    EF_sample = (B %*% Psi[x,,])[,1:nEF]
    # Perform reorder
    if(nEF > 1){
      EF_sample = EF_sample[,Order[x,]]
    }
    
    # Perform flips
    for(j in 1:nEF){
      if(Flips[x,j]){
        EF_sample[,j] = -EF_sample[,j]
      }
    }
    return(EF_sample)
  })
  
  est_df = project_mean(EF_list, nEF, domain)
  
  EF_df = map(1:length(EF_list), function(x){
    # Collate into dataframe for CI calculations
    EF_sample = EF_list[[x]]
    out_df = data.frame(EF_sample)
    colnames(out_df) = paste0("EF ", 1:nEF)
    out_df$Sample = x
    out_df$Arg = domain
    out_df = out_df %>%
      pivot_longer(-c(Sample, Arg), names_to = "EFNum", values_to = "EigFunc")
    return(out_df)
  }) %>% list_rbind()
  
  return(list(Samples = EF_df, Mean = est_df))
}

gfsr_align <- function(Weights, B, Scores, nEF, anchor = NULL){
  n_samples = dim(Weights)[1]
  
  EF_mat = map(1:n_samples, function(x){
    UDV = Scores[x,,] %*% t(B %*% t(Weights[x,,]))
    EF_sample = svd(UDV, nv = nEF)$v
    return(EF_sample)
  }) %>% abind(along = 0)
  
  # Derive negation flips by eigenfunction cross-correlation
  flips = matrix(F, nrow = n_samples, ncol = nEF)
  
  if(is.null(anchor)){
    for(i in 2:n_samples){
      for(j in 1:nEF){
        if(EF_mat[i,,j] %*% EF_mat[1,,j] < 0){
          flips[i,j] = T
        }
      }
    } 
  }
  else{
    for(i in 1:n_samples){
      for(j in 1:nEF){
        if(EF_mat[i,,j] %*% anchor[,j] < 0){
          flips[i,j] = T
        }
      }
    } 
  }
  
  return(list(Flip = flips, EF_mat = EF_mat))
}

gfsr_latent <- function(samples, B, level = NULL){
  if(is.null(level)){
    weights = samples$c
    funcs = samples$beta_psi
    samp = dim(weights)[1]
    smooths = array(0, dim = c(samp, dim(weights)[2], dim(B)[1]))
    for(x in 1:samp){
      smooths[x,,] = weights[x,,] %*% funcs[x,,] %*% t(B)
    }
  }
  else{
    if(level == 1){
      weights = samples$c1
      funcs = samples$beta_psi1
      samp = dim(weights)[1]
      smooths = array(0, dim = c(samp, dim(weights)[2], dim(B)[1]))
      for(x in 1:samp){
        smooths[x,,] = weights[x,,] %*% funcs[x,,] %*% t(B) 
      }
    }
    else {
      weights = samples$c2
      funcs = samples$beta_psi2
      samp = dim(weights)[1]
      smooths = array(0, dim = c(samp, dim(weights)[2], dim(B)[1]))
      for(x in 1:samp){
        smooths[x,,] = weights[x,,] %*% funcs[x,,] %*% t(B) 
      }
    }
  }
  return(smooths)
}

gfsr_EF <- function(EF_matrix, domain, Flips, Smooths){
  n_samples = dim(EF_matrix)[1]
  nEF = ncol(Flips)
  
  EF_df = map(1:n_samples, function(x){
    EF_sample = EF_matrix[x,,]
    # Perform flips
    for(j in 1:nEF){
      if(Flips[x,j]){
        EF_sample[,j] = -EF_sample[,j]
      }
    }
    
    # Collate into dataframe
    out_df = data.frame(EF_sample)
    colnames(out_df) = paste0("EF ", 1:nEF)
    out_df$Sample = x
    out_df$Arg = domain
    out_df = out_df %>%
      pivot_longer(-c(Sample, Arg), names_to = "EFNum", values_to = "EigFunc")
    
    return(out_df)
  }) %>% list_rbind()
  
  svd_matrix = svd(apply(Smooths, c(2,3), mean), nv = nEF)$v
  svd_df = data.frame(svd_matrix)
  colnames(svd_df) = paste0("EF ", 1:nEF)
  svd_df$Arg = domain
  svd_df = svd_df %>%
    pivot_longer(-c(Arg), names_to = "EFNum", values_to = "Estimate")
  
  return(list(Samples = EF_df, Mean = svd_df))
}

svdgp_align <- function(Scores, V, anchor = NULL){
  # Check ordering by eigenfunction
  score_vars = apply(Scores, c(1, 3), var)
  orders = do.call(rbind, map(1:nrow(score_vars), function(x){
    return(as.vector(sort(score_vars[x,], decreasing = T, index.return = T)$ix))
  }))
  
  # Extract eigenfunctions in matrix form
  n_samples = dim(V)[1]
  nEF = ncol(orders)
  
  EF_mat = map(1:n_samples, function(x){
    EF_sample = V[x,,][,1:nEF]
    if(nEF > 1){
      return(EF_sample[,orders[x,]])
    }
    else{
      return(EF_sample)
    }
  }) %>% abind(along = 0)
  
  # Derive negation flips by eigenfunction cross-correlation
  flips = matrix(F, nrow = n_samples, ncol = nEF)
  
  if(is.null(anchor)){
    for(i in 2:n_samples){
      for(j in 1:nEF){
        if(EF_mat[i,,j] %*% EF_mat[1,,j] < 0){
          flips[i,j] = T
        }
      }
    } 
  }
  else{
    for(i in 1:n_samples){
      for(j in 1:nEF){
        if(EF_mat[i,,j] %*% anchor[,j] < 0){
          flips[i,j] = T
        }
      }
    } 
  }
  
  return(list(Ord = orders, Flip = flips))
}

svdgp_EF <- function(V, domain, Order, Flips){
  n_samples = dim(V)[1]
  nEF = ncol(Order)
  
  EF_list = map(1:n_samples, function(x){
    EF_sample = V[x,,]
    
    # Perform reorder
    if(nEF > 1){
      EF_sample = EF_sample[,Order[x,]]
    }
    
    # Perform flips
    for(j in 1:nEF){
      if(Flips[x,j]){
        EF_sample[,j] = -EF_sample[,j]
      }
    }
    
    return(EF_sample)
  })
  
  # Calculate manifold mean
  est_df = project_mean(EF_list, nEF, domain)
  
  EF_df = map(1:length(EF_list), function(x){
    # Collate into dataframe
    EF_sample = EF_list[[x]]
    out_df = data.frame(EF_sample)
    colnames(out_df) = paste0("EF ", 1:nEF)
    out_df$Sample = x
    out_df$Arg = domain
    out_df = out_df %>%
      pivot_longer(-c(Sample, Arg), names_to = "EFNum", values_to = "EigFunc")
    return(out_df)
  }) %>% list_rbind()
  
  return(list(Samples = EF_df, Mean = est_df))
}

# Calculate euclidean mean and project onto the manifold
project_mean <- function(EF_list, nEF, domain){
  EFs = abind(EF_list, along = 3)
  mean_EF = apply(EFs, c(1,2), mean)
  projection_comps = svd(mean_EF)
  projection = projection_comps$u %*% t(projection_comps$v)
  mm_df = data.frame(projection)
  colnames(mm_df) = paste0("EF ", 1:nEF)
  mm_df$Arg = domain
  mm_df = mm_df %>%
    pivot_longer(-c(Arg), names_to = "EFNum", values_to = "Estimate")
  return(mm_df)
}
