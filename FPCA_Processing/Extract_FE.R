fast_FE <- function(weights, B, domain){
  n_samples = dim(weights)[1]
  
  EF_df = map(1:n_samples, function(x){
    mu_sample = B %*% weights[x,]
    
    # Collate into dataframe
    FE_df = data.frame(Mu = mu_sample, Arg = domain,
                        Sample = x)
    return(FE_df)
  }) %>% list_rbind()
  
  return(EF_df)
}

gfsr_FE <- function(weights, B, domain){
  n_samples = dim(weights)[1]
  
  FE_df = map(1:n_samples, function(x){
    # FE structure requires additional indexing
    mu_sample = as.numeric(B %*% weights[x,,]) 
    
    FE_df = data.frame(Mu = mu_sample, Arg = domain,
                       Sample = x)
    return(FE_df)
  }) %>% list_rbind()
  
  return(FE_df)
}

svdgp_FE <- function(FE_samples, domain){
  n_samples = dim(FE_samples)[1]
  
  FE_df = map(1:n_samples, function(x){
    mu_sample = FE_samples[x,]
    
    FE_df = data.frame(Mu = mu_sample, Arg = domain,
                       Sample = x)
    return(FE_df)
  }) %>% list_rbind()
  
  return(FE_df)
}