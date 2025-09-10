initialize_mixing_matrices <- function(Y, init_S, Nstar){
  
  # Get the number of components
  Q <- min(dim(init_S))
  V <- max(dim(init_S))
  if (nrow(init_S) != Q){ init_S <- t(init_S) }
  
  # Calculate projection matrix
  proj <- t(init_S) %*% solve(init_S %*% t(init_S))
  
  # Project each subject's data
  init_Ai <- matrix(0, nrow= Q, ncol = Q*Nstar)
  col_inds <- 1:Q
  new_init_S <- matrix(0, nrow = Q, ncol = V)
  Yinds <- 1:V
  for (i in 1:Nstar){
    
    Yi <- Y[, Yinds]
    
    # See if Y needs to be transposed
    #if (nrow(Yi) != Q) {Yi <- t(Yi)}
    
    # Get the singular value decomposition
    Ai_svd <- svd(Yi %*% proj)
    
    # Obtain the guess for Ai
    Ai <- Ai_svd$u %*% t(Ai_svd$v)
    
    # Store
    init_Ai[, col_inds] <- Ai
    
    # Use to work towards population map as average of unmixed data
    new_init_S <- new_init_S + t(Ai) %*% Yi / Nstar
    
    # Increment column counter
    col_inds <- col_inds + Q
    Yinds    <- Yinds + V
  }
  
  # S should be V x Q at output time
  
  return(list(
    init_S = t(new_init_S),
    init_A = init_Ai
  ))
    

}

# Function to obtain initial guess for the Dirichlet process mixture parameters
obtain_S0_DPM_guess <- function(S, Hmax){
  
  V <- nrow(S)
  Q <- ncol(S)
  
  S_flat <- c(S)
  
  # K-means
  kmeans_result <- kmeans(S_flat, centers = Hmax)
  
  mu_h       <- as.vector(kmeans_result$centers)
  Z_qv       <- matrix(kmeans_result$cluster, nrow = V)
  sigma_h_sq <- rep(0, Hmax)
  for (h in 1:Hmax){
    sigma_h_sq[h] <- var(S[Z_qv == h])
  }
  
  return(list(
    mu_h = mu_h,
    Z_qv = Z_qv,
    sigma_h_sq = sigma_h_sq
  ))
  
}




# Function to obtain initial guess for the Dirichlet process mixture parameters
# for the random effects covariances. This current implementation obtains a CRUDE guess
# as follows:
# (1) Subject/session specific networks are calculated as Aij'Yij
# (2) The current guess for the population average and betas is subtracted out
# (3) Using the residuals from step 2, we regress them onto the column space of Z
# (4) 1/Nstar * sum BB' is evaluated for each network and location
# (5) The matrices from step 4 are clustered using K-means clustering
# These clusters are used as the initial guess and weights
obtain_Rqv_DPM_guess <- function(Y_mixed_QxVNstar, A, S, Z, Prand, Hmax, Beta = NULL, X = NULL){
  
  V     <- nrow(S)
  Q     <- ncol(S)
  Nstar <- nrow(Z)

  # Check if we need to include covariate adjustment  
  include_covariate_adjustment <- FALSE
  if (!is.null(Beta) || !is.null(X)){
    if (!is.null(Beta) & !is.null(X)){
      include_covariate_adjustment <- TRUE
    } else {
      stop("Only one of Beta and X argument were required, but both are required to adjust initial Rqv DPM.")
    }
  }
  
  # Note, this X should not include the intercept, otherwise need to remove it
  covariate_effect_est <- NULL
  if (include_covariate_adjustment == TRUE){
    if (ncol(X) != nrow(Beta)){
      if (ncol(X) == (nrow(Beta) + 1)){
        X <- X[, -1] # remove intercept column
      } else {
        stop("Dimensions of X do not match Beta in Rqv initialization")
      }
    }
    covariate_effect_est   <- X %*% Beta # Nstar x QV
  }
  
  # Step 1: Subject/Session Networks
  Sij_crude <- matrix(0, nrow = Nstar, ncol = length(Y_mixed_QxVNstar) / Nstar)
  apply_unmixing_matrices(Sij_crude, Y_mixed_QxVNstar, A)
  
  # Step 2: Remove current values for the population networks and fixed effects
  Sij_residuals <- Sij_crude
  for (iq in 1:Q){
    index_set <- ((iq-1)*V+1):(iq*V)
    Sij_residuals[, index_set] <- sweep(Sij_crude[,index_set], 2, t(S[, iq]))
  }
  if (include_covariate_adjustment == TRUE){
    Sij_residuals = Sij_residuals - covariate_effect_est
  }
  
  # Step 3: Regress residuals onto column space of Z
  proj <- solve(t(Z) %*% Z) %*% t(Z) # (N x Prand) x Nstar
  # TODO - might be more storage efficient to do these steps in a loop
  Bhat <- proj %*% Sij_residuals # (N x Prand) x (Q x V)
  
  # Step 4: Estimate variance-covariance at each network/location
  # only keeping the diagonal + upper triangle for now
  covariance_utri      <- matrix(0, nrow = Prand + Prand * (Prand - 1) / 2, ncol = Q*V)
  for (icol in 1:(Q*V)){
    Bqv     <- matrix(Bhat[, icol], nrow = Prand)
    cov_hat <- cov(t(Bqv))
    covariance_utri[, icol] <- cov_hat[upper.tri(cov_hat, diag = TRUE)]
  }
  
  # Step 5: K-means clustering
  kmeans_result <- kmeans(t(covariance_utri), centers = Hmax)
  
  # Initialize DPM parameters
  Z_Rqv       <- matrix(kmeans_result$cluster, nrow = V)
  E_Rqv_inv   <- matrix(0, nrow = Prand, ncol = Prand * Hmax)
  index_set <- 1:Prand
  for (h in 1:Hmax){
    # Pull out the centroid
    cluster_utri <- kmeans_result$center[h,]
    # Turn into Prand x Prand matrix
    cluster_variance <- matrix(0, nrow = Prand, ncol = Prand)
    cluster_variance[upper.tri(cluster_variance, diag=TRUE)] <- cluster_utri
    cluster_variance[lower.tri(cluster_variance)] <- t(cluster_variance)[lower.tri(cluster_variance)]
    
    # Now store the precision
    E_Rqv_inv[, index_set] <- solve(cluster_variance)
    index_set <- index_set + Prand
  }
  
  return(list(
    Z_Rqv = Z_Rqv,
    E_Rqv_inv = E_Rqv_inv
  ))
  
}