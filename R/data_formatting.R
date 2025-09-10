format_list_time_courses_as_matrix <- function(Y){
  
  # Number of subjects (may also be subjects x visits)
  N <- length(Y)
  
  # Make sure all are the same dimension
  row_counts <- unlist(lapply(Y, nrow))
  col_counts <- unlist(lapply(Y, ncol))
  
  if (length(unique(row_counts)) != 1){
    stop("Number of data rows is not consistent across subjects.")
  }
  
  if (length(unique(col_counts)) != 1){
    stop("Number of data columns is not consistent across subjects.")
  }
  
  # Number of components
  Q <- row_counts[1]
  V <- col_counts[1]
  
  if (Q > V){
    Q <- col_counts[1]
    V <- row_counts[1]
  }
  
  # Storage for data
  data_QxNV <- matrix(0, nrow = Q, ncol = N*V)
  
  # Loop over subjects and store each subjects data
  for (i in 1:N){
    
    # Make sure Q x V
    Yi <- Y[[i]]
    if (nrow(Yi) != Q){
      Yi <- t(Yi)
    }
    
    data_QxNV[, ((i-1)*V+1):(i*V) ] <- Yi
  }
  
  return(data_QxNV)
}


format_array_time_courses_as_matrix <- function(Y){
  
  # Number of subjects (may also be subjects x visits)
  N <- dim(Y)[3]
  
  # Number of components
  Q <- dim(Y)[2]
  V <- dim(Y)[1]
  transpose_flag <- TRUE
  if (Q > V){
    #stop(paste0("Number of components, ", Q, " is greater than number of spatial locations, ", V))
    V <- dim(Y)[2]
    Q <- dim(Y)[1]
    transpose_flag <- FALSE
  }
  
  # Storage for data
  data_QxNV <- matrix(0, nrow = Q, ncol = N*V)
  
  # Loop over subjects and store each subjects data
  for (i in 1:N){
    if (transpose_flag == TRUE){
      data_QxNV[, ((i-1)*V+1):(i*V) ] <- t(Y[,,i])
    } else {
      data_QxNV[, ((i-1)*V+1):(i*V) ] <- Y[,,i]
    }
  }
  
  return(data_QxNV)
  
}

# if id_variable is provided (required for 3d array case) then can 
# check it against the Ji variable
format_time_series_input <- function(Y, longit_adjusted, 
                                     debug = FALSE,
                                     id_variable = NULL){
  
  # check how many unique IDs were provided, only needed for longitudinal
  # case
  n_id <- NULL
  if (longit_adjusted == TRUE){
    n_id <- length(unique(id_variable))
  }
  
  # Items to return
  N                <- NULL # Number of individual subjects
  Nstar            <- NULL # total number of scans across all participants, equal to N for cross-sectional case
  Jmax             <- NULL # maximum number of visits per participant
  Y_mixed_QxVNstar <- NULL # All data, is Qx(V * Nstar)
  Ji               <- NULL # number of visits per subject
  
  if (longit_adjusted == TRUE){
    
    if (is.list(Y)){

      if (debug == TRUE){print("Longitudinal adjusted list case")}
      
      # Can either be a list of list or a list of arrays
      if (is.list(Y[[1]])){
        
        if (debug == TRUE){print("Longitudinal adjusted list of lists case")}
        
        dimensions <- dim(Y[[1]][[1]])
        V          <- max(dimensions)
        Q          <- min(dimensions)
        
        # Case where each element of Y list is another list with Ji elements
        
        # Determine subjects, visits within subjects
        N  <- length(Y)
        if (N != n_id){
          stop(paste0("Time series data setup for ", N, " participants, but ID variable in data has ", n_id))
        }
        Ji <- rep(0, N)
        for (i in 1:N){
          Ji[i] <- length(Y[[i]])
        }
        Jmax  <- max(Ji)
        Nstar <- sum(Ji)
        
        # Store time series data
        Y_mixed_QxVNstar <- matrix(0.0, nrow = Q, ncol = Nstar * V)
        ending_index <- 0
        for (i in 1:N){
          starting_index <- ending_index + 1
          ending_index   <- ending_index + Ji[i]*V
          inds <- starting_index:ending_index
          # Format all visits for this subject
          Y_mixed_QxVJi <- format_list_time_courses_as_matrix(Y[[i]])
          # Add to the list of all data
          Y_mixed_QxVNstar[, inds] <- Y_mixed_QxVJi
        }
        
        
      } else {
        
        # Case where each element of the Y list is a row x col x Ji array
        if (debug == TRUE){print("Longitudinal adjusted list of arrays case")}
        
        dimensions <- dim(Y[[1]])
        V          <- max(dimensions[1:2])
        Q          <- min(dimensions[1:2])
        
        # Determine subjects, visits within subjects
        N  <- length(Y)
        if (N != n_id){
          stop(paste0("Time series data setup for ", N, " participants, but ID variable in data has ", n_id))
        }
        Ji <- rep(0, N)
        for (i in 1:N){
          Ji[i] <- dim(Y[[i]])[3]
        }
        Jmax  <- max(Ji)
        Nstar <- sum(Ji)
        
        # Store time series data
        Y_mixed_QxVNstar <- matrix(0.0, nrow = Q, ncol = Nstar * V)
        ending_index <- 0
        for (i in 1:N){
          starting_index <- ending_index + 1
          ending_index   <- ending_index + Ji[i]*V
          inds <- starting_index:ending_index
          # Format all visits for this subject
          Y_mixed_QxVJi <- format_array_time_courses_as_matrix(Y[[i]])
          # Add to the list of all data
          Y_mixed_QxVNstar[, inds] <- Y_mixed_QxVJi
        }
        
      }
      
    } else {
      # Can either be Q x V x Nstar or Q x V x N x Jmax.
      # note that Q and V might swap places here
      array_dimension <- dim(Y)
      if (length(array_dimension) == 3){
        
        dimensions <- dim(Y)
        V          <- max(dimensions[1:2])
        Q          <- min(dimensions[1:2])
        
        # For this case, we need to verify the ID information to get Ji
        N  <- n_id
        if (dimensions[3] != length(id_variable)){
          stop(paste0("Time series data setup for ", dimensions[3], " participants, but ID variable in data has ", length(id_variable)))
        }
        Ji <- rep(0, N)
        for (i in 1:N){
          Ji[i] <- sum(id_variable == unique(id_variable)[i])
        }
        Jmax  <- max(Ji)
        Nstar <- sum(Ji)
        
        # Store time series data
        Y_mixed_QxVNstar <- format_array_time_courses_as_matrix(Y)

      } else {
        
        dimensions <- dim(Y)
        V          <- max(dimensions[1:2])
        Q          <- min(dimensions[1:2])
        N          <- dimensions[3]
        Jmax       <- dimensions[4]
        if (N != n_id){
          stop(paste0("Time series data setup for ", N, " participants, but ID variable in data has ", n_id))
        }
        Ji <- rep(0, N)
        for (i in 1:N){
          Ji[i] <- sum(id_variable == unique(id_variable)[i])
        }
        Jmax  <- max(Ji)
        Nstar <- sum(Ji)
        
        # Store time series data
        Y_mixed_QxVNstar <- matrix(0.0, nrow = Q, ncol = Nstar * V)
        ending_index <- 0
        for (i in 1:N){
          starting_index <- ending_index + 1
          ending_index   <- ending_index + Ji[i]*V
          inds <- starting_index:ending_index
          # Format all visits for this subject
          Y_mixed_QxVJi <- format_array_time_courses_as_matrix(Y[,,i,1:Ji[i]])
          # Add to the list of all data
          Y_mixed_QxVNstar[, inds] <- Y_mixed_QxVJi
        }
        
      }
    }
    # BELOW is case where not longitudinal adjusted
  } else {
    Jmax <- 1
    if (is.list(Y)){
      N            <- length(Y)
      Y_mixed_QxVNstar <- format_list_time_courses_as_matrix(Y)
    } else {
      N            <- dim(Y)[3]
      Y_mixed_QxVNstar <- format_array_time_courses_as_matrix(Y)
    }
    Nstar <- N
    Ji    <- rep(1, N)
  }
  
  
  return(list(
    N                = N,
    Nstar            = Nstar,
    Jmax             = Jmax,
    Y_mixed_QxVNstar = Y_mixed_QxVNstar,
    Ji               = Ji
  ))
  
}





