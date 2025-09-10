#' Compare two sets of spatial maps
#' 
#' @param map1 A spatial map of dimension \eqn{Q \times V}.
#' @param map2 A spatial map of dimension \eqn{Q \times V} or \eqn{V \times Q}.
#' @param Q Optional argument specifying the number of components. Only needed if \code{map1} or \code{map2} are not already \eqn{Q \times V}.
#' 
#' @return A list object with the following fields:
#' 
#' \itemize{
#'   \item \code{correlation_unsorted}: Pearson correlation between the two maps using their raw ordering.
#'   \item \code{sorted_order}: Vector giving the re-indexing of \code{map1} that makes it most closely match \code{map2}.
#'   \item \code{correlation_sorted}: Pearson correlation between the two maps after re-indexing.
#' }
#' 
#' 
#' @export
compare_2d_map_correlations <- function(map1, map2, Q = NULL, cormethod = "pearson"){
  
  if (is.null(Q)){
    Q = nrow(map1)
  }
  
  if (dim(map1)[2] != Q){
    map1 <- t(map1)
  }
  
  if (dim(map2)[2] != Q){
    map2 <- t(map2)
  }
  
  cor_mat <- cor(map1, map2, method = cormethod)
  
  # Get the sorted correlations
  abs_cor <- abs(cor_mat)
  sorted_order <- rep(0, Q)
  for (i in 1:Q){
    maxcor = which(abs_cor == max(abs_cor), arr.ind = T)[1,]
    sorted_order[maxcor[2]] =  maxcor[1]
    abs_cor[,maxcor[2]] = 0
    abs_cor[maxcor[1],] = 0
  }
  
  cor_mat_sorted <- cor(map1[,sorted_order], map2, method = cormethod)
  
  return(list(
    correlation_unsorted = cor_mat,
    sorted_order = sorted_order,
    correlation_sorted = cor_mat_sorted
  ))
  
}

#' Reorder the results of a blind source separation to match a set of reference images
#' 
#' @param obj an object of \link{class} "BBSSModel"
#' @param reference_img A matrix of dimension \eqn{Q \times V} or \eqn{V \times Q} containing a set of spatial maps to use as a reference for re-ordering the contents of \code{obj}.
#' 
#' @return An object of \link{class} "BBSSModel" containing the results stored in \code{obj} re-indexed so that the components are in the order closest matching the components in \code{reference_img}.
#' 
#' @details
#' The order of the components in the BSS decomposition is not identifiable, but it is often the case that we want to compare multiple sets of components. Doing so requires re-indexing the components to match as closely as possible. This function calculates the correlation between the components stored in \code{obj} and the maps provided in \code{reference_img}. The any component-specific fields in \code{obj} are then re-ordered to obtain the closest possible matching.
#' 
#' @export
reorder_to_reference <- function(obj, reference_img){
  
  # Extract the dimensions from the obj
  Q <- obj$settings$Q
  V <- obj$settings$V
  P <- obj$settings$P
  
  # Check if reference image should be transposed
  rows_ref <- nrow(reference_img)
  cols_ref <- ncol(reference_img)
  if (rows_ref %in% c(Q, V) == FALSE){
    stop(paste0("Number of rows in reference image, ", rows_ref,
                " does not match either dimension of analysis result (V=",
                V, ", Q=", Q, ")"))
  }
  if (cols_ref %in% c(Q, V) == FALSE){
    stop(paste0("Number of columns in reference image, ", cols_ref,
                " does not match either dimension of analysis result (V=",
                V, ", Q=", Q, ")"))
  }
  if (rows_ref != V){reference_img = t(reference_img)}
  
  # Get the correlation between the population-level estimates from BBSS object
  # and the reference image
  alignment_result <- compare_2d_map_correlations(obj$population_map_results$population_spatial_maps,
                                                  reference_img, Q = Q)
  
  # Get the required sign flips and reordering
  alignment_order <- alignment_result$sorted_order
  # Currently not flipping sign since it influences the underlying DPM values
  alignment_sflip <- diag(sign(diag(alignment_result$correlation_sorted)))

  # Go through and reorder/flip all relevant attributes of the object
  new_obj <- obj
  
  ## Population Map Results ----
  new_obj$population_map_results$population_spatial_maps <- 
    obj$population_map_results$population_spatial_maps[, alignment_order]
  
  for (iq in 1:Q){
    vrange   <- ((iq-1)*V + 1):(iq * V)
    newrange <- ((alignment_order[iq]-1)*V + 1):(alignment_order[iq] * V)
    new_obj$population_map_results$cluster_probabilities[, vrange] <- 
      obj$population_map_results$cluster_probabilities[, newrange]
  }
  
  ## Covariate Effect Results ----
  for (ip in 1:P){
    new_obj$covariate_effect_results$beta_hat[,,ip] <- 
      obj$covariate_effect_results$beta_hat[,alignment_order,ip]
  }
  for (iq in 1:Q){
    vrange   <- ((iq-1)*V + 1):(iq * V)
    newrange <- ((alignment_order[iq]-1)*V + 1):(alignment_order[iq] * V)
    new_obj$covariate_effect_results$local_hs_precision[, vrange] <- 
      obj$covariate_effect_results$local_hs_precision[, newrange]
  }
  
  ## Site Effect Results ----
  if (obj$settings$site_adjusted == TRUE){
    for (is in 1:dim(obj$site_effect_results$site_effects_hat)[3]){
      new_obj$site_effect_results$site_effects_hat[,,is] <- 
        obj$site_effect_results$site_effects_hat[,alignment_order,is]
    }
  }
  
  ## Subject-Specific Maps ----
  if (!is.null(obj$subject_map_results)){
    
    for (i in 1:dim(obj$subject_map_results)[3]){
      for (j in 1:dim(obj$subject_map_results)[4]){
        new_obj$subject_map_results[,,i,j] <- obj$subject_map_results[, alignment_order, i, j] 
      }
    }
    
    
  }
  
  ## Misc Terms ----
  warning("TODO double check that unmixed Y is in correct order in reorder function (check impl)")
  for (iq in 1:Q){
    vrange   <- ((iq-1)*V + 1):(iq * V)
    newrange <- ((alignment_order[iq]-1)*V + 1):(alignment_order[iq] * V)
    new_obj$misc_terms$Y_unmixed[, vrange] <- 
      obj$misc_terms$Y_unmixed[, newrange]
  }
  
  warning("TODO - still need to realign Ai (not yet implemented)")

  return(new_obj)
}