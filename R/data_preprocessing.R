#' Prewhiten fMRI time courses
#' 
#' @param data Numeric matrix of fMRI time courses, with \eqn{T} timepoints (rows) and \eqn{V} spatial locations (columns).
#' @param Q Integer. Target dimensionality after prewhitening (must match the number of components to be estimated in blind source separation). 
#' 
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{Y}: a \eqn{V \times Q} matrix of prewhitened time courses.
#'   \item \code{eigenvalues}: numeric vector of eigenvalues of the temporal covariance matrix.
#' }
#' 
#' @examples
#' 
#' # Simulate example time courses (T=50, V=100)
#' set.seed(100)
#' simulated_time_course <- matrix(rnorm(100*50), nrow = 50, ncol = 100)
#'
#' # Prewhiten to Q = 5 components
#' prewhitened_time_courses <- prewhiten_participant(simulated_time_course, Q = 5)
#' cov(prewhitened_time_courses$Y)
#' 
#' @export
prewhiten_participant <- function(data, Q){
  
  # Check for missing values (likely bad mask step)
  if (any(is.na(data))){
    stop("`data` contains NA values.")
  }
  
  # Check for Inf values 
  if (any(is.infinite(data))){
    stop("`data` contains Inf values.")
  }
  
  # Dimensions, assumes data is Ti x V
  Ti = nrow(data)
  V  = ncol(data)
  
  # Verify input
  if (Ti < Q){
    stop(paste0("Dimension reduction to ", Q, " points requested but subject's",
                " data only contains ", Ti, " time points."))
  }
  if (V < Ti){
    warning("Number of spatial locations is less than number of time points Check if data matrix",
            " needs to be transposed.")
  }
  
  # Center each timepoint
  data <- sweep(data, 1, rowMeans(data))
  
  # PCA dimension reduction
  data_pca <- prcomp(t(data), center = FALSE)

  # Get the eigenvals
  eigvals_Q              <- data_pca$sdev[1:Q]^2

  # Eigenvectors
  U_eigen <- data_pca$rotation[, 1:Q]
  
  # Construct the transformation matrix
  tsfm_matrix <- t(U_eigen)
  tsfm_matrix <- diag( (eigvals_Q)^(-1/2) ) %*% tsfm_matrix
  
  # Construct the whitened data and flip it to be V x Q
  Y = t(tsfm_matrix %*% data)
  
  prewhitening_results <- list(Y = Y, eigenvalues = data_pca$sdev^2)
  
  # Returns whitened data
  return(prewhitening_results)
}


#' Perform PCA dimension reduction on fMRI time courses
#' 
#' @param data Numeric matrix of fMRI time courses with \eqn{T} rows (timepoints) 
#'   and \eqn{V} columns (spatial locations / voxels).
#' @param nPC Integer. Number of principal components to retain (\eqn{1 \le nPC \le T}). 
#'  Suggest to use around \eqn{2 \times } the number of components you intend to use in the blind source separation.
#' 
#' @return A numeric matrix of dimension \eqn{nPC \times V}. Each row is a spatial 
#'   principal component (a vector of loadings across voxels), scaled by the 
#'   corresponding singular value. Rows are returned (rather than columns) for 
#'   compatibility with temporal-concatenation group ICA workflows.
#'   
#' @examples
#' set.seed(100)
#' simulated_time_courses <- matrix(rnorm(100*20), nrow = 100, ncol = 20)  # 100 timepoints, 20 voxels
#' Y <- PCA_dimension_reduction(simulated_time_courses, nPC = 5)
#' dim(Y)  # 5 x 20 (5 components, each length 20)
#' 
#' @export
PCA_dimension_reduction <- function(data, nPC){
  
  # Dimensions
  Ti = nrow(data)
  V  = ncol(data)
  
  # Verify input
  if (Ti < nPC){
    stop(paste0("Dimension reduction to ", nPC, " points requested but subject's",
                " data only contains ", Ti, " time points."))
  }
  if (V < Ti){
    warning("Number of spatial locations is less than number of time points Check if data matrix",
            " needs to be transposed.")
  }
  
  # Center time points
  data <- sweep(data, 1, rowMeans(data))
  
  # PCA dimension reduction
  data_svd <- svd(t(data), nu = nPC)

  data_pc <- t(data_svd$u %*% diag(data_svd$d[1:nPC]))

  # Returns first nPC principal components
  return(data_pc)
  
}



#' Obtain an initial set of spatial components using FastICA
#' 
#' @param data Numeric matrix of dimension \eqn{Q \times V}, typically the output of 
#'   \link{PCA_dimension_reduction}, where \eqn{Q} is the number of components and 
#'   \eqn{V} is the number of spatial locations (voxels).
#' 
#' @return A numeric matrix of dimension \eqn{Q \times V}. Each row corresponds to an 
#'   estimated spatial component (loadings across voxels) obtained from the FastICA 
#'   algorithm.
#' 
#' @details
#' This function applies FastICA to the reduced data to provide an initial set of spatial 
#' components. These can be used as starting values for the population-average spatial maps 
#' in a subsequent \code{bbss} model fit.
#' 
#' @seealso 
#'    Related function: \link{PCA_dimension_reduction} for performing PCA prior to ICA.
#' 
#' @export
obtain_spatial_map_guess_fastica <- function(data){
  
  Q = nrow(data)
  V = ncol(data)
  
  # Apply fastica
  ica_decomposition <- ica::icafast(t(data), Q)
  
  S0 = ica_decomposition$S
  
  return(S0)
}



#' Create a mask
#' @param image A nifti image or 3d array to be used to generate the brain mask
#' @param threshold Everything greater than this value will be included in the mask
#' @export
create_mask <- function(image, threshold = 0.0){
  
  # Input validation
  if (inherits(image, "niftiImage")){
    
  } else {
    # Handle array case, note, nifti would also pass this, hence else statement
    if (is.array(image)){
      if (length(dim(image)) == 3){
        
      } else {
        error("If an array is provided as the image argument, it must be 3D")
      } # end of check that array is 3d
    } else {
      error("image argument must be a niftiImage or a 3d array")
    } # end of check that is array
  }
  
  mask_array         <- 1.0 * (image > threshold)
  mask_voxel_indices <- which(mask_array == 1)
  
  # Return list
  return(list(
    mask_array         = mask_array,
    mask_voxel_indices = mask_voxel_indices
  ))
  
}

#' Apply a mask to a Nifti file or array
#' @param image A 3- or 4-dimensional array to be masked
#' @param mask Mask to be applied to the image. Must be created using \link{create_mask}.
#' @param return_type either "flattened" for a matrix result or "volume" if the original volumetric structure should be preserved. Default is "flattened".
#' @export
apply_mask <- function(image,
                       mask,
                       return_type = c("flattened", "volume")){
  
  # Get the type of data to return
  return_type <- match.arg(return_type)
  
  # Determine if image is 3D or 4D
  dimension <- length(dim(image))
  if (!(dimension %in% c(3,4))){
    stop(paste0("Input image is ", dimension, " dimensional. Must be either 3d or 4d."))
  }
  
  # Verify that the mask and the image to be masked have the same size
  # (note - this does not guarantee same template space, but can rule out many problems)
  
  # Determine the number of volumes
  n_volume <- 1
  if (dimension == 4){
    n_volume <- dim(image)[4]
  }
  
  # Apply the mask
  masked_image <- NULL
  
  # Flat (matrix, T x V) version of the result
  if (return_type == "flattened"){
    # Allocate space for result
    masked_image <- matrix(0, nrow = n_volume, ncol = length(mask$mask_voxel_indices))
    
    if (dimension == 3){
      masked_image[1, ] <- image[mask$mask_voxel_indices]
    } else {
      # Loop over volumes and apply the mask
      for (i_vol in 1:n_volume){
        image_volume_i        <- image[,,,i_vol]
        masked_image[i_vol, ] <- image_volume_i[mask$mask_voxel_indices]
      }
    }
    
  } # end of flattened version
  
  # Array version of the result
  if (return_type == "volume"){
    masked_image <- image
    if (dimension == 3){
      masked_image[,,] <- image[,,] * mask$mask_array
    } else {
      for (i_vol in 1:n_volume){
        masked_image[,,,i_vol] <- image[,,,i_vol] * mask$mask_array
      }
    }
  }
  
  return(masked_image)
  
}
