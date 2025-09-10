# reference_nifti should be either a file OR a nifti object
export_as_nifti <- function(obj,
                            mask,
                            folder,
                            reference_nifti,
                            save_pattern = c("collected", "individual"),
                            include_S_maps = TRUE,
                            include_Beta_maps = TRUE,
                            include_center_effect_maps = TRUE,
                            include_subject_maps = FALSE){
  
  # Properties
  Q <- ncol(obj$population_map_results$population_spatial_maps)
  V <- nrow(obj$population_map_results$population_spatial_maps)
  P <- dim(obj$covariate_effect_results$beta_hat)[3]
  
  # Verify class of obj
  # if (!(class(obj) == "BBSSModel")){
  #   stop("obj argument must be a BBSSModel class object.")
  # }
  
  # Verify that the output folder exists
  if (!dir.exists(folder)){
    stop(paste0("Folder: ", folder, " does not exist."))
  }
  
  # Verify that the provided reference image is valid
  # if user provided a string, check if it is filepath to nifti file
  if (is.character(reference_nifti)){
    if (!file.exists(reference_nifti)){
      stop(paste0("File: ", reference_nifti, " does not exist."))
    }
    reference_nifti <- RNifti::readNifti(reference_nifti)
  } else {
    # Check if provided object is a nifti image
  }
  
  # Setup shell of nifti object that will be used (e.g. limits)
  nifti_color_limits_def <- quantile(abs(obj$population_map_results$population_spatial_maps), 0.99)
  nifti_dim              <- dim(reference_nifti)
  nifti_base_image       <- array(NA, dim = nifti_dim)
  nifti_hdr              <- RNifti::niftiHeader(reference_nifti)
  nifti_hdr$cal_min            <- -nifti_color_limits_def
  nifti_hdr$cal_max            <- nifti_color_limits_def
  
  # Make the base nifti image
  base_nifti_image <- RNifti::asNifti(nifti_base_image, reference = nifti_hdr)

  # Save population average maps
  for (iq in 1:Q){
    S0_q <- c(obj$population_map_results$population_spatial_maps[, iq])
    new_nii <- base_nifti_image
    new_nii[mask$mask_voxel_indices] <- S0_q
    new_nii_fname <- file.path(folder, paste0("population_average_network_", iq, ".nii"))
    RNifti::writeNifti(new_nii, new_nii_fname)
  }
  
  # Save covariate effect maps
  for (iq in 1:Q){
    for (ip in 1:P){
      beta_qp <- c(obj$covariate_effect_results$beta_hat[,iq,ip])
      new_nii <- base_nifti_image
      new_nii[mask$mask_voxel_indices] <- beta_qp
      new_nii_fname <- file.path(folder, paste0("covariate_", ip, "_effect_network_", iq, ".nii"))
      RNifti::writeNifti(new_nii, new_nii_fname)
    }
  }
  
  # Save center effect maps
  if (!is.null(obj$site_effect_results)){
    warning("Site effect nifti saving not yet implemented")
  }
  
  # Save subject-level maps
  if (!is.null(obj$subject_map_results)){
    warning("Subject effect nifti saving not yet implemented")
  }
  
}