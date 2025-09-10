#' Generate the model matrix for a BBSSModel based on a new set of covariate values
#' 
#' @param obj An object of class \code{"BBSSModel"} output from \link{bbss}.
#' @param newdata A \code{data.frame} (or tibble) with the covariate values at which
#'   to construct the model matrix.
#' 
#' @return A matrix containing the model matrix based on \code{newdata}.
#' 
#' @export
generate_model_matrix.BBSSModel <- function(obj, newdata){
  
  formula_input      <- obj$settings$formula
  fixed_only_formula <- lme4::nobars(formula_input)
  coding_information <- obj$settings$model_matrix_codings
  contrast_info      <- obj$settings$contrast_settings
  
  for (var in names(contrast_info)) {
    if (var %in% names(newdata)) {
      factor_levels <- rownames(attr(coding_information, "contrasts")[[var]])
      newdata[[var]] <- factor(newdata[[var]],
                               levels = factor_levels)
    }
  }
  
  model_matrix_new <- model.matrix(fixed_only_formula, data = newdata, contrasts.arg = contrast_info)
  model_matrix_new <- matrix(model_matrix_new, nrow = nrow(model_matrix_new))
  return(model_matrix_new)
}

#' Linear contrasts of covariate effects for a \code{BBSSModel}
#'
#' @details
#' Computes linear contrasts of the fixed-effect terms in a fitted \code{BBSSModel}
#' and returns location- and component-specific estimates (optionally with credible bands and
#' Savage–Dickey ratios). Supports one or many contrasts at once. NOTE - 
#' at this time the Bayes factor construction is under development. Please only 
#' make use of the credible bands. Also note that these credible bands are based
#' on the variational approximation and thus likely underestimate the actual
#' uncertainty.
#'
#' @param obj A fitted object of class \code{"BBSSModel"}.
#' @param contrast A numeric matrix specifying the linear
#'   contrast(s) over the intercept and fixed-effect coefficients. The required length
#'   (number of columns) is \eqn{1 + P}, where \eqn{P} is the number of fixed-effect
#'   terms stored in \code{obj$settings$P}. The first entry corresponds to the
#'   \emph{intercept} (population mean spatial map), and positions 2 through
#'   \eqn{P+1} correspond to the fixed-effect coefficients in the same order used at fit time.
#'   If a vector is supplied, it is treated as a single-row matrix. This argument
#'   can be directly obtained using \link{generate_model_matrix.BBSSModel}.
#' @param generate_credible_bands Logical; if \code{TRUE}, report lower/upper credible
#'   limits for each estimate.
#' @param credible_band_width Width of the credible band as a probability in \eqn{(0,1)};
#'   default \code{0.95} yields 95 percent credible bands.
#' @param component Optional integer vector of component indices to evaluate. Defaults
#'   to \code{1:obj$settings$Q}.
#' @param location Optional integer vector of voxel/vertex indices to evaluate. Defaults
#'   to \code{1:obj$settings$V}.
#' @param silent Logical; if \code{FALSE}, print progress messages.
#' @param contrast_labels Optional character vector of labels (length equal to number of
#'   contrasts) to include in the output.
#' @param calculate_savage_dickey Logical; if \code{TRUE}, compute the conditional
#'   Savage–Dickey Ratio for each contrast/location/component.
#' @param calculate_imp_reweights_savage_dickey Logical; if \code{TRUE}, compute an
#'   importance-sampling–calibrated Savage–Dickey estimate and its effective sample size
#'   (ESS). Enabling this implicitly sets \code{calculate_savage_dickey = TRUE}.
#' @param n_importance_resample Positive integer; number of importance resamples used
#'   when \code{calculate_imp_reweights_savage_dickey = TRUE}.
#'
#' @export
contrast.BBSSModel <- function(obj,
                          contrast,
                          generate_credible_bands = FALSE,
                          credible_band_width = 0.95,
                          component = NULL,
                          location  = NULL,
                          silent = TRUE,
                          contrast_labels = NULL,
                          calculate_savage_dickey = TRUE,
                          calculate_imp_reweights_savage_dickey = FALSE,
                          n_importance_resample = 1000){
  
  # Check arguments - if importance sampling requested then savage dickey ratios
  # should be returned as well
  if (calculate_imp_reweights_savage_dickey == TRUE){
    calculate_savage_dickey = TRUE
  }
  
  # Check if user has input invalid value for number of importance resamples
  if (n_importance_resample <= 0){
    stop("Number of importance samples must be > 0. Suggested value is at least 1000.")
  }
  
  if (!silent){cat("Checking contrasts are valid\n")}
  
  # Check if single contrast (vector) or multiple contrasts (matrix)
  if (is.vector(contrast)){
    contrast = matrix(contrast, nrow = 1)
  }
  
  # Store the number of contrasts as a variable for later use
  n_contrast <- nrow(contrast)
  
  # Verify the contrast dimensions work
  n_columns_of_contrast <- ncol(contrast)
  if (n_columns_of_contrast != (1 + obj$settings$P)){
    stop(paste0("Specified contrast has length of ", n_columns_of_contrast, 
                ", but expecting ", (1 + obj$settings$P)))
  }
  
  # Setup the frame we want to look at:
  ## First column is component. Second column is location. Third column is contrast.
  component_settings <- 1:obj$settings$Q
  if(!is.null(component)) {component_settings = component}
  location_settings  <- 1:obj$settings$V
  if(!is.null(location)) {location_settings = location}
  contrast_settings  <- 1:n_contrast
  # Want column order to be component, voxel, contrast, but want them
  # to alternate fastest in opposite order. This is what we build it backwards
  # and then reverse column ordering in constructing contrast grid below
  contrast_grid      <- expand.grid(contrast_settings,
                                    as.double(location_settings),
                                    as.double(component_settings))[,c(3,2,1)]
  contrast_grid      <- as.matrix(contrast_grid)
  
  # Evaluate
  if (!silent){cat("Estimating contrasts\n")}
  
  # Choice of fn depends on what terms are included in the model
  contrast_result <- NULL
  if (obj$settings$covariate_adjusted == TRUE){
    if (obj$settings$site_adjusted == TRUE){
      if (obj$settings$longit_adjusted == TRUE){
        stop("TODO cov site longit")
      } else {
        contrast_result <- eval_LC_cross_sectional_covadj_siteadj(contrast,
                                                                  contrast_grid,
                                                                  obj$misc_terms$Y_unmixed, obj$misc_terms$X_with_int,
                                                                  obj$misc_terms$site_indicators,
                                                                  obj$population_map_results$cluster_probabilities,
                                                                  obj$population_map_results$cluster_precisions,
                                                                  obj$population_map_results$cluster_means,
                                                                  obj$covariate_effect_results$local_hs_precision,
                                                                  obj$error_precision,
                                                                  obj$covariate_effect_results$global_hs_precision,
                                                                  as.double(obj$site_effect_results$site_effect_precision),
                                                                  obj$settings$V,
                                                                  generate_credible_bands,
                                                                  1.0 - credible_band_width,
                                                                  calculate_savage_dickey,
                                                                  calculate_imp_reweights_savage_dickey,
                                                                  n_importance_resample)
      }
      
    } else {
      if (obj$settings$longit_adjusted == TRUE){
        contrast_result <- eval_LC_longitudinal_covadj(contrast,
                                                       contrast_grid,
                                                       obj$misc_terms$Y_unmixed,
                                                       obj$misc_terms$X_with_int,
                                                       obj$misc_terms$Xrand,
                                                       obj$population_map_results$cluster_probabilities,
                                                       obj$longitudinal_results$RE_cluster_probs,
                                                       obj$population_map_results$cluster_precisions,
                                                       obj$population_map_results$cluster_means,
                                                       obj$covariate_effect_results$local_hs_precision,
                                                       obj$longitudinal_results$RE_precision_clusters,
                                                       obj$error_precision,
                                                       obj$covariate_effect_results$global_hs_precision,
                                                       obj$settings$V,
                                                       generate_credible_bands,
                                                       1.0 - credible_band_width,
                                                       calculate_savage_dickey,
                                                       calculate_imp_reweights_savage_dickey,
                                                       n_importance_resample)
      } else {
        contrast_result <- eval_LC_cross_sectional_covadj(contrast,
                                                          contrast_grid,
                                                          obj$misc_terms$Y_unmixed, obj$misc_terms$X_with_int,
                                                          obj$population_map_results$cluster_probabilities,
                                                          obj$population_map_results$cluster_precisions,
                                                          obj$population_map_results$cluster_means,
                                                          obj$covariate_effect_results$local_hs_precision,
                                                          obj$error_precision,
                                                          obj$covariate_effect_results$global_hs_precision,
                                                          obj$settings$V,
                                                          generate_credible_bands,
                                                          1.0 - credible_band_width,
                                                          calculate_savage_dickey,
                                                          calculate_imp_reweights_savage_dickey,
                                                          n_importance_resample)
      }
    }
    # Everything below this is NOT covariate adjusted
  } else {
    if (obj$settings$site_adjusted == TRUE){
      print("TODO")
    } else {
      print("TODO")
    }
  }
  
  # Format as tibble
  if (!silent){cat("Formatting results\n")}
  
  # (1) Format the grid
  contrast_grid           <- as_tibble(contrast_grid)
  colnames(contrast_grid) <- c("Component", "Location", "Contrast")
  if (!is.null(contrast_labels)){
    contrast_grid <- contrast_grid %>%
      mutate(Contrast_Label = contrast_labels[Contrast])
  }
  
  # (2) Format the estimates
  colname_vector <- "Estimate"
  if (generate_credible_bands == TRUE){
    colname_vector <- c(colname_vector,
                        paste0(credible_band_width*100, "% CredI Lower"),
                        paste0(credible_band_width*100, "% CredI Upper"))
  }
  if (calculate_savage_dickey == TRUE){
    colname_vector <- c(colname_vector, "Cond. Savage-Dickey")
  }
  if (calculate_imp_reweights_savage_dickey == TRUE){
    colname_vector <- c(colname_vector, "IS Calibrated Savage-Dickey", "IS ESS")
  }
  colnames(contrast_result) <- colname_vector
  
  # (3) Combine and return
  final_result <- bind_cols(contrast_grid, contrast_result)
  
  return(final_result)
}