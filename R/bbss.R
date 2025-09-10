#' Run Bayesian blind source separation on prewhitened fMRI time courses
#' 
#' @param Y The prewhitened fMRI time courses. Each time course should be prewhitened to \eqn{Q} timepoints, where \eqn{Q} is the number of spatial source signals (components) to be estimated. See details below.
#' @param initial_guess Starting values for some of the model parameters. At a minimum, this must include an initial value for the population average brain network map in a field called \code{init_S} that is of dimension \eqn{Q \times V}, where \eqn{Q} is the number of maps to estimate and \eqn{V} is the number of spatial locations. See details below.
#' @param formula A \link{formula} object specifying the covariate effects and random terms (in the case of longitudinal data). Follows the same structure as the formula from a \link{lm} or \link{lmer} function. See details below.
#' @param data A \link{tibble} or \link{data.frame} containing the covariate effects. Anything in the \code{formula} argument must also have a corresponding column in data.
#' @param maxit Number of iterations to run when optimizing the parameters of the variational distribution. Note that this currently does not involve early stopping (DEV VERSION).
#' @param Hmax Maximum number of clusters in the DPM for the spatial source signals
#' @param Hmax_Rqv  Maximum number of clusters in the DPM for the covariance of the random effects
#' @param concentration Concentration parameter for the DPM for the spatial source signals
#' @param concentration_Rqv Concentration parameter for the DPM for the covariance of the random effects
#' @param site_variable Name of the column of \code{data} containing the site variable indicators. If this is included, a site-adjusted BBSS model will be fit that incorporates random intercepts for site.
#' @param coding_scheme Specify the coding scheme to apply to any categorical variables in the data. Options are "effects" coding (DEFAULT) or "reference" cell coding (IN DEV)
#' @param base_measure_lambda Hyperparameter for lambda in the base measure for the covariance of the random effects
#' @param base_measure_mu Hyperparameter for the mean in the base measure for the covariance of the random effects. Default is 0.
#' @param base_measure_alpha Hyperparameter for the shape in the base measure for the covariance of the random effects
#' @param base_measure_beta Hyperparameter for scale in the base measure for the covariance of the random effects
#' @param print_freq How often summary information about the algorithm should be printed to the console. Default is every 5 iterations.
#' @param store_subject_maps Boolean variable storing whether the individual subject-level map estimates should be stored. Default is FALSE, since these can use up a substantial amount of memory.
#' 
#' @details
#' The function fits a variational approximation to the SparseBayes ICA (Lukemire, 2023)  and REMBRAiNDT models (under review). Models for \code{bbss} are specified symbolically, similar to \link{lm} and \link{lmer}.
#' 
#' @section Preparing the \code{Y} Argument:
#' \code{bbss} expects the \code{Y} argument to contain the prewhitened time courses for each participant and scan (in the event of multiple scans per subject). The function will accept either a list or array argument for \code{Y} in the following forms:
#' 
#' \itemize{
#'   \item A 3D Array that is \eqn{Q \times V \times N^*}, where \eqn{N^*} is the total number of scans across all subjects. This is equivalent to concatenating the data along the 3rd dimension. Note that for datasets without repeated measures this is just \eqn{Q \times V \times N}.
#'   \item A 4D Array that is \eqn{Q \times V \times N \times J^*}, where \eqn{J^*} is the maximum number of scans any subject has. Note that this can be memory inefficient if some subjects have many scans and others only have a few.
#'   \item A list of length \eqn{N}, where list element \eqn{i} is another list containing \eqn{J_i} matrices of dimension \eqn{Q \times V}. 
#'   \item A list of length \eqn{N}, where list element \eqn{i} is an array of size \eqn{Q \times V \times J_i}.
#' }
#' 
#' The prewhitening can be carried out using the \link{prewhiten_participant} function provided by this library. The most convenient format will depend on how your data is stored and whether or not you have multiple scans per subject. Finally, the internal functions that parse the input will attempt to fix some orientation errors, for example if a \eqn{V \times Q \times N} array is provided, the library will handle permuting the indices as needed.
#' 
#' @section Preparing Initial Guess:
#' An initial guess for the model parameters must be provided in the \code{initial_guess} argument. The \code{initial_guess} argument should be structured as a list with the following elements:
#' 
#' \itemize{
#'   \item \code{init_S}: (Required) A matrix (\eqn{Q \times V} or \eqn{V \times Q}) containing the initial values for the spatial maps (components). This can be obtained using FASTICA via the \link{obtain_spatial_map_guess_fastica} function provided in this library. \bold{This is the only required element of the \code{initial_guess} list}. 
#'   \item \code{init_Ai}: (Optional) An array (\eqn{Q \times Q \times N^*}) containing the initial values for each subject's mixing matrix. Note that if this is not provided, \link{bbss} will calculate a guess for you based on a regression of the prewhitened time courses onto \code{init_S}.
#' }
#' 
#' The remaining parameters will be initialized within the \link{bbss} function. 
#' 
#' @section Preparing the \code{data} Argument:
#' The \code{data} argument should be either a \link{tibble} or \link{data.frame} containing the values of any clinical or demographic covariates. This should be prepared in the same way as the \code{data} argument from \link{lm} or \link{lmer}. It can also be useful to include a column containing the filepath to each nifti file, as this can be used to help ensure that the order of \code{Y} is the same as the order of the covariates by applying the preprocessing one-row-at-a-time using this column.
#' 
#' 
#' @section Specifying the model formula:
#' The model \code{formula} should be provided with no term on the left-hand side. Some example forms:
#' \itemize{
#'   \item \code{~ sex + age}: Will fit a model with two estimated beta maps: one for sex and one for age. Requires \code{data} to have columns named "sex" and "age".
#'   \item \code{~ sex*age}: Will fit a model with three estimated beta maps: one for sex, one for age, and one for a sex by age interaction. Requires \code{data} to have columns named "sex" and "age".
#'   \item \code{~ sex + age + (1 | ID)}: Will fit a model with two estimated beta maps: one for sex and one for age. Additionally, this model will include a random intercept corresponding to \code{ID}. Requires \code{data} to have columns named "sex", "age", and "ID".
#' }
#' 
#' Note that the column names can be any valid column name, and the above uses of "sex", "age", "ID" are just examples.
#' 
#' @return An object of class BBSSModel containing parameter estimates, fitted spatial maps, covariate effects, and algorithm diagnostic information.
#' 
#' @references
#' Lukemire, J., Pagnoni, G., & Guo, Y. (2023). \emph{Sparse Bayesian modeling of hierarchical independent component analysis: Reliable estimation of individual differences in brain networks}.
#' Biometrics, 79(4), 3599-3611.
#' 
#' @seealso \link{prewhiten_participant}, \link{obtain_spatial_map_guess_fastica}
#' 
#' @export
bbss <- function(Y, initial_guess,
                 formula = NULL, data = NULL,
                 maxit = 100, Hmax = 10, Hmax_Rqv = 5,
                 concentration = 1.0, concentration_Rqv = 1.0,
                 site_variable = NULL,
                 coding_scheme = "effects",
                 base_measure_lambda = 1.0, base_measure_mu = 0,
                 base_measure_alpha = 2.0, base_measure_beta = 1.0,
                 print_freq = 5,
                 store_subject_maps = FALSE){
  
  # Y should be a list, where each element is a Q x V matrix
  # initial_guess must contain a field called S0 that is Q x V
  
  # Tracking which type of model we are fitting
  covariate_adjusted <- FALSE
  site_adjusted      <- ifelse(is.null(site_variable), FALSE, TRUE)
  longit_adjusted    <- FALSE
  
  # Evaluate user-specified formula
  # TODO check for incorrectly specified random term
  X_no_int <- NULL
  X        <- NULL
  Xrand    <- NULL
  P <- 0
  Prand <- 0
  model_matrix_codings <- NULL
  id_variable        <- NULL # name of the id variable in the covariates
  id_variable_values <- NULL # actual value for the id variable in the covariates
  if (!is.null(formula)){
    
    required_variables       <- all.vars(formula)

    # Check that variables exist in data object
    # if not, give detailed error
    # TODO - this does not handle interactions correctly yet!
    print("TODO Interaction checks")
    check_variables_in_data <- required_variables %in% names(data)
    if (any(check_variables_in_data == FALSE)){
      stop(paste0("The provided data object is missing columns for: ",
                  required_variables[check_variables_in_data==FALSE],
                  collapse = ", "))
    }
    
    # Check if the user has specified any random effects
    all_formula_names      <- all.names(formula)
    contains_random_effect <- "|" %in% all.names(formula)
    if (contains_random_effect == TRUE){
      
      # First check that too many terms were not input
      n_random_term <- sum(all.names(formula) == "|")
      if (n_random_term > 1){
        stop("Only one random effect (corresponding to subject) should be present in the equation.
             To specify a center clustering variable, use the site_variable argument.")
      }
      
      # Split up the random effect into 3 pieces
      # (1) the |, (2) the lhs terms (3) the id variable
      random_effect_specification <- lme4::findbars(formula)[[1]]
      
      # Second, check that the user specified AT least a random intercept
      # the main potential problem would be if the user input a 0,
      # we do not allow this
      print("TODO stop and correct user if they tried to remove random intercept")
      random_effect_terms <- random_effect_specification[[2]]
      
      # Third, check that any random slopes were included in the data object
      required_random_variables <- all.vars(random_effect_terms)
      check_random_variables_in_data <- required_random_variables %in% names(data)
      if (any(check_random_variables_in_data == FALSE)){
        stop(paste0("The provided data object is missing columns for: ",
                    required_random_variables[check_random_variables_in_data==FALSE],
                    collapse = ", "))
      }
      
      # Fourth, check that the id variable was included in the data object
      # additionally, make sure user did not set up a nested structure
      id_variable <- paste0(random_effect_specification[[3]])
      id_included_check <- id_variable %in% names(data)
      if (id_included_check == FALSE){
        stop(paste0(id_variable, " was listed as the subject id variable, but is not present in the provided data argument."))
      }
      id_variable_values = data %>% pull(id_variable)
      
      # At this point all checks have passed, we still need to:
      # 1) Set the longit_adjusted flag
      # 2) Construct the RE model matrix (done below after constructing regular model matrix first)
      
      longit_adjusted <- TRUE
     
    }

    # Use user-specified coding scheme for categorical variables
    ## First figure out which variables are character variables
    # nobars is used to remove any random terms for this part:
    fixed_terms_formula  <- lme4::nobars(formula)
    
    # Now get the coding
    factor_vars <- names(data)[sapply(data, is.factor)]
    char_vars   <- names(data)[sapply(data, is.character)]
    combined_cat_vars <- c(factor_vars, char_vars)
    combined_cat_vars <- combined_cat_vars[combined_cat_vars %in% all.vars(fixed_terms_formula)]
    contrast_settings <- setNames(rep(list(contr.sum), length(combined_cat_vars)), combined_cat_vars)
    
    # Construct model matrix from variables
    model_matrix_codings <- model.matrix(data = data,
                                             object = fixed_terms_formula,
                                             contrasts = contrast_settings)
    
    X_no_int <- as.matrix(model_matrix_codings)[,-1]
    
    # Set the flag that this model is covariate adjusted
    covariate_adjusted <- TRUE
    
    # Count the number of beta terms
    P <- ncol(X_no_int)
    
    # Create version of model matrix that includes intercept
    X <- cbind(1, X_no_int)
    
    # If longitudinal data was included, also setup the model matrix for
    # the random effects
    if (longit_adjusted == TRUE){
      combined_random_cat_vars <- combined_cat_vars[combined_cat_vars %in% required_random_variables]
      random_contrast_settings <- setNames(rep(list(contr.sum),
                                               length(combined_random_cat_vars)),
                                           combined_random_cat_vars)

      # Create a temporary version of the codings. Note that this does not account
      # for the fact that we need separate columns for each subject
      random_model_matrix_codings <- model.matrix(data = data,
                                           object = formula(paste0("~ ", deparse(random_effect_specification[[2]]))),
                                           contrasts = random_contrast_settings)
      
      # If no slopes were specified, the above will only have a single row. We can 
      # tweak this here
      if (nrow(random_model_matrix_codings) == 1){
        random_model_matrix_codings <- matrix(1.0, nrow = nrow(X_no_int))
      }
      
      # Now force this to be the correct dimensions using the id variable
      random_block_size <- ncol(random_model_matrix_codings)
      unique_subj       <- data %>% pull(!! id_variable) %>% unique
      N_unique_subj     <- length(unique_subj)
      new_model_matrix  <- matrix(0, nrow = nrow(random_model_matrix_codings),
                                  ncol = random_block_size * N_unique_subj)
      
      for (i_subj in 1:N_unique_subj){
        current_id   <- unique_subj[i_subj]
        row_index_id <- which(data %>% pull(!! id_variable) == current_id)
        col_index_id <- (random_block_size * (i_subj-1) + 1):(random_block_size*i_subj)
        new_model_matrix[row_index_id, col_index_id] <- random_model_matrix_codings[row_index_id, ]
      }
      
      Prand <- random_block_size
      
      Xrand <- new_model_matrix
    } # end of longitudinal
  }
  
  
  
  ## Time Series Formatting ----
  # Format the input data as a wide matrix
  # Here we want Q x NV or Q x NNstarV, since the main use will be unmixing this
  # object, which requires multiplication with QxQ A transpose
  formatting_time_series <- format_time_series_input(Y, longit_adjusted,
                                                     debug = FALSE,
                                                     id_variable = id_variable_values)
  N                <- formatting_time_series$N # Number of individual subjects
  Nstar            <- formatting_time_series$Nstar # total number of scans across all participants, equal to N for cross-sectional case
  Jmax             <- formatting_time_series$Jmax # maximum number of visits per participant
  Y_mixed_QxVNstar <- formatting_time_series$Y_mixed_QxVNstar # All data, is Qx(V * Nstar)
  Ji               <- formatting_time_series$Ji # number of visits per subject
  
  ## Initial Guess ----
  # Check if some initialization is required
  
  # First, check if a population average map has been provided. This is
  # required!
  if ( ("init_S" %in% names(initial_guess)) == FALSE ){
    stop(paste0("At a minimum, bbss requires the initial_guess argument to contain an init_S field."))
  }
  
  # Check if we need to initialize the mixing matrices, if so, we run a simple
  # dual regression procedure
  if ( ("init_Ai" %in% names(initial_guess)) == FALSE ){
    
    print("Obtaining Initial Mixing Matrices")
    
    initial_guess <- initialize_mixing_matrices(Y_mixed_QxVNstar,
                                                initial_guess$init_S,
                                                Nstar)
    
  }
  
  V <- nrow(initial_guess$init_S)
  Q <- ncol(initial_guess$init_S)
  
  if (V < Q){
    stop("An error appears to have occured, the number of locations is larger than the number of components")
  }
  
  # Allocate storage for the unmixed version of Y, which will be N x QV
  Y_unmixed <- matrix(0, nrow = Nstar, ncol = length(Y_mixed_QxVNstar) / Nstar)

  # Use the initial values of A to set Y_unmixed
  apply_unmixing_matrices(Y_unmixed, Y_mixed_QxVNstar, initial_guess$init_A)
  
  # Also determine the trace of YY' across all subjects - needed for later
  sum_trace_YYt <- sum(diag(Y_mixed_QxVNstar %*% t(Y_mixed_QxVNstar) ))
  
  # Allocate Variables ----
  
  ### Mixing Matrices ----
  E_Ai <- initial_guess$init_A + 0 
  
  ### Population Source Signals ----
  #### Expectations
  E_S           <- initial_guess$init_S + 0
  E_S_sq        <- E_S^2
  #### Free Variational Parameters TODO - do I store V in case of covariates?
  VariHyper_M_s_qv <- initial_guess$init_S + 0
  VariHyper_V_s_qv <- initial_guess$init_S + 0
  
  # Note - this is just a linear combination of other variables (S, Beta)
  # I'm just giving this a name for storage and reuse later.
  E_Si          <- matrix(0, nrow = Nstar, ncol = Q*V)
  trace_E_Si_sq <- matrix(0, nrow = Nstar); # only needed for covadj case
  
  ### Covariate Effects ----
  #### Terms below correspond to covariate effects and relevant horseshoe
  #### parameters
  #### Expectations
  E_Beta         = NULL
  E_Beta_sq      = NULL
  E_lambda_sq    = NULL
  E_tau_sq           = NULL
  E_tau_sq_inv       = NULL
  E_xi_sq            = NULL
  E_xi_sq_inv        = NULL
  if (covariate_adjusted == TRUE){
    E_Beta      <- matrix(0, nrow = P, ncol = Q * V)
    E_Beta_sq   <- matrix(0, nrow = P, ncol = Q * V)
    E_lambda_sq <- matrix(100.0, nrow = P, ncol = Q * V)
    E_tau_sq        <- 1.0
    E_tau_sq_inv    <- 1.0
    E_xi_sq         <- 1.0
    E_xi_sq_inv     <- 1.0
  }
  
  ### Free variational parameters
  ### Note - we don't store all of these, instead opting to store
  ### some items just as expectations, since some of the variance terms
  ### will be (P+1)x(P+1)xQxV
  
  ### S DPM Stick Breaking ----
  #### Expectations
  E_nu      <- rbeta(Hmax, 1, concentration)
  E_log_nu_h         <- rep(0, Hmax)
  E_log_1_minus_nu_h <- rep(0, Hmax)
  #### Free Variational Parameters
  VarHyper_A_nu_h    <- rep(0, Hmax)
  VarHyper_B_nu_h    <- rep(0, Hmax)
  
  ### S DPM Clusters ----
  S0_DPM_guess <- obtain_S0_DPM_guess(initial_guess$init_S, Hmax)
  
  #### Expectations
  E_mu_h           <- S0_DPM_guess$mu_h
  E_mu_h_sq        <- E_mu_h^2
  E_gamma_h_sq     <- S0_DPM_guess$sigma_h_sq
  E_log_gamma_h_sq <- log(E_gamma_h_sq)
  E_gamma_h_sq_inv <- 1 / E_gamma_h_sq
  #E_indicator_H <- array(0.01 / (Hmax - 1), dim = c(Hmax)) # pi(Hqv = h)
  pi_Zqv_eq_h <- matrix(0.01 / (Hmax-1), nrow = Hmax, ncol = V*Q)
  column_index <- 0
  for (iq in 1:Q){
    for (iv in 1:V){
      column_index <- column_index + 1
      h_qv <- S0_DPM_guess$Z_qv[iv, iq]
      pi_Zqv_eq_h[h_qv, column_index] <- 0.99
    }
  }
  #### Free Variational Parameters
  VarHyper_M_h    <- rep(0, Hmax)
  VarHyper_L_h    <- rep(0, Hmax)
  VarHyper_A_h    <- rep(0, Hmax)
  VarHyper_B_h    <- rep(0, Hmax)
  
  ### Longitudinal Parameters ----
  longit_varlist <- c("E_Rqv_nu", "E_Rqv_log_nu_h", "E_Rqv_log_1_minus_nu_h",
                      "VarHyper_A_Rqv_nu_h", "VarHyper_B_Rqv_nu_h",
                      "E_Rqv_inv", "E_log_det_Rqv_inv", "pi_Ztildeqv_eq_h")
    for (variable in longit_varlist) {
      assign(variable, NULL, envir = environment())
    }
  
  if (longit_adjusted == TRUE){
    ### Random Eff DPM Stick Breaking ----
    Rqv_DPM_guess <- obtain_Rqv_DPM_guess(Y_mixed_QxVNstar, E_Ai,
                                          E_S, Xrand, Prand, Hmax_Rqv,
                                          Beta = E_Beta, X = X)
    #### Expectations
    E_Rqv_nu               <- rbeta(Hmax_Rqv, 1, concentration_Rqv)
    E_Rqv_log_nu_h         <- rep(0, Hmax_Rqv)
    E_Rqv_log_1_minus_nu_h <- rep(0, Hmax_Rqv)
    #### Free Variational Parameters
    VarHyper_A_Rqv_nu_h    <- rep(0, Hmax_Rqv)
    VarHyper_B_Rqv_nu_h    <- rep(0, Hmax_Rqv)
    
    ### Random Eff DPM Clusters ----
    #### Expectations
    base_measure_scale <- diag(10.0, nrow = Prand)
    base_measure_df    <- Prand + 5
    E_Rqv_inv          <- Rqv_DPM_guess$E_Rqv_inv + 0
    E_log_det_Rqv_inv <- rep(0, Hmax_Rqv)
    E_tr_BiBit_qv     <- matrix(0, nrow = Prand, ncol = Prand *Q*V) # big, but req for cluster membership
    pi_Ztildeqv_eq_h  <- matrix(0.01 / (Hmax_Rqv-1), nrow = Hmax_Rqv, ncol = V*Q)
    column_index <- 0
    for (iq in 1:Q){
      for (iv in 1:V){
        column_index <- column_index + 1
        h_qv <- Rqv_DPM_guess$Z_Rqv[iv, iq]
        pi_Ztildeqv_eq_h[h_qv, column_index] <- 0.99
      }
    }
    #### Free Variational Parameters
    VarHyper_Rqv_df    <- rep(0.0, Hmax_Rqv)
    VarHyper_Rqv_scale <- matrix(0.0, nrow = Prand, ncol = Prand * Hmax_Rqv)
  }

  ### Site Effects ----
  E_b_site            <- NULL
  sum_E_b_site_sq     <- NULL
  E_sigma_site_sq     <- NULL
  E_sigma_site_sq_inv <- NULL
  site_indicators     <- NULL
  n_site              <- NULL
  if (site_adjusted == TRUE){
    
    # Build the coded version of the site matrix (1 if belongs to site)
    site_column <- data %>% dplyr::select(!! site_variable)
    
    # TODO warning if not enough sites
    
    # Apply reference cell coding to the site variable to obtain
    # a set of indicators
    site_counters <- as.numeric(factor(site_column %>% pull(!! site_variable)))
    n_site        <- max(site_counters)
    site_indicators <- matrix(0, nrow = Nstar, ncol = max(site_counters))
    for (i in 1:N){
      site_indicators[i, site_counters[i]] <- 1
    }
    
    # Initialize required variables
    E_b_site            <- matrix(0, nrow = n_site, ncol = Q*V)
    sum_E_b_site_sq     <- matrix(0, nrow = 1, ncol = 1)
    E_sigma_site_sq     <- 1.0
    E_sigma_site_sq_inv <- 1.0
  }
  
  ### Error Variance ----
  #### Expectations
  E_sigma_sq <- 1.0
  E_sigma_sq_inv <- 1 / E_sigma_sq
  #### Free Variational Parameters
  VariHyper_A_sigma_sq <- Nstar*Q*V/2 + Q*V/2
  if (covariate_adjusted == TRUE){
    VariHyper_A_sigma_sq <- Nstar*Q*V/2 + Q*V/2 + P*Q*V/2
  }
  VariHyper_B_sigma_sq <- E_sigma_sq * (VariHyper_A_sigma_sq - 1)
  
  # Misc Tracking Setup ----
  
  ### Relative Changes ----
  relchange_sigma_sq_error <- rep(0, maxit)
  relchange_sigma_sq_site  <- rep(0, maxit)
  relchange_tau_sq         <- rep(0, maxit)
  relchange_mu_h           <- matrix(0, nrow = Hmax, ncol = maxit)
  relchange_sigma_h_sq     <- matrix(0, nrow = Hmax, ncol = maxit)
  
  ### Computation Time ----
  time_S <- 0.0
  time_S_cluster_probs <- 0.0
  time_S_sticks <- 0.0
  time_S_cluster_param <- 0.0
  time_mixing_matrix <- 0.0
  time_unmixing_Y_update <- 0.0
  time_error_variance <- 0.0
  time_horseshoe_hyper <- 0.0
  time_site_var <- 0.0
  
  # CAVI Loop ----
  iter <- 1
  for (iter in 1:maxit){
    
    # Print current status
    if (iter/print_freq == round(iter/print_freq)){
      cat("Iteration", iter, "\n")
      cat("Current error variance: ", E_sigma_sq, "\n")
      
      if (site_adjusted){
        cat("Current site effect variance: ", E_sigma_site_sq, "\n")
      }
      
      if (longit_adjusted == TRUE){
        cat("Current Top RE Precision Terms: ", "\n")
        print_index_set <- 1:Prand 
        for (hr in 1:min(3, Hmax_Rqv)){
          cat(E_Rqv_inv[, print_index_set], "\n")
          cat("Corresponding df in Variational Approx: ", VarHyper_Rqv_df[hr], "\n")
          print_index_set <- print_index_set + Prand
        }
      }
    }
    
    ## Primary Updates ----
    
    ## Update the S latent cluster memberships ----
    spatial_cluster_probabilities_start_time = Sys.time()
    calibrate_spatial_map_cluster_probabilities(pi_Zqv_eq_h,
                                                E_log_nu_h,
                                                E_log_1_minus_nu_h,
                                                E_log_gamma_h_sq,
                                                E_gamma_h_sq_inv,
                                                E_mu_h,
                                                E_mu_h_sq,
                                                E_S,
                                                E_S_sq,
                                                E_sigma_sq_inv)
    
    if (any(is.na( pi_Zqv_eq_h ))){
      stop("NA in spatial cluster membership probabilities")
    }
    
    # print("CHEATING")
    # pi_Zqv_eq_h[,] <- 0.0001 / (Hmax-1)
    # for (iq in 1:Q){
    #   index_set <- ((iq-1)*V+1):(iq*V)
    #   cm_q <- c(true_values$true_s_cluster_membership[iq, ])
    #   for (h in 1:Hmax){
    #     pi_Zqv_eq_h[h, index_set[cm_q == h]] <- 0.9999
    #   }
    # }
    # print("DONE CHEATING")
    
    ### Timing
    elapsed_time <- Sys.time() - spatial_cluster_probabilities_start_time
    time_S_cluster_probs <- time_S_cluster_probs + as.numeric(elapsed_time)
    
    ## Update Spatial Maps ----
    spatial_map_start_time = Sys.time()
    ## TODO eventually these depend on covariate or not, and repeated measures or not
    ## update depends on whether or not covariates are included
    if (covariate_adjusted == FALSE){
      if (longit_adjusted == TRUE){
        if (site_adjusted == TRUE){
          stop("TODO not cov, yes logit, yes site")
        } else {
          stop("TODO not cov, yes logit, no site")
        }
      }
      else {
        if (site_adjusted == TRUE){
          cross_sectional_siteadj_spatial_map_updates(
            E_S, E_S_sq, E_b_site, sum_E_b_site_sq,
            E_Si, trace_E_Si_sq, Y_unmixed, site_indicators,
            pi_Zqv_eq_h, E_gamma_h_sq_inv, E_mu_h, 
            E_sigma_sq_inv, E_sigma_site_sq_inv
          )
        } else {
          cross_sectional_spatial_map_updates(E_S, E_S_sq, E_Si,
                                              Y_unmixed, pi_Zqv_eq_h,
                                              E_gamma_h_sq_inv,
                                              E_mu_h, E_sigma_sq_inv)
        }
      }
    } else {
      if (longit_adjusted == TRUE){
        if (site_adjusted == TRUE){
          # st <- Sys.time()
          longitudinal_covadj_siteadj_spatial_map_updates(E_S, E_S_sq,
                                                          E_Beta, E_Beta_sq,
                                                          E_b_site, sum_E_b_site_sq,
                                                          E_Si, trace_E_Si_sq,
                                                          E_tr_BiBit_qv,
                                                          Y_unmixed,
                                                          X, Xrand, site_indicators,
                                                          pi_Zqv_eq_h, pi_Ztildeqv_eq_h,
                                                          E_gamma_h_sq_inv, E_Rqv_inv,
                                                          E_mu_h, Ji,
                                                          E_lambda_sq,
                                                          E_sigma_sq_inv,
                                                          E_tau_sq_inv,
                                                          E_sigma_site_sq_inv)
          # te <- difftime(Sys.time(), st)
          # print(paste0("Old version took: ", te))
          # print(E_S[1:3, 1:3])
          # 
          # st <- Sys.time()
          # longitudinal_covadj_siteadj_spatial_map_updates_fast(E_S, E_S_sq,
          #                                                 E_Beta, E_Beta_sq,
          #                                                 E_b_site, sum_E_b_site_sq,
          #                                                 E_Si, trace_E_Si_sq,
          #                                                 E_tr_BiBit_qv,
          #                                                 Y_unmixed,
          #                                                 X, Xrand, site_indicators,
          #                                                 pi_Zqv_eq_h, pi_Ztildeqv_eq_h,
          #                                                 E_gamma_h_sq_inv, E_Rqv_inv,
          #                                                 E_mu_h, Ji,
          #                                                 E_lambda_sq,
          #                                                 E_sigma_sq_inv,
          #                                                 E_tau_sq_inv,
          #                                                 E_sigma_site_sq_inv)
          # te <- difftime(Sys.time(), st)
          # print(paste0("New version took: ", te))
          # print(E_S[1:3, 1:3])
        } else {
          
          longitudinal_covadj_spatial_map_updates(E_S, E_S_sq,
                                                  E_Beta, E_Beta_sq,
                                                          E_Si, trace_E_Si_sq,
                                                          E_tr_BiBit_qv,
                                                          Y_unmixed,
                                                          X, Xrand,
                                                          pi_Zqv_eq_h, pi_Ztildeqv_eq_h,
                                                          E_gamma_h_sq_inv, E_Rqv_inv,
                                                          E_mu_h, Ji,
                                                          E_lambda_sq,
                                                          E_sigma_sq_inv,
                                                          E_tau_sq_inv)
          
        }
      } else {
        if (site_adjusted == TRUE){
          cross_sectional_covadj_siteadj_spatial_map_updates(
            E_S, E_S_sq, E_Beta, E_Beta_sq, E_b_site, sum_E_b_site_sq,
            E_Si, trace_E_Si_sq, Y_unmixed, X, site_indicators,
            pi_Zqv_eq_h, E_gamma_h_sq_inv, E_mu_h, 
            E_lambda_sq, E_sigma_sq_inv, E_tau_sq_inv, E_sigma_site_sq_inv) 
        } else {
          cross_sectional_covadj_spatial_map_updates(E_S, E_S_sq,
                                                     E_Beta, E_Beta_sq, E_Si,
                                                     trace_E_Si_sq,
                                                     Y_unmixed, X, pi_Zqv_eq_h,
                                                     E_gamma_h_sq_inv, E_mu_h,
                                                     E_lambda_sq, E_sigma_sq_inv,
                                                     E_tau_sq_inv) 
        }
      }
    }
    
    if (any(is.na(E_S))){
      stop("NA values found in source signals.")
    }
    if (any(is.na(E_Beta))){
      stop("NA values found in covariate effect maps.")
    }
    if (any(is.na(E_Beta_sq))){
      stop("NA values found in covariate effect squared expectations.")
    }
    elapsed_time <- Sys.time() - spatial_map_start_time
    time_S <- time_S + as.numeric(elapsed_time)

    ## Update S Stick Breaking Weights ----
    S_sticks_start_time <- Sys.time()
    probability_sums                     <- apply(pi_Zqv_eq_h, 1, sum) 
    cumulative_probability_sums_reversed <- rev(cumsum(c(0, rev(probability_sums))))[-1]
    VarHyper_A_nu_h <- probability_sums  + 1
    VarHyper_B_nu_h <- cumulative_probability_sums_reversed + concentration
    # Store the required expectations
    E_log_nu_h         <- digamma(VarHyper_A_nu_h) - digamma(VarHyper_A_nu_h + VarHyper_B_nu_h)
    E_log_1_minus_nu_h <- digamma(VarHyper_B_nu_h) - digamma(VarHyper_A_nu_h + VarHyper_B_nu_h)
    ### Timing
    elapsed_time <- Sys.time() - S_sticks_start_time
    time_S_sticks <- time_S_sticks + as.numeric(elapsed_time)
    
    ## Update S Cluster Parameters ----
    S_cluster_param_start_time <- Sys.time()
    calibrate_spatial_map_cluster_parameters(VarHyper_M_h, VarHyper_L_h,
                                             VarHyper_A_h, VarHyper_B_h,
                                             pi_Zqv_eq_h,
                                             E_S, E_S_sq, 
                                             E_sigma_sq_inv, base_measure_lambda,
                                             base_measure_mu, base_measure_alpha, 
                                             base_measure_beta)
    
    
    # Store the required expectations
    prev_mu_h        <- E_mu_h
    E_mu_h           <- VarHyper_M_h + 0 # keep it from direct assignment to M_h since updated in place
    E_mu_h_sq        <- VarHyper_B_h / ((VarHyper_A_h - 1)*VarHyper_L_h) + VarHyper_M_h^2
    prev_E_gamma_h_sq <- E_gamma_h_sq
    E_gamma_h_sq     <- VarHyper_B_h / (VarHyper_A_h - 1)
    # TODO CHECK IF ALREADY DOCUMENT I HAD SIGN WRONG HERE!
    E_log_gamma_h_sq <- -(digamma(VarHyper_A_h) - log(VarHyper_B_h)) 
    E_gamma_h_sq_inv <- VarHyper_A_h / VarHyper_B_h #1 / E_gamma_h_sq
    
    ### Relative Changes
    relchange_sigma_h_sq[, iter] <- abs(E_gamma_h_sq - prev_E_gamma_h_sq) / prev_E_gamma_h_sq
    relchange_mu_h[, iter]       <- abs(E_mu_h - prev_mu_h) / abs(prev_mu_h)
    
    print("Cluster means:")
    print(E_mu_h)
    print("Cluster vars:")
    print(E_gamma_h_sq)
    
    if (any(is.na( log(VarHyper_B_h) ))){
      n_na_SSq <- sum(is.na(E_S_sq))
      print(VarHyper_L_h) 
      print(VarHyper_M_h)
      print(VarHyper_A_h)
      print(VarHyper_B_h)
      stop("NA in VarHyper_B_h")
    }
    
    ### Timing
    elapsed_time <- Sys.time() - S_cluster_param_start_time
    time_S_cluster_param <- time_S_cluster_param + as.numeric(elapsed_time)
    
    ## Update RandEff Stick Breaking Weights ----
    if (longit_adjusted == TRUE){
      Rqv_sticks_start_time <- Sys.time()
      probability_sums                     <- apply(pi_Ztildeqv_eq_h, 1, sum) 
      cumulative_probability_sums_reversed <- rev(cumsum(c(0, rev(probability_sums))))[-1]
      VarHyper_A_Rqv_nu_h <- probability_sums  + 1
      VarHyper_B_Rqv_nu_h <- cumulative_probability_sums_reversed + concentration_Rqv
      # Store the required expectations
      E_Rqv_log_nu_h         <- digamma(VarHyper_A_Rqv_nu_h) - digamma(VarHyper_A_Rqv_nu_h + VarHyper_B_Rqv_nu_h)
      E_Rqv_log_1_minus_nu_h <- digamma(VarHyper_B_Rqv_nu_h) - digamma(VarHyper_A_Rqv_nu_h + VarHyper_B_Rqv_nu_h)
      ### Timing
      elapsed_time <- Sys.time() - Rqv_sticks_start_time
      # TODO store
      #time_S_sticks <- time_S_sticks + as.numeric(elapsed_time)
    }
    
    ## Update RandEff Cluster Parameters ----
    if (longit_adjusted == TRUE){
      calibrate_Rqv_cluster_parameters(VarHyper_Rqv_df, VarHyper_Rqv_scale,
                                       pi_Ztildeqv_eq_h,
                                       E_tr_BiBit_qv, base_measure_scale,
                                       E_sigma_sq_inv,
                                       base_measure_df, V, N)
      
      ### Update expected value of the log determinant of each cluster's 
      #   precision matrix
      for (h in 1:Hmax_Rqv){
        h_inds    <- ((h-1)*Prand + 1):(h*Prand)
        scale_mat <- VarHyper_Rqv_scale[, h_inds]
        
        # log determinant
        log_det_scale <- NA
        if (Prand == 1){
          log_det_scale <- log(scale_mat)
        } else {
          log_det_scale <- determinant(scale_mat, logarithm=TRUE)$modulus
        }
        
        # multivariate digamma function
        digamma_multivar <- sum(digamma((VarHyper_Rqv_df[h] + 1 - (1:Prand)) / 2))
        
        # # Only track the relevant part of the expression
        # if (Prand == 1){
        #   E_log_det_Rqv_inv[h] <- log((VarHyper_Rqv_scale[,h_inds]))
        # } else {
        #   E_log_det_Rqv_inv[h] <- log(det(VarHyper_Rqv_scale[,h_inds]))
        # }
        
        E_log_det_Rqv_inv[h] <- digamma_multivar + Prand * log(2) - log_det_scale
        
        # Plug in new expectations
        #E_Rqv_inv[,h_inds] = VarHyper_Rqv_scale[,h_inds] / (VarHyper_Rqv_df[h] - Prand - 1)
        #E_Rqv_inv[,h_inds] = VarHyper_Rqv_scale[,h_inds] / (VarHyper_Rqv_df[h] - Prand - 1)
        #E_Rqv_inv[,h_inds] = solve(VarHyper_Rqv_scale[,h_inds]) * (VarHyper_Rqv_df[h] - Prand - 1)
        E_Rqv_inv[,h_inds] = solve(VarHyper_Rqv_scale[,h_inds]) * (VarHyper_Rqv_df[h])
        
      }
      print("Estimates for precisions:")
      print(E_Rqv_inv)
      print("dfs:")
      print(VarHyper_Rqv_df)
    }

    
    ## Update the Rqv latent cluster memberships ----
    if (longit_adjusted == TRUE){
      calibrate_Rqv_cluster_probabilities(pi_Ztildeqv_eq_h, E_tr_BiBit_qv,  E_log_det_Rqv_inv,
                                          E_Rqv_inv, E_Rqv_log_nu_h, E_Rqv_log_1_minus_nu_h,
                                          E_sigma_sq_inv, V, N)
    }
    
    # print("Cheating for random effects clusters")
    # pi_Ztildeqv_eq_h[,] <- 0.001 / (nrow(pi_Ztildeqv_eq_h)-1)
    # for (iq in 1:Q){
    #   index_set <- ((iq-1)*V+1):(iq*V)
    #   cm_q <- c(true_values$true_Rqv_cluster_membership[iq, ])
    #   for (h in 1:nrow(pi_Ztildeqv_eq_h)){
    #     pi_Ztildeqv_eq_h[h, index_set[cm_q == h]] <- 0.999
    #   }
    # }
    # print("done cheating")
    # 
    ## Update Mixing Matrices ----
    mixing_matrix_start_time = Sys.time()
    update_all_mixing_matrices(E_Ai, Y_mixed_QxVNstar, E_Si, E_sigma_sq_inv)
    ### Timing
    elapsed_time <- Sys.time() - mixing_matrix_start_time
    time_mixing_matrix <- time_mixing_matrix + elapsed_time
    
    ### Re-demix the data using the newest mixing matrix
    apply_unmixing_matrix_start_time = Sys.time()
    apply_unmixing_matrices(Y_unmixed, Y_mixed_QxVNstar, E_Ai)
    ### Timing
    elapsed_time <- Sys.time() - apply_unmixing_matrix_start_time
    time_unmixing_Y_update <- time_unmixing_Y_update + elapsed_time
    
    
    ## Update Error Variance ----
    error_variance_start_time <- Sys.time()
    sum_trace_E_Si_sq <- sum(trace_E_Si_sq)
    if (covariate_adjusted == FALSE){
      VariHyper_B_sigma_sq <- update_error_variance_variational_scale(Y_unmixed,
                                                                      pi_Zqv_eq_h,
                                                                      E_S,
                                                                      E_Si,
                                                                      E_S_sq,
                                                                      E_mu_h,
                                                                      E_mu_h_sq,
                                                                      E_gamma_h_sq_inv,
                                                                      sum_trace_YYt,
                                                                      sum_trace_E_Si_sq)
    } else {
      VariHyper_B_sigma_sq <- update_error_variance_variational_scale_covadj(
        Y_unmixed,
        pi_Zqv_eq_h,
        E_S,
        E_Si,
        E_S_sq,
        E_mu_h,
        E_mu_h_sq,
        E_gamma_h_sq_inv,
        E_Beta_sq, E_lambda_sq,
        sum_trace_YYt, sum_trace_E_Si_sq,
        E_tau_sq_inv
      )
    }
    
    prev_E_sigma_sq <- E_sigma_sq
    E_sigma_sq_inv <- VariHyper_A_sigma_sq / VariHyper_B_sigma_sq
    E_sigma_sq     <- VariHyper_B_sigma_sq / (VariHyper_A_sigma_sq - 1)
    E_log_sigma_sq <- -(digamma(VariHyper_A_sigma_sq) - log(VariHyper_B_sigma_sq))
    ### Update relative change
    relchange_sigma_sq_error[iter] <- abs(E_sigma_sq-prev_E_sigma_sq)/prev_E_sigma_sq
    
    ### Timing for Error Variance Update
    elapsed_time <- Sys.time() - error_variance_start_time
    time_error_variance <- time_error_variance + elapsed_time
    
    ## Update Horseshoe Params ----
    
    horseshoe_hyper_start_time <- Sys.time()
    
    ### Global Shrinkage ----
    if (covariate_adjusted == TRUE){
      
      # Global Shrinkage
      VarHyper_A_tausq <- (Q*V*P + 1) / 2
      VarHyper_B_tausq <- 0.5 * E_sigma_sq_inv * sum(E_lambda_sq * E_Beta_sq) + E_xi_sq_inv
      prev_E_tau_sq <- E_tau_sq
      E_tau_sq     <- VarHyper_B_tausq / (VarHyper_A_tausq - 1)
      E_tau_sq_inv <- VarHyper_A_tausq / VarHyper_B_tausq
      relchange_tau_sq[iter] <- abs(E_tau_sq-prev_E_tau_sq)/prev_E_tau_sq
      
      # Hyperparameter for Global Shrinkage (xi)
      E_xi_sq_inv = 1 / (E_tau_sq_inv + 1)
      
      if (any(is.na( E_tau_sq_inv ))){
        stop("NA in E_tau_sq_inv")
      }
      
    }
    
    ### Local Shrinkage (inverted) ----
    if (covariate_adjusted == TRUE){
      Lqvp <- 0.5 * E_sigma_sq_inv * E_tau_sq_inv * E_Beta_sq
      E_lambda_sq <- 1 / (Lqvp * exp(Lqvp) * gsl::expint_E1(Lqvp)) - 1
      
      # exp * expint_E1 can run into numeric problems, here we resolve those
      na_lambda <- which(is.na(E_lambda_sq))
      if (length(na_lambda) > 0){
        print(paste0("Switching to CFA for ExpE1 for: ", length(na_lambda), " cases."))
        # evaluate expectation using continuous fraction approximation
        for (case in na_lambda){
          denom_term <- lentz_CFA(Lqvp[case])
          E_lambda_sq[case] <- 1 / (Lqvp[case] * denom_term) - 1
        }
        print(E_lambda_sq[na_lambda])
        # Check secondary failure cases
        zero_lambda <- which(E_lambda_sq == 0)
        if (length(zero_lambda) > 0){
          print(paste0("  CFA failed for: ", length(zero_lambda), " cases, using underflow fallback."))
          E_lambda_sq[zero_lambda] <- 1e-16
        }
      }
      
      if (any(is.na( E_lambda_sq ))){
        
        prop_na_E_Beta_sq <- sum(is.na(E_Beta_sq)) / length(E_Beta_sq)
        prop_na_E_lambda_sq <- sum(is.na(E_lambda_sq)) / length(E_lambda_sq)
        cat("Error variance: ", E_sigma_sq_inv, "\n")
        cat("Proportion of E beta sq that are na: ", prop_na_E_Beta_sq, "\n")
        cat("Proportion of E lambda sq that are na: ", prop_na_E_lambda_sq, "\n")
        cat("Summary of the Lqvp terms: \n")
        print(summary(c(Lqvp)))
        
        stop("NA in E_lambda_sq")
      }
      
    }
    
    ### timing for horseshoe hyperparameters
    elapsed_time <- Sys.time() - horseshoe_hyper_start_time
    time_horseshoe_hyper <- time_horseshoe_hyper + elapsed_time
    
    
    ## Update Site Variance ----
    site_var_start_time <- Sys.time()
    if (site_adjusted == TRUE){
      VarHyper_B_sigma_site_sq = sum_E_b_site_sq / 2 + 1.0
      VarHyper_A_sigma_site_sq = (Q*V*n_site+1)/2
      prev_E_sigma_site_sq <- E_sigma_site_sq
      E_sigma_site_sq     <- VarHyper_B_sigma_site_sq / (VarHyper_A_sigma_site_sq - 1)
      E_sigma_site_sq_inv <- VarHyper_A_sigma_site_sq / VarHyper_B_sigma_site_sq
      relchange_sigma_sq_site[iter] <- abs(E_sigma_site_sq-prev_E_sigma_site_sq) / prev_E_sigma_site_sq
    }
    elapsed_time <- Sys.time() - site_var_start_time
    time_site_var <- time_site_var + elapsed_time
    
    
    
  } # CAVI iterations
  
  # Format Output ----
  
  # Format beta as an array
  covariate_effect_estimates <- NULL
  local_precision_parameters <- NULL
  if (covariate_adjusted == TRUE){
    covariate_effect_estimates <- array(0, dim = c(V, Q, P))
    local_precision_parameters <- array(0, dim = c(V, Q, P))
    for (iq in 1:Q){
      for (ip in 1:P){
        covariate_effect_estimates[, iq, ip] <- E_Beta[ip, ((iq-1)*V+1):(iq * V) ]
        local_precision_parameters[, iq, ip] <- E_lambda_sq[ip, ((iq-1)*V+1):(iq * V) ]
      }
    }
  }
  
  
  # Format site effects as an array
  site_effect_estimates <- NULL
  if (site_adjusted == TRUE){
    site_effect_estimates <- array(0, dim = c(V, Q, n_site))
    for (iq in 1:Q){
      for (is in 1:n_site){
        site_effect_estimates[, iq, is] <- E_b_site[is, ((iq-1)*V+1):(iq * V) ]
      }
    }
  }
  
  # Split the results into lists containing different types of variables
  population_map_results = list(
    population_spatial_maps = E_S,
    cluster_probabilities   = pi_Zqv_eq_h,
    cluster_means           = E_mu_h,
    cluster_variances       = E_gamma_h_sq,
    cluster_precisions      = E_gamma_h_sq_inv
  )
  
  subject_map_results <- NULL
  if (store_subject_maps == TRUE){
    row_index  <- 1
    subject_map_results = array(0, dim = c(V, Q, N, Jmax))
    for (i in 1:N){
      for (j in 1:Ji[i]){
        column_set <- 1:V
        for (iq in 1:Q){
          subject_map_results[, iq, i, j] <- E_Si[row_index, column_set]
          column_set <- column_set + V # incr columns of E_Si
        }
        row_index  <- row_index + 1 # incr row of E_Si
      }
    }
  }
  
  covariate_effect_results = list(
    beta_hat            = covariate_effect_estimates,
    local_hs_precision  = E_lambda_sq,
    global_hs_variance  = E_tau_sq,
    global_hs_precision = E_tau_sq_inv
  )
  
  site_effect_results <- list(
    site_effects_hat = site_effect_estimates,
    site_effect_var = E_sigma_site_sq,
    site_effect_precision = E_sigma_site_sq_inv
  )
  
  longitudinal_results <- list(
    RE_precision_clusters = E_Rqv_inv,
    RE_cluster_probs = pi_Ztildeqv_eq_h
  )
  
  misc_terms = list(
    Y_unmixed  = Y_unmixed,
    X_with_int = X,
    Xrand      = Xrand,
    site_indicators = site_indicators
  )
  
  settings <- list(
    N = N,
    P = P,
    Q = Q,
    V = V,
    covariate_adjusted = covariate_adjusted,
    site_adjusted = site_adjusted,
    longit_adjusted = longit_adjusted,
    formula = formula,
    model_matrix_codings = model_matrix_codings,
    contrast_settings = contrast_settings
  )
  
  relative_change_tracking <- list(
    relchange_sigma_sq_error = relchange_sigma_sq_error,
    relchange_sigma_sq_site  = relchange_sigma_sq_site,
    relchange_tau_sq         = relchange_tau_sq,
    relchange_mu_h           = relchange_mu_h,
    relchange_sigma_h_sq     = relchange_sigma_h_sq
  )

  timing = tibble(
    update = c("Spatial Maps", "S Cluster Probs", "S Sticks", "S Cluster Params",
               "Mixing Matrix", "Unmixing Data", "Error Variance", "Horseshoe Hyperparameters",
               "Site Variance"),
    time   = c(time_S, time_S_cluster_probs, time_S_sticks, time_S_cluster_param,
               time_mixing_matrix, time_unmixing_Y_update, time_error_variance, 
               time_horseshoe_hyper, time_site_var)
  )
  
  result <- list(
    population_map_results   = population_map_results,
    subject_map_results      = subject_map_results,
    covariate_effect_results = covariate_effect_results,
    site_effect_results      = site_effect_results,
    longitudinal_results     = longitudinal_results,
    error_variance           = E_sigma_sq,
    error_precision          = E_sigma_sq_inv,
    settings = settings,
    misc_terms = misc_terms,
    relative_change_tracking = relative_change_tracking,
    timing = timing
  )
  
  class(result) <- "BBSSModel"

  return(result)
  
}