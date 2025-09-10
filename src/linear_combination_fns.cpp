// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// 
// double estimate_prior_marginal_density(const Eigen::VectorXd& c_vec, int S = 10000, double tau_scale = 1.0, double lambda_scale = 1.0) {
//   int P = c_vec.size();
//   Eigen::VectorXd tau_samples(S);
//   Eigen::MatrixXd lambda_samples(S, P);
//   Eigen::VectorXd sigma2(S);
//   Eigen::VectorXd result(S);
//   double pi_val = M_PI;
//   
//   // Generate half-Cauchy samples using inverse transform
//   for (int s = 0; s < S; ++s) {
//     tau_samples[s] = tau_scale * tan(pi_val * R::runif(0, 0.5));
//     for (int j = 0; j < P; ++j) {
//       lambda_samples(s, j) = lambda_scale * tan(pi_val * R::runif(0, 0.5));
//     }
//   }
//   
//   for (int s = 0; s < S; ++s) {
//     double lambda_squared_weighted = 0.0;
//     for (int j = 0; j < P; ++j) {
//       lambda_squared_weighted += c_vec[j] * c_vec[j] * std::pow(lambda_samples(s, j), 2);
//     }
//     sigma2[s] = std::pow(tau_samples[s], 2) * lambda_squared_weighted;
//     result[s] = 1.0 / std::sqrt(2.0 * M_PI * sigma2[s]);
//   }
//   
//   // Monte Carlo average
//   return result.mean();
// }
//E_sigma_sq_inv
double estimate_prior_marginal_density(
    const Eigen::VectorXd& contrast_vector,
    const double sigma_sq,
    const int n_samples = 50000) {
  
  using Eigen::VectorXd;
  int P = contrast_vector.size();
  VectorXd z_samples(n_samples);
  
  for (int i = 0; i < n_samples; ++i) {
    // Sample xi^2 ~ IG(1/2, 1)
    double xi_sq = 1.0 / R::rgamma(0.5, 1.0);
    
    // Sample tau^2 | xi^2 ~ IG(1/2, 1/xi^2)
    double tau_sq = 1.0 / R::rgamma(0.5, 1.0 / xi_sq);
    
    // Sample lambda^2_j ~ horseshoe slab: pi(lam2) ∝ (lam2)^(-1/2)(1+lam2)^(-1)
    double variance_sum = 0.0;
    for (int j = 0; j < P; ++j) {
      double u = R::runif(0.0, 1.0);
      double lambda_j_sq = std::tan(M_PI * u / 2.0);  // inverse CDF of half-Cauchy^2
      variance_sum += contrast_vector(j) * contrast_vector(j) / lambda_j_sq;
    }
    
    double prior_var = sigma_sq * tau_sq * variance_sum;
    double dens = R::dnorm(0.0, 0.0, std::sqrt(prior_var), false);
    z_samples(i) = dens;
  }
  
  return z_samples.mean();
}

// [[Rcpp::export]]
Eigen::MatrixXd eval_LC_cross_sectional_covadj(const Eigen::Map<Eigen::MatrixXd> & contrasts,
                                               const Eigen::Map<Eigen::MatrixXd> & contrast_specification,
                                               const Eigen::Map<Eigen::MatrixXd> & Y,
                                               const Eigen::Map<Eigen::MatrixXd> & X,
                                               const Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                                               const Eigen::Map<Eigen::VectorXd> & E_gamma_h_sq_inv,
                                               const Eigen::Map<Eigen::VectorXd> & E_mu_h,
                                               const Eigen::Map<Eigen::MatrixXd> & E_lambda_sq,
                                               double E_sigma_sq_inv,
                                               double E_tau_sq_inv,
                                               int    V,
                                               bool calculate_credible_interval,
                                               double credible_interval_width_pct,
                                               bool calculate_savage_dickey,
                                               bool imp_sample_savage_dickey,
                                               int n_importance_resample) {
  
  // Dimensions
  size_t Q = Y.cols() / V;
  size_t P = contrasts.cols() - 1;
  size_t H = E_mu_h.size();
  size_t N = Y.rows();
  size_t n_contrast = contrasts.rows();
  
  // If this is true, will test evidence of contrast being zero
  //bool calculate_savage_dickey = true;
  // If this is true, will also use importance sampling to calculate IS version of ratio
  //bool imp_sample_savage_dickey = true;
  
  // XtX for posterior precision
  Eigen::MatrixXd XtX = X.transpose() * X;
  
  // Storage for intermediate quantities
  double s_prior_mean = 0.0;
  double s_prior_precision = 0.0;
  Eigen::MatrixXd prior_precision     = Eigen::MatrixXd::Zero(P + 1, P + 1);
  Eigen::MatrixXd prior_variance      = Eigen::MatrixXd::Zero(P + 1, P + 1); // added for savage-dickey
  double prior_variance_combination   = 0.0; // added for savage-dickey
  double prior_mean_combination       = 0.0; // added for savage-dickey
  double savage_dickey_ratio_denom    = 0.0; // added for savage-dickey
  double savage_dickey_ratio_numer    = 0.0; // added for savage-dickey
  Eigen::MatrixXd posterior_precision = Eigen::MatrixXd::Zero(P + 1, P + 1);
  Eigen::MatrixXd posterior_variance  = Eigen::MatrixXd::Zero(P + 1, P + 1);
  Eigen::VectorXd posterior_mean      = Eigen::VectorXd::Zero(P + 1);
  
  // Storage for final result
  int n_result_column = 1;
  if (calculate_credible_interval){
    n_result_column += 2;
  }
  if (calculate_savage_dickey){
    n_result_column += 1; // direct version
  }
  if (imp_sample_savage_dickey){
    n_result_column += 2; // adj imp samp version and corresponding ESS
  }
  int n_result_row = contrast_specification.rows();
  
  // This gets populated as we look through locations
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n_result_row, n_result_column);
  
  int iq   = 0;
  int iv   = 0;
  int irow = 0;
  double point_estimate = 0.0;
  double var_hat        = 0.0;
  double CI_llim = 0.0;
  double CI_ulim = 0.0;
  double Z_crit  = R::qnorm(1.0 - credible_interval_width_pct / 2.0, 0.0, 1.0, 1, 0);
  
  
  // Get the marginal prior for each contrast
  Eigen::VectorXd marginal_prior_contrast_probs = Eigen::VectorXd::Zero(n_contrast);
  if (calculate_savage_dickey == true){
    for (int icontrast = 0; icontrast < n_contrast; icontrast++){
      std::cout << "WARNING SAVAGE DICKEY CALCULATION WILL BE WRONG FOR CONTRAST WITH S0" << std::endl;
      marginal_prior_contrast_probs(icontrast) = estimate_prior_marginal_density(contrasts.row(icontrast), 1.0/E_sigma_sq_inv);
    }
    std::cout << "Marginal prior probabilities:" << std::endl;
    std::cout << marginal_prior_contrast_probs.transpose() << std::endl;
  }
  

  while (irow < n_result_row){
    
    // Get the settings
    // THESE NEED TO BE 0-based, this is why we subtract 1
    iq = contrast_specification(irow, 0) - 1 ;
    iv = contrast_specification(irow, 1) - 1;
    
    // Get the prior mean and variance based on the DPM cluster weights and parameters
    s_prior_mean      = 0.0;
    s_prior_precision = 0.0;
    for (int h = 0; h < H; h++){
      s_prior_mean      += pi_Zqv_eq_h(h, iq*V + iv) * E_mu_h(h);
      s_prior_precision += pi_Zqv_eq_h(h, iq*V + iv) * E_gamma_h_sq_inv(h);
    }
    
    // Evaluate the prior contribution to the posterior precision
    prior_precision(0, 0) = E_sigma_sq_inv * s_prior_precision;
    // Covariate Effects
    for (int p = 0; p < P; p++){
      prior_precision(p + 1, p + 1) = E_sigma_sq_inv * E_tau_sq_inv * E_lambda_sq(p, iq*V + iv);
    }
    
    if (calculate_savage_dickey){
      prior_variance = prior_precision.inverse();
    }
    
    // Evaluate the variational distn precision
    posterior_precision = E_sigma_sq_inv * XtX + prior_precision;
    
    // Invert to obtain variational distn variance
    posterior_variance = posterior_precision.inverse();
    
    // Get the variational mean
    posterior_mean = X.transpose() * Y.col(iq*V + iv) * E_sigma_sq_inv;
    
    // Increment with term coming from DPM on population source signals
    posterior_mean(0) += E_sigma_sq_inv * s_prior_precision * s_prior_mean;
    
    // Finally pre-multiply by posterior variance
    posterior_mean = posterior_variance * posterior_mean;
    
    // Update all contrasts
    for (int icontrast = 0; icontrast < n_contrast; icontrast++){
      
      // Evaluate the contrast
      point_estimate = (contrasts.row(icontrast) * posterior_mean.segment(0, P+1))(0);
      
      // Store the point estimate
      result(irow, 0) = point_estimate;
      
      if (calculate_credible_interval){
        
        // Evaluate the credible interval (if requested)
        var_hat = (contrasts.row(icontrast) *
          posterior_variance.block(0, 0, P+1, P+1) *
          contrasts.row(icontrast).transpose())(0);
        
        CI_llim = point_estimate - Z_crit * sqrt(var_hat);
        CI_ulim = point_estimate + Z_crit * sqrt(var_hat);
        
        result(irow, 1) = CI_llim;
        result(irow, 2) = CI_ulim;
        
      }  // end of CI check
      
      if (calculate_savage_dickey){
        
        // Only need to do this here if we did not do it earlier
        // for the credible interval calculation
        if (calculate_credible_interval == false){
          var_hat = (contrasts.row(icontrast) *
            posterior_variance.block(0, 0, P+1, P+1) *
            contrasts.row(icontrast).transpose())(0);
        }
        
        prior_variance_combination = (contrasts.row(icontrast) *
          prior_variance.block(0, 0, P+1, P+1) *
          contrasts.row(icontrast).transpose())(0);
        
        // Only element of prior mean that can be non-zero is the part
        // corresponding to S0
        prior_mean_combination = contrasts(icontrast, 0) * s_prior_mean;
        
        // Denominator for savage dickey ratio
        // savage_dickey_ratio_denom = R::dnorm4(0.0, prior_mean_combination,
        //                                       sqrt(prior_variance_combination), false);
        // std::cout << savage_dickey_ratio_denom << "   versus   " << marginal_prior_contrast_probs(icontrast) << std::endl;
        savage_dickey_ratio_denom = marginal_prior_contrast_probs(icontrast);
        
        // Numerator for savage dickey ratio
        savage_dickey_ratio_numer = R::dnorm4(0.0, point_estimate,
                                              sqrt(var_hat), false);
        
        
        result(irow, 3) = savage_dickey_ratio_numer / savage_dickey_ratio_denom;
        
        // Now importance-reweighted version
        
        if (imp_sample_savage_dickey == true){
          
          // this involves drawing samples from posterior, so going to pre-calculate
          // cholesky pieces
          Eigen::MatrixXd chol_Sigma = posterior_variance.llt().matrixL();
          Eigen::VectorXd theta_IS = Eigen::VectorXd::Zero(P+1);
          Eigen::VectorXd theta_IS_raw = Eigen::VectorXd::Zero(P+1); // before apply mean and var
          Eigen::VectorXd theta_centered = Eigen::VectorXd::Zero(P+1);
          Eigen::VectorXd resid_IS       = Eigen::VectorXd::Zero(N);
          int n_importance_resample = 1000;
          Eigen::VectorXd importance_weights = Eigen::VectorXd::Zero(n_importance_resample);
          // Need all because have to store them for bandwidth estimate
          Eigen::VectorXd all_IS_contrast_realizations = Eigen::VectorXd::Zero(n_importance_resample);
          Eigen::VectorXd all_IS_contrast_realizations_sorted = Eigen::VectorXd::Zero(n_importance_resample);
          double contrib_log_prior = 0.0;
          double contrib_log_lik   = 0.0;
          double contrib_log_vari  = 0.0;
          // Precalculate leading terms that do not depend on theta
          double log_prior_lead_term = log(prior_variance.determinant()) * (-1.0);
          double log_lik_lead_term = N/2.0 * log(E_sigma_sq_inv);
          double log_vari_lead_term  = log(posterior_variance.determinant()) * (-1.0);
          // Terms for kernel density estimate
          double ESS = 0.0;
          double IQR_Z = 0.0;
          double bandwidth = 0.0;
          double sd_Z = 0.0; // Standard deviation of contrast realizations
          double KDE_0 = 0.0;
          
          for (int i_IS = 0; i_IS < n_importance_resample; i_IS++){
            
            // Draw a sample from the variational posterior
            for (int p = 0; p < P + 1; p++){
              theta_IS_raw(p) = R::rnorm(0.0, 1.0);
            }
            theta_IS = posterior_mean + chol_Sigma * theta_IS_raw;
            
            // Calculate importance weights
            
            //// Log Prior
            theta_centered    = theta_IS; // Under prior, mean is 0, except element 11
            theta_centered(0) =  theta_centered(0) - s_prior_mean;
            contrib_log_prior = log_prior_lead_term - 0.5 * (
              theta_centered.transpose() * prior_precision * theta_centered)(0);
            
            //// Log Likelihood
            //X.transpose() * Y.col(iq*V + iv)
            resid_IS = Y.col(iq*V + iv) - X * theta_IS;
            contrib_log_lik   = log_lik_lead_term - 0.5 * E_sigma_sq_inv * (
              resid_IS.pow(2).sum()
            );
            
            //// Log Variational Approximation
            // Can just use std norm in exponential term
            contrib_log_vari  = log_vari_lead_term - 0.5 * theta_IS_raw.pow(2).sum();
            
            // Raw weights on log scale
            importance_weights(i_IS) = (contrib_log_prior + contrib_log_lik) - contrib_log_vari;
            
            // Store the contrast value for this realization
            all_IS_contrast_realizations(i_IS) = contrasts.row(icontrast) * theta_IS;
          }
          
          // Move the weights from log scale to regular scale
          double max_weight = importance_weights.maxCoeff();
          for (int iw = 0; iw < n_importance_resample; iw++){
            importance_weights(iw) -= max_weight;
            importance_weights(iw) = exp(importance_weights(iw));
          }
          double sum_weights = importance_weights.sum();
          for (int iw = 0; iw < n_importance_resample; iw++){
            importance_weights(iw) /= sum_weights;
          }
          
          // KDE
          ESS = pow(importance_weights.array().sum(), 2) / importance_weights.array().pow(2).sum();
          double weighted_mean = (importance_weights.array() * all_IS_contrast_realizations.array()).sum();
          double weighted_var = ((importance_weights.array() *
                                 (all_IS_contrast_realizations.array() - weighted_mean).square()).sum()) /
                                   (1.0 - importance_weights.array().square().sum()); // Effective denominator
          sd_Z = sqrt(weighted_var);
          //sd_Z = (all_IS_contrast_realizations.array() - all_IS_contrast_realizations.mean()).pow(2).sum() / (all_IS_contrast_realizations.size()-1);
          std::sort(all_IS_contrast_realizations.data(), 
                    all_IS_contrast_realizations.data() + all_IS_contrast_realizations.size()); 
          IQR_Z = all_IS_contrast_realizations(75) - all_IS_contrast_realizations(25);
          bandwidth = 0.9 * std::min(sd_Z, IQR_Z / 1.34) * pow(ESS, -0.2);
          KDE_0 = 0.0;
          for (int i_IS = 0; i_IS < n_importance_resample; i_IS++){
            KDE_0 += importance_weights(i_IS) * R::dnorm4(all_IS_contrast_realizations(i_IS), 0.0, bandwidth, false);
          }
          
          // From KDE we have posterior density, we evaluated prior dens earlier
          result(irow, 4) = KDE_0 / savage_dickey_ratio_denom;
          result(irow, 5) = ESS;
        } // end of importance resampling check
      } //  end savage dickey check
      
      // increment rows
      irow++;
      
    }// end of loop over contrasts for this location
    
  } // rows
  
  return(result);
  
}

// [[Rcpp::export]]
Eigen::MatrixXd eval_LC_cross_sectional_covadj_siteadj(const Eigen::Map<Eigen::MatrixXd> & contrasts,
                                                       const Eigen::Map<Eigen::MatrixXd> & contrast_specification,
                                                       const Eigen::Map<Eigen::MatrixXd> & Y,
                                                       const Eigen::Map<Eigen::MatrixXd> & X,
                                                       const Eigen::Map<Eigen::MatrixXd> & site_indicators,
                                                       const Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                                                       const Eigen::Map<Eigen::VectorXd> & E_gamma_h_sq_inv,
                                                       const Eigen::Map<Eigen::VectorXd> & E_mu_h,
                                                       const Eigen::Map<Eigen::MatrixXd> & E_lambda_sq,
                                                       double E_sigma_sq_inv,
                                                       double E_tau_sq_inv,
                                                       double E_sigma_site_sq_inv,
                                                       int    V,
                                                       bool calculate_credible_interval,
                                                       double credible_interval_width_pct,
                                                       bool calculate_savage_dickey,
                                                       bool imp_sample_savage_dickey,
                                                       int n_importance_resample) {
  
  // Dimensions
  size_t Q = Y.cols() / V;
  size_t P = contrasts.cols() - 1;
  size_t H = E_mu_h.size();
  size_t N = Y.rows();
  size_t nSite = site_indicators.cols();
  size_t n_contrast = contrasts.rows();
  
  // Attach the site indicators to the rest of the design matrix
  Eigen::MatrixXd XC(X.rows(), X.cols() + site_indicators.cols());
  XC.leftCols(X.cols())                = X;
  XC.rightCols(site_indicators.cols()) = site_indicators;
  
  // XtX for posterior precision
  Eigen::MatrixXd XtX = XC.transpose() * XC;
  
  // Storage for intermediate quantities
  double s_prior_mean = 0.0;
  double s_prior_precision = 0.0;
  Eigen::MatrixXd prior_precision     = Eigen::MatrixXd::Zero(P + 1 + nSite, P + 1 + nSite);
  Eigen::MatrixXd prior_variance      = Eigen::MatrixXd::Zero(P + 1 + nSite, P + 1 + nSite); // added for savage-dickey
  double prior_variance_combination   = 0.0; // added for savage-dickey
  double prior_mean_combination       = 0.0; // added for savage-dickey
  double savage_dickey_ratio_denom    = 0.0; // added for savage-dickey
  double savage_dickey_ratio_numer    = 0.0; // added for savage-dickey
  Eigen::MatrixXd posterior_precision = Eigen::MatrixXd::Zero(P + 1 + nSite, P + 1 + nSite);
  Eigen::MatrixXd posterior_variance  = Eigen::MatrixXd::Zero(P + 1 + nSite, P + 1 + nSite);
  Eigen::VectorXd posterior_mean      = Eigen::VectorXd::Zero(P + 1 + nSite);
  
  // Storage for final result
  int n_result_column = 1;
  if (calculate_credible_interval){
    n_result_column += 2;
  }
  if (calculate_savage_dickey){
    n_result_column += 1; // direct version
  }
  if (imp_sample_savage_dickey){
    n_result_column += 2; // adj imp samp version and corresponding ESS
  }
  int n_result_row = contrast_specification.rows();
  
  // This gets populated as we look through locations
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n_result_row, n_result_column);
  
  int iq   = 0;
  int iv   = 0;
  int irow = 0;
  double point_estimate = 0.0;
  double var_hat        = 0.0;
  double CI_llim = 0.0;
  double CI_ulim = 0.0;
  double Z_crit  = R::qnorm(1.0 - credible_interval_width_pct / 2.0, 0.0, 1.0, 1, 0);
  
  while (irow < n_result_row){
    
    // Get the settings
    // THESE NEED TO BE 0-based
    iq = contrast_specification(irow, 0) - 1 ;
    iv = contrast_specification(irow, 1) - 1;
    
    // Get the prior mean and variance based on the DPM cluster weights and parameters
    s_prior_mean      = 0.0;
    s_prior_precision = 0.0;
    for (int h = 0; h < H; h++){
      s_prior_mean      += pi_Zqv_eq_h(h, iq*V + iv) * E_mu_h(h);
      s_prior_precision += pi_Zqv_eq_h(h, iq*V + iv) * E_gamma_h_sq_inv(h);
    }
    
    // Evaluate the prior contribution to the posterior precision
    prior_precision(0, 0) = E_sigma_sq_inv * s_prior_precision;
    // Covariate Effects
    for (int p = 0; p < P; p++){
      prior_precision(p + 1, p + 1) = E_sigma_sq_inv * E_tau_sq_inv * E_lambda_sq(p, iq*V + iv);
    }
    // Site Effects
    for (int s = 0; s < nSite; s++){
      // TODO consider also including sigma sq error similar to beta terms
      // TODO speed up by caching lower right block and use block inverses
      prior_precision(P + 1 + s, P + 1 + s) = E_sigma_site_sq_inv;
    }
    
    // Calculate the prior variance for savage dickey ratio
    if (calculate_savage_dickey){
      prior_variance = prior_precision.inverse();
    }
    
    // Evaluate the variational distn precision
    posterior_precision = E_sigma_sq_inv * XtX + prior_precision;
    
    // Invert to obtain variational distn variance
    posterior_variance = posterior_precision.inverse();
    
    // Get the variational mean
    posterior_mean = XC.transpose() * Y.col(iq*V + iv) * E_sigma_sq_inv;
    
    // Increment with term coming from DPM on population source signals
    posterior_mean(0) += E_sigma_sq_inv * s_prior_precision * s_prior_mean;
    
    // Finally pre-multiply by posterior variance
    posterior_mean = posterior_variance * posterior_mean;
    
    // Update all contrasts
    for (int icontrast = 0; icontrast < n_contrast; icontrast++){
      
      // Evaluate the contrast
      point_estimate = (contrasts.row(icontrast) * posterior_mean.segment(0, P+1))(0);
      
      // if (irow > (n_result_row - 6) ){
      //   std::cout << "On row " << irow << " for contrast: " << icontrast << std::endl;
      //   std::cout << "PE " << point_estimate << std::endl;
      // }
      
      // Store the point estimate
      result(irow, 0) = point_estimate;
      
      if (calculate_credible_interval){
        
        // Evaluate the credible interval (if requested)
        var_hat = (contrasts.row(icontrast) *
          posterior_variance.block(0, 0, P+1, P+1) *
          contrasts.row(icontrast).transpose())(0);
        
        CI_llim = point_estimate - Z_crit * sqrt(var_hat);
        CI_ulim = point_estimate + Z_crit * sqrt(var_hat);
        
        result(irow, 1) = CI_llim;
        result(irow, 2) = CI_ulim;
        
      }  // end of CI check
      
      if (calculate_savage_dickey){
        
        // Only need to do this here if we did not do it earlier
        // for the credible interval calculation
        if (calculate_credible_interval == false){
          var_hat = (contrasts.row(icontrast) *
            posterior_variance.block(0, 0, P+1, P+1) *
            contrasts.row(icontrast).transpose())(0);
        }
        
        prior_variance_combination = (contrasts.row(icontrast) *
          prior_variance.block(0, 0, P+1, P+1) *
          contrasts.row(icontrast).transpose())(0);
        
        // Only element of prior mean that can be non-zero is the part
        // corresponding to S0
        prior_mean_combination = contrasts(icontrast, 0) * s_prior_mean;
        
        // Denominator for savage dickey ratio
        savage_dickey_ratio_denom = R::dnorm4(0.0, prior_mean_combination,
                                              sqrt(prior_variance_combination), false);
        
        // Numerator for savage dickey ratio
        savage_dickey_ratio_numer = R::dnorm4(0.0, point_estimate,
                                              sqrt(var_hat), false);
        
        result(irow, 3) = savage_dickey_ratio_numer / savage_dickey_ratio_denom;
        
        if ( std::isnan(result(irow, 3)) ){
          std::cout << "Num: " << savage_dickey_ratio_numer << 
            " and denom " << savage_dickey_ratio_denom <<
              " and log version of numer: " << R::dnorm4(0.0, point_estimate,
                              sqrt(var_hat), true) << std::endl;
          std::cout << "    when looking at the denominator we have a prior mean of " <<
            prior_mean_combination << " and a prior variance of : " << prior_variance_combination <<
              std::endl;
          
          std::cout << "The prior precision here was:" << std::endl;
          std::cout << prior_precision << std::endl;
          
          std::cout << "The lambda for this location is:" << std::endl;
          std::cout << E_lambda_sq.col(iq*V + iv) << std::endl;
        }
        
        // Now importance-reweighted version
        
        if (imp_sample_savage_dickey == true){
          
          // this involves drawing samples from posterior, so going to pre-calculate
          // cholesky pieces
          Eigen::MatrixXd chol_Sigma = posterior_variance.llt().matrixL();
          Eigen::VectorXd theta_IS = Eigen::VectorXd::Zero(P+1+nSite);
          Eigen::VectorXd theta_IS_raw = Eigen::VectorXd::Zero(P+1+nSite); // before apply mean and var
          Eigen::VectorXd theta_centered = Eigen::VectorXd::Zero(P+1+nSite);
          Eigen::VectorXd resid_IS       = Eigen::VectorXd::Zero(N);
          int n_importance_resample = 1000;
          Eigen::VectorXd importance_weights = Eigen::VectorXd::Zero(n_importance_resample);
          // Need all because have to store them for bandwidth estimate
          Eigen::VectorXd all_IS_contrast_realizations = Eigen::VectorXd::Zero(n_importance_resample);
          Eigen::VectorXd all_IS_contrast_realizations_sorted = Eigen::VectorXd::Zero(n_importance_resample);
          double contrib_log_prior = 0.0;
          double contrib_log_lik   = 0.0;
          double contrib_log_vari  = 0.0;
          // Precalculate leading terms that do not depend on theta
          double log_prior_lead_term = log(prior_variance.determinant()) * (-1.0);
          double log_lik_lead_term = N/2.0 * log(E_sigma_sq_inv);
          double log_vari_lead_term  = log(posterior_variance.determinant()) * (-1.0);
          // Terms for kernel density estimate
          double ESS = 0.0;
          double IQR_Z = 0.0;
          double bandwidth = 0.0;
          double sd_Z = 0.0; // Standard deviation of contrast realizations
          double KDE_0 = 0.0;
          
          for (int i_IS = 0; i_IS < n_importance_resample; i_IS++){
            
            // Draw a sample from the variational posterior
            for (int p = 0; p < P + 1 + nSite; p++){
              theta_IS_raw(p) = R::rnorm(0.0, 1.0);
            }
            theta_IS = posterior_mean + chol_Sigma * theta_IS_raw;
            
            // Calculate importance weights
            
            //// Log Prior
            theta_centered    = theta_IS; // Under prior, mean is 0, except element 11
            theta_centered(0) =  theta_centered(0) - s_prior_mean;
            contrib_log_prior = log_prior_lead_term - 0.5 * (
              theta_centered.transpose() * prior_precision * theta_centered)(0);
            
            //// Log Likelihood
            resid_IS = Y.col(iq*V + iv) - XC * theta_IS;
            contrib_log_lik   = log_lik_lead_term - 0.5 * E_sigma_sq_inv * (
              resid_IS.pow(2).sum()
            );
            
            //// Log Variational Approximation
            // Can just use std norm in exponential term
            contrib_log_vari  = log_vari_lead_term - 0.5 * theta_IS_raw.pow(2).sum();
            
            // Raw weights on log scale
            importance_weights(i_IS) = (contrib_log_prior + contrib_log_lik) - contrib_log_vari;
            
            // Store the contrast value for this realization
            // segment is used since the contrasts are only in terms of S and Beta
            all_IS_contrast_realizations(i_IS) = (contrasts.row(icontrast) * theta_IS.segment(0, P+1))(0);
          }
          
          // Move the weights from log scale to regular scale
          double max_weight = importance_weights.maxCoeff();
          for (int iw = 0; iw < n_importance_resample; iw++){
            importance_weights(iw) -= max_weight;
            importance_weights(iw) = exp(importance_weights(iw));
          }
          double sum_weights = importance_weights.sum();
          for (int iw = 0; iw < n_importance_resample; iw++){
            importance_weights(iw) /= sum_weights;
          }
          
          // KDE
          ESS = pow(importance_weights.array().sum(), 2) / importance_weights.array().pow(2).sum();
          double weighted_mean = (importance_weights.array() * all_IS_contrast_realizations.array()).sum();
          double weighted_var = ((importance_weights.array() *
                                 (all_IS_contrast_realizations.array() - weighted_mean).square()).sum()) /
                                   (1.0 - importance_weights.array().square().sum()); // Effective denominator
          sd_Z = sqrt(weighted_var);
          //sd_Z = (all_IS_contrast_realizations.array() - all_IS_contrast_realizations.mean()).pow(2).sum() / (all_IS_contrast_realizations.size()-1);
          std::sort(all_IS_contrast_realizations.data(), 
                    all_IS_contrast_realizations.data() + all_IS_contrast_realizations.size()); 
          IQR_Z = all_IS_contrast_realizations(75) - all_IS_contrast_realizations(25);
          bandwidth = 0.9 * std::min(sd_Z, IQR_Z / 1.34) * pow(ESS, -0.2);
          KDE_0 = 0.0;
          for (int i_IS = 0; i_IS < n_importance_resample; i_IS++){
            KDE_0 += importance_weights(i_IS) * R::dnorm4(all_IS_contrast_realizations(i_IS), 0.0, bandwidth, false);
          }
          
          // From KDE we have posterior density, we evaluated prior dens earlier
          result(irow, 4) = KDE_0 / savage_dickey_ratio_denom;
          result(irow, 5) = ESS;
        } // end of importance resampling check
      } //  end savage dickey check
      
      // increment rows
      irow++;
      
    }// end of loop over contrasts for this location
    
  } // rows
  
  return(result);
  
}








// [[Rcpp::export]]
Eigen::MatrixXd eval_LC_longitudinal_covadj(const Eigen::Map<Eigen::MatrixXd> & contrasts,
                                            const Eigen::Map<Eigen::MatrixXd> & contrast_specification,
                                            const Eigen::Map<Eigen::MatrixXd> & Y,
                                            const Eigen::Map<Eigen::MatrixXd> & X,
                                            const Eigen::Map<Eigen::MatrixXd> & Xrand,
                                            const Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                                            const Eigen::Map<Eigen::MatrixXd> & pi_Ztildeqv_eq_h,
                                            const Eigen::Map<Eigen::VectorXd> & E_gamma_h_sq_inv,
                                            const Eigen::Map<Eigen::VectorXd> & E_mu_h,
                                            const Eigen::Map<Eigen::MatrixXd> & E_lambda_sq,
                                            const Eigen::Map<Eigen::MatrixXd> & E_Rqv_h_inv,
                                            double E_sigma_sq_inv,
                                            double E_tau_sq_inv,
                                            int    V,
                                            bool calculate_credible_interval,
                                            double credible_interval_width_pct,
                                            bool calculate_savage_dickey,
                                            bool imp_sample_savage_dickey,
                                            int n_importance_resampl) {
  
  // Dimensions
  size_t Q          = Y.cols() / V;
  size_t P          = contrasts.cols() - 1;
  size_t Prand      = E_Rqv_h_inv.rows();
  size_t H          = E_mu_h.size();
  size_t NJ         = Y.rows();
  size_t N          = Xrand.cols() / Prand;
  size_t H_re       = E_Rqv_h_inv.cols() / Prand;
  size_t n_contrast = contrasts.rows();
  
  // Attach the random effects to the fixed effects design matrix
  Eigen::MatrixXd XCXtilde(X.rows(), X.cols() + Xrand.cols());
  XCXtilde.leftCols(X.cols())     = X;
  XCXtilde.rightCols(Xrand.cols()) = Xrand;
  
  // XtX for posterior precision
  Eigen::MatrixXd XtX = XCXtilde.transpose() * XCXtilde;
  
  // Storage for intermediate quantities
  double s_prior_mean = 0.0;
  double s_prior_precision = 0.0;
  Eigen::MatrixXd prior_precision     = Eigen::MatrixXd::Zero(P + 1  + N*Prand, P + 1  + N*Prand);
  Eigen::MatrixXd prior_variance      = Eigen::MatrixXd::Zero(P + 1 + N*Prand, P + 1 + N*Prand); // added for savage-dickey
  double prior_variance_combination   = 0.0; // added for savage-dickey
  double prior_mean_combination       = 0.0; // added for savage-dickey
  double savage_dickey_ratio_denom    = 0.0; // added for savage-dickey
  double savage_dickey_ratio_numer    = 0.0; // added for savage-dickey
  Eigen::MatrixXd posterior_precision = Eigen::MatrixXd::Zero(P + 1  + N*Prand, P + 1  + N*Prand);
  Eigen::MatrixXd posterior_variance  = Eigen::MatrixXd::Zero(P + 1  + N*Prand, P + 1  + N*Prand);
  Eigen::VectorXd posterior_mean      = Eigen::VectorXd::Zero(P + 1  + N*Prand);
  
  // Intermediate quantity storing the expected precision of the random effects for 
  // the specific network/location combination
  Eigen::MatrixXd random_eff_E_precision = Eigen::MatrixXd::Zero(Prand, Prand);
  
  // Storage for final result
  int n_result_column = 1;
  if (calculate_credible_interval){
    n_result_column += 2;
  }
  if (calculate_savage_dickey){
    n_result_column += 1; // direct version
  }
  if (imp_sample_savage_dickey){
    n_result_column += 2; // adj imp samp version and corresponding ESS
  }
  int n_result_row = contrast_specification.rows();
  
  // This gets populated as we look through locations
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n_result_row, n_result_column);
  
  int iq   = 0;
  int iv   = 0;
  int irow = 0;
  double point_estimate = 0.0;
  double var_hat        = 0.0;
  double CI_llim = 0.0;
  double CI_ulim = 0.0;
  double Z_crit  = R::qnorm(1.0 - credible_interval_width_pct / 2.0, 0.0, 1.0, 1, 0);
  
  while (irow < n_result_row){
    
    // Get the settings
    // THESE NEED TO BE 0-based
    iq = contrast_specification(irow, 0) - 1 ;
    iv = contrast_specification(irow, 1) - 1;
    
    // Get the prior mean and variance based on the DPM cluster weights and parameters
    s_prior_mean      = 0.0;
    s_prior_precision = 0.0;
    for (int h = 0; h < H; h++){
      s_prior_mean      += pi_Zqv_eq_h(h, iq*V + iv) * E_mu_h(h);
      s_prior_precision += pi_Zqv_eq_h(h, iq*V + iv) * E_gamma_h_sq_inv(h);
    }
    
    // Get the B_{qv} prior precision based on the DPM cluster weights and parameters
    random_eff_E_precision *= 0.0;
    for (int h_re = 0; h_re < H_re; h_re++){
      random_eff_E_precision += pi_Ztildeqv_eq_h(h_re, iq*V + iv) * 
        E_Rqv_h_inv.block(0, h_re*Prand, Prand, Prand);
    }
    
    // Evaluate the prior contribution to the posterior precision
    prior_precision(0, 0) = E_sigma_sq_inv * s_prior_precision;
    // Covariate Effects
    for (int p = 0; p < P; p++){
      prior_precision(p + 1, p + 1) = E_sigma_sq_inv * E_tau_sq_inv * E_lambda_sq(p, iq*V + iv);
    }
    for (int i = 0; i < N; i++){
      prior_precision.block(P + 1  + Prand*i,
                            P + 1  + Prand*i,
                            Prand, Prand) = random_eff_E_precision;
    }
    
    // Calculate the prior variance for savage dickey ratio
    if (calculate_savage_dickey){
      prior_variance = prior_precision.inverse();
    }
    
    // Evaluate the variational distn precision
    posterior_precision = E_sigma_sq_inv * XtX + prior_precision;
    
    // Invert to obtain variational distn variance
    posterior_variance = posterior_precision.inverse();
    
    // Get the variational mean
    posterior_mean = XCXtilde.transpose() * Y.col(iq*V + iv) * E_sigma_sq_inv;
    
    // Increment with term coming from DPM on population source signals
    posterior_mean(0) += E_sigma_sq_inv * s_prior_precision * s_prior_mean;
    
    // Finally pre-multiply by posterior variance
    posterior_mean = posterior_variance * posterior_mean;
    
    // Update all contrasts
    for (int icontrast = 0; icontrast < n_contrast; icontrast++){
      
      // Evaluate the contrast
      point_estimate = (contrasts.row(icontrast) * posterior_mean.segment(0, P+1))(0);
      
      // Store the point estimate
      result(irow, 0) = point_estimate;
      
      if (calculate_credible_interval){
        
        // Evaluate the credible interval (if requested)
        var_hat = (contrasts.row(icontrast) *
          posterior_variance.block(0, 0, P+1, P+1) *
          contrasts.row(icontrast).transpose())(0);
        
        CI_llim = point_estimate - Z_crit * sqrt(var_hat);
        CI_ulim = point_estimate + Z_crit * sqrt(var_hat);
        
        result(irow, 1) = CI_llim;
        result(irow, 2) = CI_ulim;
        
      }  // end of CI check
      
      if (calculate_savage_dickey){
        
        // Only need to do this here if we did not do it earlier
        // for the credible interval calculation
        if (calculate_credible_interval == false){
          var_hat = (contrasts.row(icontrast) *
            posterior_variance.block(0, 0, P+1, P+1) *
            contrasts.row(icontrast).transpose())(0);
        }
        
        prior_variance_combination = (contrasts.row(icontrast) *
          prior_variance.block(0, 0, P+1, P+1) *
          contrasts.row(icontrast).transpose())(0);
        
        // Only element of prior mean that can be non-zero is the part
        // corresponding to S0
        prior_mean_combination = contrasts(icontrast, 0) * s_prior_mean;
        
        // Denominator for savage dickey ratio
        savage_dickey_ratio_denom = R::dnorm4(0.0, prior_mean_combination,
                                              sqrt(prior_variance_combination), false);
        
        // Numerator for savage dickey ratio
        savage_dickey_ratio_numer = R::dnorm4(0.0, point_estimate,
                                              sqrt(var_hat), false);
        
        result(irow, 3) = savage_dickey_ratio_numer / savage_dickey_ratio_denom;
        
        if ( std::isnan(result(irow, 3)) ){
          std::cout << "Num: " << savage_dickey_ratio_numer << 
            " and denom " << savage_dickey_ratio_denom <<
              " and log version of numer: " << R::dnorm4(0.0, point_estimate,
                              sqrt(var_hat), true) << std::endl;
          std::cout << "    when looking at the denominator we have a prior mean of " <<
            prior_mean_combination << " and a prior variance of : " << prior_variance_combination <<
              std::endl;
          
          std::cout << "The prior precision here was:" << std::endl;
          std::cout << prior_precision << std::endl;
          
          std::cout << "The lambda for this location is:" << std::endl;
          std::cout << E_lambda_sq.col(iq*V + iv) << std::endl;
        }
        
        // Now importance-reweighted version
        
        if (imp_sample_savage_dickey == true){
          
          // this involves drawing samples from posterior, so going to pre-calculate
          // cholesky pieces
          Eigen::MatrixXd chol_Sigma = posterior_variance.llt().matrixL();
          Eigen::VectorXd theta_IS = Eigen::VectorXd::Zero(P+1+N*Prand);
          Eigen::VectorXd theta_IS_raw = Eigen::VectorXd::Zero(P+1+N*Prand); // before apply mean and var
          Eigen::VectorXd theta_centered = Eigen::VectorXd::Zero(P+1+N*Prand);
          Eigen::VectorXd resid_IS       = Eigen::VectorXd::Zero(N);
          int n_importance_resample = 1000;
          Eigen::VectorXd importance_weights = Eigen::VectorXd::Zero(n_importance_resample);
          // Need all because have to store them for bandwidth estimate
          Eigen::VectorXd all_IS_contrast_realizations = Eigen::VectorXd::Zero(n_importance_resample);
          Eigen::VectorXd all_IS_contrast_realizations_sorted = Eigen::VectorXd::Zero(n_importance_resample);
          double contrib_log_prior = 0.0;
          double contrib_log_lik   = 0.0;
          double contrib_log_vari  = 0.0;
          // Precalculate leading terms that do not depend on theta
          double log_prior_lead_term = log(prior_variance.determinant()) * (-1.0);
          double log_lik_lead_term = N/2.0 * log(E_sigma_sq_inv);
          double log_vari_lead_term  = log(posterior_variance.determinant()) * (-1.0);
          // Terms for kernel density estimate
          double ESS = 0.0;
          double IQR_Z = 0.0;
          double bandwidth = 0.0;
          double sd_Z = 0.0; // Standard deviation of contrast realizations
          double KDE_0 = 0.0;
          
          for (int i_IS = 0; i_IS < n_importance_resample; i_IS++){
            
            // Draw a sample from the variational posterior
            for (int p = 0; p < P + 1 + N*Prand; p++){
              theta_IS_raw(p) = R::rnorm(0.0, 1.0);
            }
            theta_IS = posterior_mean + chol_Sigma * theta_IS_raw;
            
            // Calculate importance weights
            
            //// Log Prior
            theta_centered    = theta_IS; // Under prior, mean is 0, except element 11
            theta_centered(0) =  theta_centered(0) - s_prior_mean;
            contrib_log_prior = log_prior_lead_term - 0.5 * (
              theta_centered.transpose() * prior_precision * theta_centered)(0);
            
            //// Log Likelihood
            resid_IS = Y.col(iq*V + iv) - XCXtilde * theta_IS;
            contrib_log_lik   = log_lik_lead_term - 0.5 * E_sigma_sq_inv * (
              resid_IS.pow(2).sum()
            );
            
            //// Log Variational Approximation
            // Can just use std norm in exponential term
            contrib_log_vari  = log_vari_lead_term - 0.5 * theta_IS_raw.pow(2).sum();
            
            // Raw weights on log scale
            importance_weights(i_IS) = (contrib_log_prior + contrib_log_lik) - contrib_log_vari;
            
            // Store the contrast value for this realization
            // segment is used since the contrasts are only in terms of S and Beta
            all_IS_contrast_realizations(i_IS) = (contrasts.row(icontrast) * theta_IS.segment(0, P+1))(0);
          }
          
          // Move the weights from log scale to regular scale
          double max_weight = importance_weights.maxCoeff();
          for (int iw = 0; iw < n_importance_resample; iw++){
            importance_weights(iw) -= max_weight;
            importance_weights(iw) = exp(importance_weights(iw));
          }
          double sum_weights = importance_weights.sum();
          for (int iw = 0; iw < n_importance_resample; iw++){
            importance_weights(iw) /= sum_weights;
          }
          
          // KDE
          ESS = pow(importance_weights.array().sum(), 2) / importance_weights.array().pow(2).sum();
          double weighted_mean = (importance_weights.array() * all_IS_contrast_realizations.array()).sum();
          double weighted_var = ((importance_weights.array() *
                                 (all_IS_contrast_realizations.array() - weighted_mean).square()).sum()) /
                                   (1.0 - importance_weights.array().square().sum()); // Effective denominator
          sd_Z = sqrt(weighted_var);
          //sd_Z = (all_IS_contrast_realizations.array() - all_IS_contrast_realizations.mean()).pow(2).sum() / (all_IS_contrast_realizations.size()-1);
          std::sort(all_IS_contrast_realizations.data(), 
                    all_IS_contrast_realizations.data() + all_IS_contrast_realizations.size()); 
          IQR_Z = all_IS_contrast_realizations(75) - all_IS_contrast_realizations(25);
          bandwidth = 0.9 * std::min(sd_Z, IQR_Z / 1.34) * pow(ESS, -0.2);
          KDE_0 = 0.0;
          for (int i_IS = 0; i_IS < n_importance_resample; i_IS++){
            KDE_0 += importance_weights(i_IS) * R::dnorm4(all_IS_contrast_realizations(i_IS), 0.0, bandwidth, false);
          }
          
          // From KDE we have posterior density, we evaluated prior dens earlier
          result(irow, 4) = KDE_0 / savage_dickey_ratio_denom;
          result(irow, 5) = ESS;
        } // end of importance resampling check
      } //  end savage dickey check
      
      // increment rows
      irow++;
      
    }// end of loop over contrasts for this location
    
  } // rows
  
  return(result);
  
}
