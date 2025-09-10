// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// Base version
// [[Rcpp::export]]
void cross_sectional_spatial_map_updates(Eigen::Map<Eigen::MatrixXd> & E_S,
                                         Eigen::Map<Eigen::MatrixXd> & E_S_sq,
                                         Eigen::Map<Eigen::MatrixXd> & E_Si,
                                         const Eigen::Map<Eigen::MatrixXd> & Y,
                                         const Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                                         const Eigen::Map<Eigen::VectorXd> & E_gamma_h_sq_inv,
                                         const Eigen::Map<Eigen::VectorXd> & E_mu_h,
                                         double E_sigma_sq_inv){
  
  // Dimensions
  size_t V = E_S.rows();
  size_t Q = E_S.cols();
  size_t N = Y.rows();
  size_t H = E_mu_h.size();
  
  // Initialize
  double prior_mean            = 0.0;
  double prior_precision       = 0.0;
  double variational_precision = 0.0;
  double variational_variance  = 0.0;
  double variational_mean      = 0.0;
  
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      
      // Get the prior mean and variance based on the DPM cluster weights and parameters
      prior_mean      = 0.0;
      prior_precision = 0.0;
      for (int h = 0; h < H; h++){
        prior_mean      += pi_Zqv_eq_h(h, iq*V + iv) * E_mu_h(h);
        prior_precision += pi_Zqv_eq_h(h, iq*V + iv) * E_gamma_h_sq_inv(h);
      }
      
      
      // Evaluate the variational distn precision
      variational_precision = E_sigma_sq_inv * N + prior_precision * E_sigma_sq_inv;
      
      // Invert to obtain variational distn variance
      variational_variance = 1.0 / variational_precision;
      
      // Get the variational mean
      variational_mean = variational_variance * 
        (Y.col(iq*V + iv).sum() * E_sigma_sq_inv + prior_mean * prior_precision * E_sigma_sq_inv);
      
      // Update S, S^2
      E_S(iv, iq)           = variational_mean;
      E_S_sq(iv, iq)        = pow(variational_mean, 2) + variational_variance;
      
      for (int i = 0; i < N; i++){
        E_Si(i, iq*V + iv) = E_S(iv, iq);
      }
      
    } // v
  } // q
  
}


// SITE         - Yes
// COVARIATE    - No
// Longitudinal - No
// [[Rcpp::export]]
void cross_sectional_siteadj_spatial_map_updates(Eigen::Map<Eigen::MatrixXd> & E_S,
                                                        Eigen::Map<Eigen::MatrixXd> & E_S_sq,
                                                        Eigen::Map<Eigen::MatrixXd> & E_b_site,
                                                        Eigen::Map<Eigen::VectorXd> & sum_E_b_site_sq,
                                                        Eigen::Map<Eigen::MatrixXd> & E_Si,
                                                        Eigen::Map<Eigen::MatrixXd> & trace_E_Si_sq,
                                                        const Eigen::Map<Eigen::MatrixXd> & Y,
                                                        const Eigen::Map<Eigen::MatrixXd> & site_indicators,
                                                        const Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                                                        const Eigen::Map<Eigen::VectorXd> & E_gamma_h_sq_inv,
                                                        const Eigen::Map<Eigen::VectorXd> & E_mu_h,
                                                        double E_sigma_sq_inv,
                                                        double E_sigma_site_sq_inv) {
  
  // Dimensions
  size_t V = E_S.rows();
  size_t Q = E_S.cols();
  size_t H = E_mu_h.size();
  size_t N = Y.rows();
  size_t nSite = E_b_site.rows();
  
  // Attach the site indicators to the rest of the design matrix
  Eigen::MatrixXd C_with_int = Eigen::MatrixXd::Zero(N, nSite + 1);
  for (int i = 0; i < N; i++){
    C_with_int(i,0) = 1.0;
  }
  C_with_int.rightCols(site_indicators.cols()) = site_indicators;
  
  // XtX for posterior precision
  Eigen::MatrixXd CtC = C_with_int.transpose() * C_with_int;
  
  // Zero out E_Si_sq term. This gets accumulated over in the loop below
  for (int i = 0; i < N; i++){
    trace_E_Si_sq(i) = 0.0;
  }
  
  // Zero out the term storing the sums of the bi_squared terms
  sum_E_b_site_sq(0) = 0.0;
  
  // Storage for intermediate quantities
  double s_prior_mean = 0.0;
  double s_prior_precision = 0.0;
  Eigen::MatrixXd prior_precision     = Eigen::MatrixXd::Zero(1 + nSite, 1 + nSite);
  Eigen::MatrixXd posterior_precision = Eigen::MatrixXd::Zero(1 + nSite, 1 + nSite);
  Eigen::MatrixXd posterior_variance  = Eigen::MatrixXd::Zero(1 + nSite, 1 + nSite);
  Eigen::VectorXd posterior_mean      = Eigen::VectorXd::Zero(1 + nSite);
  
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      
      // Get the prior mean and variance based on the DPM cluster weights and parameters
      s_prior_mean      = 0.0;
      s_prior_precision = 0.0;
      for (int h = 0; h < H; h++){
        s_prior_mean      += pi_Zqv_eq_h(h, iq*V + iv) * E_mu_h(h);
        s_prior_precision += pi_Zqv_eq_h(h, iq*V + iv) * E_gamma_h_sq_inv(h);
      }
      
      // Evaluate the prior contribution to the posterior precision
      prior_precision(0, 0) = E_sigma_sq_inv * s_prior_precision;
      // Site Effects
      for (int s = 0; s < nSite; s++){
        prior_precision(1 + s, 1 + s) = E_sigma_site_sq_inv;
      }
      
      // Evaluate the variational distn precision
      posterior_precision = E_sigma_sq_inv * CtC + prior_precision;
      
      // Invert to obtain variational distn variance
      posterior_variance = posterior_precision.inverse();
      
      // Get the variational mean
      posterior_mean = C_with_int.transpose() * Y.col(iq*V + iv) * E_sigma_sq_inv;
      
      // Increment with term coming from DPM on population source signals
      posterior_mean(0) += E_sigma_sq_inv * s_prior_precision * s_prior_mean;
      
      // Finally pre-multiply by posterior variance
      posterior_mean = posterior_variance * posterior_mean;
      
      // Update S, S^2, and Beta expectations
      E_S(iv, iq)              = posterior_mean(0);
      E_S_sq(iv, iq)           = pow(posterior_mean(0), 2) + posterior_variance(0, 0);
      for (int s = 0; s < nSite; s++){
        E_b_site(s, iq*V + iv)    = posterior_mean(1 + s);
        sum_E_b_site_sq(0) += pow(posterior_mean(1 + s), 2) + 
          posterior_variance(1 + s, 1 + s);
      }
      
      // Update subject-specific maps (linear combination of other variables)
      for (int i = 0; i < N; i++){
        E_Si(i, iq*V + iv) = (C_with_int.row(i) * posterior_mean)(0);
        trace_E_Si_sq(i)  += (C_with_int.row(i) * 
          (posterior_variance + posterior_mean*posterior_mean.transpose()) *
          C_with_int.row(i).transpose())(0);
      }
      
    } // v
  } // q
  
}



// E_S is V x Q
// E_Beta is P x QV one q at a time (IC 1, P x V) (IC 2, P x V) etc
// E_S_sq is V x Q
// Y is N x VQ
// covadj means covariates are included
// lambda_sq is a precision term, not variance. tau sq and sigma and gamma are all variance
// [[Rcpp::export]]
void cross_sectional_covadj_spatial_map_updates(Eigen::Map<Eigen::MatrixXd> & E_S,
                                                Eigen::Map<Eigen::MatrixXd> & E_S_sq,
                                                Eigen::Map<Eigen::MatrixXd> & E_Beta,
                                                Eigen::Map<Eigen::MatrixXd> & E_Beta_sq,
                                                Eigen::Map<Eigen::MatrixXd> & E_Si,
                                                Eigen::Map<Eigen::MatrixXd> & trace_E_Si_sq,
                                                const Eigen::Map<Eigen::MatrixXd> & Y,
                                                const Eigen::Map<Eigen::MatrixXd> & X,
                                                const Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                                                const Eigen::Map<Eigen::VectorXd> & E_gamma_h_sq_inv,
                                                const Eigen::Map<Eigen::VectorXd> & E_mu_h,
                                                const Eigen::Map<Eigen::MatrixXd> & E_lambda_sq,
                                                double E_sigma_sq_inv,
                                                double E_tau_sq_inv) {
  
  // Dimensions
  size_t V = E_S.rows();
  size_t Q = E_S.cols();
  size_t P = E_Beta.rows();
  size_t H = E_mu_h.size();
  size_t N = Y.rows();
  
  // XtX for posterior precision
  Eigen::MatrixXd XtX = X.transpose() * X;
  
  // Zero out E_Si_sq term. This gets accumulated over in the loop below
  for (int i = 0; i < N; i++){
    trace_E_Si_sq(i) = 0.0;
  }
  
  // Storage for intermediate quantities
  double s_prior_mean = 0.0;
  double s_prior_precision = 0.0;
  Eigen::MatrixXd prior_precision     = Eigen::MatrixXd::Zero(P + 1, P + 1);
  Eigen::MatrixXd posterior_precision = Eigen::MatrixXd::Zero(P + 1, P + 1);
  Eigen::MatrixXd posterior_variance  = Eigen::MatrixXd::Zero(P + 1, P + 1);
  Eigen::VectorXd posterior_mean      = Eigen::VectorXd::Zero(P + 1);
  
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      
      // Get the prior mean and variance based on the DPM cluster weights and parameters
      s_prior_mean      = 0.0;
      s_prior_precision = 0.0;
      for (int h = 0; h < H; h++){
        s_prior_mean      += pi_Zqv_eq_h(h, iq*V + iv) * E_mu_h(h);
        s_prior_precision += pi_Zqv_eq_h(h, iq*V + iv) * E_gamma_h_sq_inv(h);
      }
      
      // Evaluate the prior contribution to the posterior precision
      prior_precision(0, 0) = E_sigma_sq_inv * s_prior_precision;
      for (int p = 0; p < P; p++){
        prior_precision(p + 1, p + 1) = E_sigma_sq_inv * E_tau_sq_inv * E_lambda_sq(p, iq*V + iv);
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
      
      // Update S, S^2, and Beta expectations
      E_S(iv, iq)              = posterior_mean(0);
      E_S_sq(iv, iq)           = pow(posterior_mean(0), 2) + posterior_variance(0, 0);
      for (int p = 0; p < P; p++){
        E_Beta(p, iq*V + iv)    = posterior_mean(p + 1);
        E_Beta_sq(p, iq*V + iv) = pow(posterior_mean(p + 1), 2) + 
          posterior_variance(p + 1, p + 1);
      }
      
      // Update subject-specific maps (linear combination of other variables)
      for (int i = 0; i < N; i++){
        E_Si(i, iq*V + iv) = (X.row(i) * posterior_mean)(0);
        trace_E_Si_sq(i)  += (X.row(i) * 
          (posterior_variance + posterior_mean*posterior_mean.transpose()) *
          X.row(i).transpose())(0);
      }
      
    } // v
  } // q
  
}






// [[Rcpp::export]]
void cross_sectional_covadj_siteadj_spatial_map_updates(Eigen::Map<Eigen::MatrixXd> & E_S,
                                                        Eigen::Map<Eigen::MatrixXd> & E_S_sq,
                                                        Eigen::Map<Eigen::MatrixXd> & E_Beta,
                                                        Eigen::Map<Eigen::MatrixXd> & E_Beta_sq,
                                                        Eigen::Map<Eigen::MatrixXd> & E_b_site,
                                                        Eigen::Map<Eigen::VectorXd> & sum_E_b_site_sq,
                                                        Eigen::Map<Eigen::MatrixXd> & E_Si,
                                                        Eigen::Map<Eigen::MatrixXd> & trace_E_Si_sq,
                                                        const Eigen::Map<Eigen::MatrixXd> & Y,
                                                        const Eigen::Map<Eigen::MatrixXd> & X,
                                                        const Eigen::Map<Eigen::MatrixXd> & site_indicators,
                                                        const Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                                                        const Eigen::Map<Eigen::VectorXd> & E_gamma_h_sq_inv,
                                                        const Eigen::Map<Eigen::VectorXd> & E_mu_h,
                                                        const Eigen::Map<Eigen::MatrixXd> & E_lambda_sq,
                                                        double E_sigma_sq_inv,
                                                        double E_tau_sq_inv,
                                                        double E_sigma_site_sq_inv) {
  
  // Dimensions
  size_t V = E_S.rows();
  size_t Q = E_S.cols();
  size_t P = E_Beta.rows();
  size_t H = E_mu_h.size();
  size_t N = Y.rows();
  size_t nSite = E_b_site.rows();
  
  // Attach the site indicators to the rest of the design matrix
  Eigen::MatrixXd XC(X.rows(), X.cols() + site_indicators.cols());
  XC.leftCols(X.cols())                = X;
  XC.rightCols(site_indicators.cols()) = site_indicators;
  
  // XtX for posterior precision
  Eigen::MatrixXd XtX = XC.transpose() * XC;
  
  // Zero out E_Si_sq term. This gets accumulated over in the loop below
  for (int i = 0; i < N; i++){
    trace_E_Si_sq(i) = 0.0;
  }
  
  // Zero out the term storing the sums of the bi_squared terms
  sum_E_b_site_sq(0) = 0.0;
  
  // Storage for intermediate quantities
  double s_prior_mean = 0.0;
  double s_prior_precision = 0.0;
  Eigen::MatrixXd prior_precision     = Eigen::MatrixXd::Zero(P + 1 + nSite, P + 1 + nSite);
  Eigen::MatrixXd posterior_precision = Eigen::MatrixXd::Zero(P + 1 + nSite, P + 1 + nSite);
  Eigen::MatrixXd posterior_variance  = Eigen::MatrixXd::Zero(P + 1 + nSite, P + 1 + nSite);
  Eigen::VectorXd posterior_mean      = Eigen::VectorXd::Zero(P + 1 + nSite);
  
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      
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
      
      // Update S, S^2, and Beta expectations
      E_S(iv, iq)              = posterior_mean(0);
      E_S_sq(iv, iq)           = pow(posterior_mean(0), 2) + posterior_variance(0, 0);
      for (int p = 0; p < P; p++){
        E_Beta(p, iq*V + iv)    = posterior_mean(p + 1);
        E_Beta_sq(p, iq*V + iv) = pow(posterior_mean(p + 1), 2) + 
          posterior_variance(p + 1, p + 1);
      }
      for (int s = 0; s < nSite; s++){
        E_b_site(s, iq*V + iv)    = posterior_mean(P + 1 + s);
        sum_E_b_site_sq(0) += pow(posterior_mean(P + 1 + s), 2) + 
          posterior_variance(P + 1 + s, P + 1 + s);
      }
      
      // Update subject-specific maps (linear combination of other variables)
      for (int i = 0; i < N; i++){
        E_Si(i, iq*V + iv) = (XC.row(i) * posterior_mean)(0);
        trace_E_Si_sq(i)  += (XC.row(i) * 
          (posterior_variance + posterior_mean*posterior_mean.transpose()) *
          XC.row(i).transpose())(0);
      }
      
    } // v
  } // q
  
}


// This version blocks some of the inverse steps to speed up the computation by
// avoiding recomputing the parts of the inverse that are not changing from
// location to location
// [[Rcpp::export]]
void cross_sectional_covadj_siteadj_spatial_map_blocked_updates(Eigen::Map<Eigen::MatrixXd> & E_S,
                                                        Eigen::Map<Eigen::MatrixXd> & E_S_sq,
                                                        Eigen::Map<Eigen::MatrixXd> & E_Beta,
                                                        Eigen::Map<Eigen::MatrixXd> & E_Beta_sq,
                                                        Eigen::Map<Eigen::MatrixXd> & E_b_site,
                                                        Eigen::Map<Eigen::VectorXd> & sum_E_b_site_sq,
                                                        Eigen::Map<Eigen::MatrixXd> & E_Si,
                                                        Eigen::Map<Eigen::MatrixXd> & trace_E_Si_sq,
                                                        const Eigen::Map<Eigen::MatrixXd> & Y,
                                                        const Eigen::Map<Eigen::MatrixXd> & X,
                                                        const Eigen::Map<Eigen::MatrixXd> & site_indicators,
                                                        const Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                                                        const Eigen::Map<Eigen::VectorXd> & E_gamma_h_sq_inv,
                                                        const Eigen::Map<Eigen::VectorXd> & E_mu_h,
                                                        const Eigen::Map<Eigen::MatrixXd> & E_lambda_sq,
                                                        double E_sigma_sq_inv,
                                                        double E_tau_sq_inv,
                                                        double E_sigma_site_sq_inv) {
  
  // Dimensions
  size_t V = E_S.rows();
  size_t Q = E_S.cols();
  size_t P = E_Beta.rows();
  size_t H = E_mu_h.size();
  size_t N = Y.rows();
  size_t nSite = E_b_site.rows();
  
  // Attach the site indicators to the rest of the design matrix
  Eigen::MatrixXd XC(X.rows(), X.cols() + site_indicators.cols());
  XC.leftCols(X.cols())                = X;
  XC.rightCols(site_indicators.cols()) = site_indicators;
  
  // XtX for posterior precision
  Eigen::MatrixXd XptXp = X.transpose() * X;
  Eigen::MatrixXd XtX   = XC.transpose() * XC;
  
  
  
  // Zero out E_Si_sq term. This gets accumulated over in the loop below
  for (int i = 0; i < N; i++){
    trace_E_Si_sq(i) = 0.0;
  }
  
  // Zero out the term storing the sums of the bi_squared terms
  sum_E_b_site_sq(0) = 0.0;
  
  // Storage for intermediate quantities
  double s_prior_mean = 0.0;
  double s_prior_precision = 0.0;
  Eigen::MatrixXd prior_precision     = Eigen::MatrixXd::Zero(P + 1 + nSite, P + 1 + nSite);
  Eigen::MatrixXd posterior_precision = Eigen::MatrixXd::Zero(P + 1 + nSite, P + 1 + nSite);
  Eigen::MatrixXd posterior_variance  = Eigen::MatrixXd::Zero(P + 1 + nSite, P + 1 + nSite);
  Eigen::VectorXd posterior_mean      = Eigen::VectorXd::Zero(P + 1 + nSite);
  
  // Intermediate quantities related to block updates
  // Define the blocks
  Eigen::MatrixXd block_A = Eigen::MatrixXd::Zero(P + 1, P + 1);
  Eigen::MatrixXd block_A_minus_BDiC_inverse = Eigen::MatrixXd::Zero(P + 1, P + 1);
  Eigen::MatrixXd block_B = X.transpose() * site_indicators * E_sigma_sq_inv;
  Eigen::MatrixXd block_C = block_B.transpose();
  Eigen::MatrixXd block_D = site_indicators.transpose() * site_indicators * E_sigma_sq_inv;
  for (int s = 0; s < nSite; s++){
    block_D(s, s) += E_sigma_site_sq_inv;
  }
  // Pre-calculate terms for quick inverse
  Eigen::MatrixXd block_D_inverse = block_D.inverse();
  Eigen::MatrixXd block_BDiC      = block_B * block_D_inverse * block_C;
  Eigen::MatrixXd block_BDi       = block_B * block_D_inverse;
  Eigen::MatrixXd block_DiC       = block_D_inverse * block_C;
  
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      
      // Get the prior mean and variance based on the DPM cluster weights and parameters
      s_prior_mean      = 0.0;
      s_prior_precision = 0.0;
      for (int h = 0; h < H; h++){
        s_prior_mean      += pi_Zqv_eq_h(h, iq*V + iv) * E_mu_h(h);
        s_prior_precision += pi_Zqv_eq_h(h, iq*V + iv) * E_gamma_h_sq_inv(h);
      }
      
      // Evaluate the prior contribution to the posterior precision
      block_A *= 0.0;
      block_A(0, 0) = E_sigma_sq_inv * s_prior_precision;
      // Covariate Effects
      for (int p = 0; p < P; p++){
        block_A(p + 1, p + 1) = E_sigma_sq_inv * E_tau_sq_inv * E_lambda_sq(p, iq*V + iv);
      }
      block_A = block_A + XptXp * E_sigma_sq_inv;
      
      // Block inversion
      block_A_minus_BDiC_inverse = (block_A - block_BDiC).inverse();
      
      // Invert to obtain variational distn variance
      posterior_variance.block(0, 0, P+1, P+1)     = block_A_minus_BDiC_inverse;
      posterior_variance.block(0, P+1, P+1, nSite) = -1.0 * block_A_minus_BDiC_inverse * block_BDi;
      posterior_variance.block(P+1, 0, nSite, P+1) = posterior_variance.block(0, P+1, P+1, nSite).transpose();
      posterior_variance.block(P+1, P+1, nSite, nSite) = block_D_inverse + 
        block_DiC * block_A_minus_BDiC_inverse * block_BDi;
      
      // Get the variational mean
      posterior_mean = XC.transpose() * Y.col(iq*V + iv) * E_sigma_sq_inv;
      
      // Increment with term coming from DPM on population source signals
      posterior_mean(0) += E_sigma_sq_inv * s_prior_precision * s_prior_mean;
      
      // Finally pre-multiply by posterior variance
      posterior_mean = posterior_variance * posterior_mean;
      
      // Update S, S^2, and Beta expectations
      E_S(iv, iq)              = posterior_mean(0);
      E_S_sq(iv, iq)           = pow(posterior_mean(0), 2) + posterior_variance(0, 0);
      for (int p = 0; p < P; p++){
        E_Beta(p, iq*V + iv)    = posterior_mean(p + 1);
        E_Beta_sq(p, iq*V + iv) = pow(posterior_mean(p + 1), 2) + 
          posterior_variance(p + 1, p + 1);
      }
      for (int s = 0; s < nSite; s++){
        E_b_site(s, iq*V + iv)    = posterior_mean(P + 1 + s);
        sum_E_b_site_sq(0) += pow(posterior_mean(P + 1 + s), 2) + 
          posterior_variance(P + 1 + s, P + 1 + s);
      }
      
      // Update subject-specific maps (linear combination of other variables)
      for (int i = 0; i < N; i++){
        E_Si(i, iq*V + iv) = (XC.row(i) * posterior_mean)(0);
        trace_E_Si_sq(i)  += (XC.row(i) * 
          (posterior_variance + posterior_mean*posterior_mean.transpose()) *
          XC.row(i).transpose())(0);
      }
      
    } // v
  } // q
  
}
