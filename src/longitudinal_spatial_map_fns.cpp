// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

//pi_Ztildeqv_eq_h is cluster memberships for the subject effects
// [[Rcpp::export]]
void longitudinal_covadj_siteadj_spatial_map_updates(Eigen::Map<Eigen::MatrixXd> & E_S,
                                                     Eigen::Map<Eigen::MatrixXd> & E_S_sq,
                                                     Eigen::Map<Eigen::MatrixXd> & E_Beta,
                                                     Eigen::Map<Eigen::MatrixXd> & E_Beta_sq,
                                                     Eigen::Map<Eigen::MatrixXd> & E_b_site,
                                                     Eigen::Map<Eigen::VectorXd> & sum_E_b_site_sq,
                                                     Eigen::Map<Eigen::MatrixXd> & E_Si,
                                                     Eigen::Map<Eigen::MatrixXd> & trace_E_Si_sq,
                                                     Eigen::Map<Eigen::MatrixXd> & E_tr_BiBit_qv,
                                                     const Eigen::Map<Eigen::MatrixXd> & Y,
                                                     const Eigen::Map<Eigen::MatrixXd> & X,
                                                     const Eigen::Map<Eigen::MatrixXd> & Xrand,
                                                     const Eigen::Map<Eigen::MatrixXd> & site_indicators,
                                                     const Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                                                     const Eigen::Map<Eigen::MatrixXd> & pi_Ztildeqv_eq_h,
                                                     const Eigen::Map<Eigen::VectorXd> & E_gamma_h_sq_inv,
                                                     const Eigen::Map<Eigen::MatrixXd> & E_Rqv_h_inv,
                                                     const Eigen::Map<Eigen::VectorXd> & E_mu_h,
                                                     const Eigen::Map<Eigen::VectorXd> & Ji,
                                                     const Eigen::Map<Eigen::MatrixXd> & E_lambda_sq,
                                                     double E_sigma_sq_inv,
                                                     double E_tau_sq_inv,
                                                     double E_sigma_site_sq_inv) {
  
  // Dimensions
  size_t V     = E_S.rows();
  size_t Q     = E_S.cols();
  size_t P     = E_Beta.rows();
  size_t Prand = E_Rqv_h_inv.rows();
  size_t H     = E_mu_h.size();
  size_t NJ    = Y.rows();
  size_t N     = Xrand.cols() / Prand;
  size_t nSite = E_b_site.rows();
  size_t H_re  = E_Rqv_h_inv.cols() / Prand;
  
  //std::cout << "N = " << N << std::endl;

  // Attach the site indicators to the rest of the design matrix
  Eigen::MatrixXd XC(X.rows(), X.cols() + site_indicators.cols());
  XC.leftCols(X.cols())                = X;
  XC.rightCols(site_indicators.cols()) = site_indicators;
  
  // Attach the random effects to the fixed + site effects design matrix
  Eigen::MatrixXd XCXtilde(XC.rows(), XC.cols() + Xrand.cols());
  XCXtilde.leftCols(XC.cols())     = XC;
  XCXtilde.rightCols(Xrand.cols()) = Xrand;
  
  // XtX for posterior precision
  Eigen::MatrixXd XtX = XCXtilde.transpose() * XCXtilde;
  
  // Zero out E_Si_sq term. This gets accumulated over in the loop below
  for (int i = 0; i < NJ; i++){
    trace_E_Si_sq(i) = 0.0;
  }
  
  // Zero out the term storing the sums of the bi_squared terms
  sum_E_b_site_sq(0) = 0.0;
  
  // This version is needed for the cluster membership updates
  for (int j = 0; j < E_tr_BiBit_qv.cols(); j++){
    for (int i = 0; i < E_tr_BiBit_qv.rows(); i++){
      E_tr_BiBit_qv(i, j) = 0.0;
    }
  }

  // Storage for intermediate quantities
  double s_prior_mean = 0.0;
  double s_prior_precision = 0.0;
  Eigen::MatrixXd prior_precision     = Eigen::MatrixXd::Zero(P + 1 + nSite + N*Prand, P + 1 + nSite + N*Prand);
  Eigen::MatrixXd posterior_precision = Eigen::MatrixXd::Zero(P + 1 + nSite + N*Prand, P + 1 + nSite + N*Prand);
  Eigen::MatrixXd posterior_variance  = Eigen::MatrixXd::Zero(P + 1 + nSite + N*Prand, P + 1 + nSite + N*Prand);
  Eigen::VectorXd posterior_mean      = Eigen::VectorXd::Zero(P + 1 + nSite + N*Prand);
  
  // Intermediate quantity storing the expected precision of the random effects for 
  // the specific network/location combination
  Eigen::MatrixXd random_eff_E_precision = Eigen::MatrixXd::Zero(Prand, Prand);
  Eigen::MatrixXd E_BiBit                = Eigen::MatrixXd::Zero(Prand, Prand);
  
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      
      // Get the S prior mean and variance based on the DPM cluster weights and parameters
      s_prior_mean      = 0.0;
      s_prior_precision = 0.0;
      for (int h = 0; h < H; h++){
        s_prior_mean      += pi_Zqv_eq_h(h, iq*V + iv) * E_mu_h(h);
        s_prior_precision += pi_Zqv_eq_h(h, iq*V + iv) * E_gamma_h_sq_inv(h);
      }
      
      // if (iq == 0){
      //   if (iv == 0){
      //     std::cout << "S prior mean and precision: " << std::endl;
      //     std::cout << s_prior_mean << "   and  " << s_prior_precision << std::endl;
      //   }
      // }
      
      
      // Get the B_{qv} prior precision based on the DPM cluster weights and parameters
      random_eff_E_precision *= 0.0;
      for (int h_re = 0; h_re < H_re; h_re++){
        random_eff_E_precision += pi_Ztildeqv_eq_h(h_re, iq*V + iv) * 
          E_Rqv_h_inv.block(0, h_re*Prand, Prand, Prand);
      }
      
      // if (iq == 0){
      //   if (iv == 0){
      //     std::cout << "The precision of the random effects for this element is: " << std::endl;
      //     std::cout << random_eff_E_precision << std::endl;
      //   }
      // }
      
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
      for (int i = 0; i < N; i++){
        prior_precision.block(P + 1 + nSite + Prand*i,
                              P + 1 + nSite + Prand*i,
                              Prand, Prand) = random_eff_E_precision;
      }
      
      // if (iq == 0){
      //   if (iv == 0){
      //     std::cout << "The prior precision for this element is: " << std::endl;
      //     std::cout << prior_precision.block(0, 0, P + 1 + nSite + Prand*3, P + 1 + nSite + Prand*3) << std::endl;
      //   }
      // }
      
      // Evaluate the variational distn precision
      posterior_precision = E_sigma_sq_inv * XtX + prior_precision;
      // 
      // if (iq == 0){
      //   if (iv == 0){
      //     std::cout << "The posterior precision for this element is: " << std::endl;
      //     std::cout << posterior_precision.block(0, 0, P + 1 + nSite + Prand*2, P + 1 + nSite + Prand*2) << std::endl;
      //   }
      // }
      
      // Invert to obtain variational distn variance
      posterior_variance = posterior_precision.inverse();
      
      // if (iq == 0){
      //   if (iv == 0){
      //     std::cout << "The posterior variance for this element is: " << std::endl;
      //     std::cout << posterior_variance.block(0, 0, P + 1 + nSite + Prand*2, P + 1 + nSite + Prand*2) << std::endl;
      //   }
      // }
      
      // Get the variational mean
      posterior_mean = XCXtilde.transpose() * Y.col(iq*V + iv) * E_sigma_sq_inv;
      
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
      // We don't store Biqv, but we do need some summary statistics based on it
      for (int i = 0; i < N; i++){
        E_BiBit = posterior_variance.block(P + 1 + nSite + i*Prand,
                                           P + 1 + nSite + i*Prand,
                                           Prand, Prand) + 
                                             posterior_mean.segment(P + 1 + nSite + i*Prand, Prand) *
                                             posterior_mean.segment(P + 1 + nSite + i*Prand, Prand).transpose();
        E_tr_BiBit_qv.block(0, iq*V*Prand + iv*Prand, Prand, Prand) += E_BiBit;
      }
      
      // Update subject-specific maps (linear combination of other variables)
      // TODO precalc posterior_variance + posterior_mean*posterior_mean.transpose()
      Eigen::MatrixXd middle_term = (posterior_variance + posterior_mean*posterior_mean.transpose());
      for (int ij = 0; ij < NJ; ij++){
        E_Si(ij, iq*V + iv) = (XCXtilde.row(ij) * posterior_mean)(0);
        trace_E_Si_sq(ij)  += (XCXtilde.row(ij) * 
          middle_term *
          XCXtilde.row(ij).transpose())(0);
      }
      
    } // v
  } // q
  
}





// [[Rcpp::export]]
void longitudinal_covadj_siteadj_spatial_map_updates_fast(Eigen::Map<Eigen::MatrixXd> & E_S,
                                                     Eigen::Map<Eigen::MatrixXd> & E_S_sq,
                                                     Eigen::Map<Eigen::MatrixXd> & E_Beta,
                                                     Eigen::Map<Eigen::MatrixXd> & E_Beta_sq,
                                                     Eigen::Map<Eigen::MatrixXd> & E_b_site,
                                                     Eigen::Map<Eigen::VectorXd> & sum_E_b_site_sq,
                                                     Eigen::Map<Eigen::MatrixXd> & E_Si,
                                                     Eigen::Map<Eigen::MatrixXd> & trace_E_Si_sq,
                                                     Eigen::Map<Eigen::MatrixXd> & E_tr_BiBit_qv,
                                                     const Eigen::Map<Eigen::MatrixXd> & Y,
                                                     const Eigen::Map<Eigen::MatrixXd> & X,
                                                     const Eigen::Map<Eigen::MatrixXd> & Xrand,
                                                     const Eigen::Map<Eigen::MatrixXd> & site_indicators,
                                                     const Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                                                     const Eigen::Map<Eigen::MatrixXd> & pi_Ztildeqv_eq_h,
                                                     const Eigen::Map<Eigen::VectorXd> & E_gamma_h_sq_inv,
                                                     const Eigen::Map<Eigen::MatrixXd> & E_Rqv_h_inv,
                                                     const Eigen::Map<Eigen::VectorXd> & E_mu_h,
                                                     const Eigen::Map<Eigen::VectorXd> & Ji,
                                                     const Eigen::Map<Eigen::MatrixXd> & E_lambda_sq,
                                                     double E_sigma_sq_inv,
                                                     double E_tau_sq_inv,
                                                     double E_sigma_site_sq_inv) {
  
  // Dimensions
  size_t V     = E_S.rows();
  size_t Q     = E_S.cols();
  size_t P     = E_Beta.rows();
  size_t Prand = E_Rqv_h_inv.rows();
  size_t H     = E_mu_h.size();
  size_t NJ    = Y.rows();
  size_t N     = Xrand.cols() / Prand;
  size_t nSite = E_b_site.rows();
  size_t H_re  = E_Rqv_h_inv.cols() / Prand;
  
  //std::cout << "N = " << N << std::endl;
  
  // Attach the site indicators to the rest of the design matrix
  Eigen::MatrixXd XC(X.rows(), X.cols() + site_indicators.cols());
  XC.leftCols(X.cols())                = X;
  XC.rightCols(site_indicators.cols()) = site_indicators;
  
  // Attach the random effects to the fixed + site effects design matrix
  Eigen::MatrixXd XCXtilde(XC.rows(), XC.cols() + Xrand.cols());
  XCXtilde.leftCols(XC.cols())     = XC;
  XCXtilde.rightCols(Xrand.cols()) = Xrand;
  
  // XtX for posterior precision
  Eigen::MatrixXd XtX = XCXtilde.transpose() * XCXtilde;
  
  // Zero out E_Si_sq term. This gets accumulated over in the loop below
  for (int i = 0; i < NJ; i++){
    trace_E_Si_sq(i) = 0.0;
  }
  
  // Zero out the term storing the sums of the bi_squared terms
  sum_E_b_site_sq(0) = 0.0;
  
  // This version is needed for the cluster membership updates
  for (int j = 0; j < E_tr_BiBit_qv.cols(); j++){
    for (int i = 0; i < E_tr_BiBit_qv.rows(); i++){
      E_tr_BiBit_qv(i, j) = 0.0;
    }
  }
  
  // Storage for intermediate quantities
  double s_prior_mean = 0.0;
  double s_prior_precision = 0.0;
  Eigen::MatrixXd prior_precision     = Eigen::MatrixXd::Zero(P + 1 + nSite + N*Prand, P + 1 + nSite + N*Prand);
  Eigen::MatrixXd posterior_precision = Eigen::MatrixXd::Zero(P + 1 + nSite + N*Prand, P + 1 + nSite + N*Prand);
  Eigen::MatrixXd posterior_variance  = Eigen::MatrixXd::Zero(P + 1 + nSite + N*Prand, P + 1 + nSite + N*Prand);
  Eigen::VectorXd posterior_mean      = Eigen::VectorXd::Zero(P + 1 + nSite + N*Prand);
  
  Eigen::MatrixXd random_effect_var  = Eigen::MatrixXd::Zero(N*Prand, N*Prand);
  
  // Intermediate quantity storing the expected precision of the random effects for 
  // the specific network/location combination
  Eigen::MatrixXd random_eff_E_precision = Eigen::MatrixXd::Zero(Prand, Prand);
  Eigen::MatrixXd random_eff_E_precision_inv = Eigen::MatrixXd::Zero(Prand, Prand);
  Eigen::MatrixXd E_BiBit                = Eigen::MatrixXd::Zero(Prand, Prand);
  
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      
      // Get the S prior mean and variance based on the DPM cluster weights and parameters
      s_prior_mean      = 0.0;
      s_prior_precision = 0.0;
      for (int h = 0; h < H; h++){
        s_prior_mean      += pi_Zqv_eq_h(h, iq*V + iv) * E_mu_h(h);
        s_prior_precision += pi_Zqv_eq_h(h, iq*V + iv) * E_gamma_h_sq_inv(h);
      }
      
      // if (iq == 0){
      //   if (iv == 0){
      //     std::cout << "S prior mean and precision: " << std::endl;
      //     std::cout << s_prior_mean << "   and  " << s_prior_precision << std::endl;
      //   }
      // }
      
      
      // Get the B_{qv} prior precision based on the DPM cluster weights and parameters
      random_eff_E_precision *= 0.0;
      for (int h_re = 0; h_re < H_re; h_re++){
        random_eff_E_precision += pi_Ztildeqv_eq_h(h_re, iq*V + iv) * 
          E_Rqv_h_inv.block(0, h_re*Prand, Prand, Prand);
      }
      random_eff_E_precision_inv = random_eff_E_precision.inverse();
      
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
      for (int i = 0; i < N; i++){
        prior_precision.block(P + 1 + nSite + Prand*i,
                              P + 1 + nSite + Prand*i,
                              Prand, Prand) = random_eff_E_precision;
        random_effect_var.block(Prand*i, Prand*i, Prand, Prand) = random_eff_E_precision_inv;
      }
 
      // Evaluate the variational distn precision
      posterior_precision = E_sigma_sq_inv * XtX + prior_precision;

      // **Block Inversion Using Schur Complement**
      size_t block_size = P + 1 + nSite;
      size_t rand_size = N * Prand;
      
      Eigen::MatrixXd A11 = posterior_precision.block(0, 0, block_size, block_size);
      Eigen::MatrixXd A12 = posterior_precision.block(0, block_size, block_size, rand_size);
      Eigen::MatrixXd A21 = posterior_precision.block(block_size, 0, rand_size, block_size);
      Eigen::MatrixXd A22 = posterior_precision.block(block_size, block_size, rand_size, rand_size);
      
      // Compute A11_inv
      Eigen::MatrixXd A11_inv = A11.llt().solve(Eigen::MatrixXd::Identity(A11.rows(), A11.cols()));
      
      // Compute block diagonal A22^{-1}
      Eigen::MatrixXd A22_inv = Eigen::MatrixXd::Zero(rand_size, rand_size);
      for (int i = 0; i < N; i++) {
        Eigen::MatrixXd block_inv = A22.block(i * Prand, i * Prand, Prand, Prand).inverse();
        A22_inv.block(i * Prand, i * Prand, Prand, Prand) = block_inv;
      }
      
      // Compute Schur Complement S = A11 - A12 * A22^{-1} * A21
      Eigen::MatrixXd A12_A22_inv = A12 * A22_inv;
      Eigen::MatrixXd S = A11 - A12_A22_inv * A21;
      
      // Compute S^{-1}
      Eigen::MatrixXd S_inv = S.llt().solve(Eigen::MatrixXd::Identity(S.rows(), S.cols()));
      
      // Compute A^{-1} using block inversion
      //posterior_variance.setZero(P + N * Prand, P + N * Prand);
      posterior_variance.block(0, 0, block_size, block_size) = S_inv;
      posterior_variance.block(0, block_size, block_size, rand_size) = -S_inv * A12_A22_inv;
      posterior_variance.block(block_size, 0, rand_size, block_size) = -A22_inv * A21 * S_inv;
      posterior_variance.block(block_size, block_size, rand_size, rand_size) = A22_inv + A22_inv * A21 * S_inv * A12_A22_inv;
      
      
      
      // Invert to obtain variational distn variance
      //posterior_variance = posterior_precision.inverse();
      
      // if (iq == 0){
      //   if (iv == 0){
      //     std::cout << "The posterior variance for this element is: " << std::endl;
      //     std::cout << posterior_variance.block(0, 0, P + 1 + nSite + Prand*2, P + 1 + nSite + Prand*2) << std::endl;
      //   }
      // }
      
      // Get the variational mean
      posterior_mean = XCXtilde.transpose() * Y.col(iq*V + iv) * E_sigma_sq_inv;
      
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
      // We don't store Biqv, but we do need some summary statistics based on it
      for (int i = 0; i < N; i++){
        E_BiBit = posterior_variance.block(P + 1 + nSite + i*Prand,
                                           P + 1 + nSite + i*Prand,
                                           Prand, Prand) +
                                             posterior_mean.segment(P + 1 + nSite + i*Prand, Prand) *
                                             posterior_mean.segment(P + 1 + nSite + i*Prand, Prand).transpose();
        E_tr_BiBit_qv.block(0, iq*V*Prand + iv*Prand, Prand, Prand) += E_BiBit;
      }

      // Update subject-specific maps (linear combination of other variables)
      // TODO precalc posterior_variance + posterior_mean*posterior_mean.transpose()
      Eigen::MatrixXd middle_term = (posterior_variance + posterior_mean*posterior_mean.transpose());
      for (int ij = 0; ij < NJ; ij++){
        E_Si(ij, iq*V + iv) = (XCXtilde.row(ij) * posterior_mean)(0);
        trace_E_Si_sq(ij)  += (XCXtilde.row(ij) *
          middle_term *
          XCXtilde.row(ij).transpose())(0);
      }
      
    } // v
  } // q
  
}









// [[Rcpp::export]]
void longitudinal_covadj_spatial_map_updates(Eigen::Map<Eigen::MatrixXd> & E_S,
                                                     Eigen::Map<Eigen::MatrixXd> & E_S_sq,
                                                     Eigen::Map<Eigen::MatrixXd> & E_Beta,
                                                     Eigen::Map<Eigen::MatrixXd> & E_Beta_sq,
                                                     Eigen::Map<Eigen::MatrixXd> & E_Si,
                                                     Eigen::Map<Eigen::MatrixXd> & trace_E_Si_sq,
                                                     Eigen::Map<Eigen::MatrixXd> & E_tr_BiBit_qv,
                                                     const Eigen::Map<Eigen::MatrixXd> & Y,
                                                     const Eigen::Map<Eigen::MatrixXd> & X,
                                                     const Eigen::Map<Eigen::MatrixXd> & Xrand,
                                                     const Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                                                     const Eigen::Map<Eigen::MatrixXd> & pi_Ztildeqv_eq_h,
                                                     const Eigen::Map<Eigen::VectorXd> & E_gamma_h_sq_inv,
                                                     const Eigen::Map<Eigen::MatrixXd> & E_Rqv_h_inv,
                                                     const Eigen::Map<Eigen::VectorXd> & E_mu_h,
                                                     const Eigen::Map<Eigen::VectorXd> & Ji,
                                                     const Eigen::Map<Eigen::MatrixXd> & E_lambda_sq,
                                                     double E_sigma_sq_inv,
                                                     double E_tau_sq_inv) {
  
  // Dimensions
  size_t V     = E_S.rows();
  size_t Q     = E_S.cols();
  size_t P     = E_Beta.rows();
  size_t Prand = E_Rqv_h_inv.rows();
  size_t H     = E_mu_h.size();
  size_t NJ    = Y.rows();
  size_t N     = Xrand.cols() / Prand;
  size_t H_re  = E_Rqv_h_inv.cols() / Prand;
  
  //std::cout << "N = " << N << std::endl;
  
  // Attach the random effects to the fixed effects design matrix
  Eigen::MatrixXd XCXtilde(X.rows(), X.cols() + Xrand.cols());
  XCXtilde.leftCols(X.cols())     = X;
  XCXtilde.rightCols(Xrand.cols()) = Xrand;
  
  // XtX for posterior precision
  Eigen::MatrixXd XtX = XCXtilde.transpose() * XCXtilde;
  
  // Zero out E_Si_sq term. This gets accumulated over in the loop below
  for (int i = 0; i < NJ; i++){
    trace_E_Si_sq(i) = 0.0;
  }
  
  // This version is needed for the cluster membership updates
  for (int j = 0; j < E_tr_BiBit_qv.cols(); j++){
    for (int i = 0; i < E_tr_BiBit_qv.rows(); i++){
      E_tr_BiBit_qv(i, j) = 0.0;
    }
  }
  
  // Storage for intermediate quantities
  double s_prior_mean = 0.0;
  double s_prior_precision = 0.0;
  Eigen::MatrixXd prior_precision     = Eigen::MatrixXd::Zero(P + 1  + N*Prand, P + 1  + N*Prand);
  Eigen::MatrixXd posterior_precision = Eigen::MatrixXd::Zero(P + 1  + N*Prand, P + 1  + N*Prand);
  Eigen::MatrixXd posterior_variance  = Eigen::MatrixXd::Zero(P + 1  + N*Prand, P + 1  + N*Prand);
  Eigen::VectorXd posterior_mean      = Eigen::VectorXd::Zero(P + 1  + N*Prand);
  
  // Intermediate quantity storing the expected precision of the random effects for 
  // the specific network/location combination
  Eigen::MatrixXd random_eff_E_precision = Eigen::MatrixXd::Zero(Prand, Prand);
  Eigen::MatrixXd E_BiBit                = Eigen::MatrixXd::Zero(Prand, Prand);
  
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      
      // Get the S prior mean and variance based on the DPM cluster weights and parameters
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
      
      // Update S, S^2, and Beta expectations
      E_S(iv, iq)              = posterior_mean(0);
      E_S_sq(iv, iq)           = pow(posterior_mean(0), 2) + posterior_variance(0, 0);
      for (int p = 0; p < P; p++){
        E_Beta(p, iq*V + iv)    = posterior_mean(p + 1);
        E_Beta_sq(p, iq*V + iv) = pow(posterior_mean(p + 1), 2) + 
          posterior_variance(p + 1, p + 1);
      }
      // We don't store Biqv, but we do need some summary statistics based on it
      for (int i = 0; i < N; i++){
        E_BiBit = posterior_variance.block(P + 1  + i*Prand,
                                           P + 1  + i*Prand,
                                           Prand, Prand) + 
                                             posterior_mean.segment(P + 1  + i*Prand, Prand) *
                                             posterior_mean.segment(P + 1  + i*Prand, Prand).transpose();
        E_tr_BiBit_qv.block(0, iq*V*Prand + iv*Prand, Prand, Prand) += E_BiBit;
      }
      
      // Update subject-specific maps (linear combination of other variables)
      // TODO precalc posterior_variance + posterior_mean*posterior_mean.transpose()
      Eigen::MatrixXd middle_term = (posterior_variance + posterior_mean*posterior_mean.transpose());
      for (int ij = 0; ij < NJ; ij++){
        E_Si(ij, iq*V + iv) = (XCXtilde.row(ij) * posterior_mean)(0);
        trace_E_Si_sq(ij)  += (XCXtilde.row(ij) * 
          middle_term *
          XCXtilde.row(ij).transpose())(0);
      }
      
    } // v
  } // q
  
}
