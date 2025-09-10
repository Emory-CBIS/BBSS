// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// This version is for the case where there are no covariate effects. See covadj 
// version if covariates are included
// [[Rcpp::export]]
double update_error_variance_variational_scale(const Eigen::Map<Eigen::MatrixXd> & Y_unmixed,
                             const Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                             const Eigen::Map<Eigen::MatrixXd> & E_S,
                             const Eigen::Map<Eigen::MatrixXd> & E_Si,
                             const Eigen::Map<Eigen::MatrixXd> & E_S_sq,
                             const Eigen::Map<Eigen::VectorXd> & E_mu_h,
                             const Eigen::Map<Eigen::VectorXd> & E_mu_h_sq,
                             const Eigen::Map<Eigen::VectorXd> & E_gamma_h_sq_inv,
                             double sum_trace_YYt,
                             double sum_trace_E_Si_sq){
  
  size_t N = Y_unmixed.rows();
  size_t V = E_S.rows();
  size_t Q = E_S.cols();
  size_t H = E_mu_h.size();
  
  double VariHyper_B_sigma_sq     = 0.0;
  double sum_pi_h_times_mu        = 0.0;
  double sum_pi_h_times_mu_sq     = 0.0;
  double sum_pi_h_times_gamma_inv = 0.0;
  
  // Term 1: Contribution from errors at the (unmixed) time series level (Y)
  double sse_term = 0.0;
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      for (int i = 0; i < N; i++){
        sse_term -= Y_unmixed(i, iq*V + iv) * E_Si(i, iq*V + iv);
      }
    }
  }
  sse_term += sum_trace_YYt / 2.0 + sum_trace_E_Si_sq / 2.0;
  
  // Term 2: Contribution from S
  double s_term = 0.0;
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      
      // Reset terms
      sum_pi_h_times_mu        = 0.0;
      sum_pi_h_times_mu_sq     = 0.0;
      sum_pi_h_times_gamma_inv = 0.0;
      
      for (int h = 0; h < H; h++){
        sum_pi_h_times_mu        += pi_Zqv_eq_h(h, iq*V + iv) * E_mu_h(h);
        sum_pi_h_times_mu_sq     += pi_Zqv_eq_h(h, iq*V + iv) * E_mu_h_sq(h);
        sum_pi_h_times_gamma_inv += pi_Zqv_eq_h(h, iq*V + iv) * E_gamma_h_sq_inv(h);
      } // h
      
      s_term += 0.5 * (
        sum_pi_h_times_gamma_inv * (E_S_sq(iv, iq) - 2 * E_S(iv, iq) * sum_pi_h_times_mu + sum_pi_h_times_mu_sq)
      );
    }
  }
  
  VariHyper_B_sigma_sq = sse_term + s_term;
  
  return(VariHyper_B_sigma_sq);
  
}






// [[Rcpp::export]]
double update_error_variance_variational_scale_covadj(const Eigen::Map<Eigen::MatrixXd> & Y_unmixed,
                                               const Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                                               const Eigen::Map<Eigen::MatrixXd> & E_S,
                                               const Eigen::Map<Eigen::MatrixXd> & E_Si,
                                               const Eigen::Map<Eigen::MatrixXd> & E_S_sq,
                                               const Eigen::Map<Eigen::VectorXd> & E_mu_h,
                                               const Eigen::Map<Eigen::VectorXd> & E_mu_h_sq,
                                               const Eigen::Map<Eigen::VectorXd> & E_gamma_h_sq_inv,
                                               const Eigen::Map<Eigen::MatrixXd> & E_Beta_sq,
                                               const Eigen::Map<Eigen::MatrixXd> & E_lambda_sq,
                                               double sum_trace_YYt,
                                               double sum_trace_E_Si_sq,
                                               double E_tau_sq_inv){
  
  size_t N = Y_unmixed.rows();
  size_t V = E_S.rows();
  size_t Q = E_S.cols();
  size_t H = E_mu_h.size();
  size_t P = E_Beta_sq.rows();
  
  double VariHyper_B_sigma_sq     = 0.0;
  double sum_pi_h_times_mu        = 0.0;
  double sum_pi_h_times_mu_sq     = 0.0;
  double sum_pi_h_times_gamma_inv = 0.0;
  
  // Term 1: Contribution from errors at the (unmixed) time series level (Y)
  double sse_term = 0.0;
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      for (int i = 0; i < N; i++){
        sse_term -= Y_unmixed(i, iq*V + iv) * E_Si(i, iq*V + iv);
      }
    }
  }
  sse_term += sum_trace_YYt / 2.0 + sum_trace_E_Si_sq / 2.0;
  
  // Term 2: Contribution from S
  double s_term = 0.0;
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      
      // Reset terms
      sum_pi_h_times_mu        = 0.0;
      sum_pi_h_times_mu_sq     = 0.0;
      sum_pi_h_times_gamma_inv = 0.0;
      
      for (int h = 0; h < H; h++){
        sum_pi_h_times_mu        += pi_Zqv_eq_h(h, iq*V + iv) * E_mu_h(h);
        sum_pi_h_times_mu_sq     += pi_Zqv_eq_h(h, iq*V + iv) * E_mu_h_sq(h);
        sum_pi_h_times_gamma_inv += pi_Zqv_eq_h(h, iq*V + iv) * E_gamma_h_sq_inv(h);
      } // h
      
      s_term += 0.5 * (
        sum_pi_h_times_gamma_inv * (E_S_sq(iv, iq) - 2 * E_S(iv, iq) * sum_pi_h_times_mu + sum_pi_h_times_mu_sq)
      );
    }
  }
  
  
  // Term 3: Contribution from beta
  double beta_term = 0.0;
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      for (int ip = 0; ip < P; ip++){
        beta_term += 0.5 * E_tau_sq_inv * E_lambda_sq(ip, iq*V + iv) * E_Beta_sq(ip, iq*V + iv);
      }
    }
  }
  
  
  VariHyper_B_sigma_sq = sse_term + s_term + beta_term;
  
  return(VariHyper_B_sigma_sq);
  
}
