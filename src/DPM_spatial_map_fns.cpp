// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// pi_Zqv_eq_h is H x QV
// [[Rcpp::export]]
void calibrate_spatial_map_cluster_probabilities(Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                                                 const Eigen::Map<Eigen::VectorXd> & E_log_nu_h,
                                                 const Eigen::Map<Eigen::VectorXd> & E_log_1_minus_nu_h,
                                                 const Eigen::Map<Eigen::VectorXd> & E_log_gamma_h_sq,
                                                 const Eigen::Map<Eigen::VectorXd> & E_gamma_h_sq_inv,
                                                 const Eigen::Map<Eigen::VectorXd> & E_mu_h,
                                                 const Eigen::Map<Eigen::VectorXd> & E_mu_h_sq,
                                                 const Eigen::Map<Eigen::MatrixXd> & E_S,
                                                 const Eigen::Map<Eigen::MatrixXd> & E_S_sq,
                                                 double E_sigma_sq_inv){
  
  // Dimensions
  size_t V = E_S.rows();
  size_t Q = E_S.cols();
  size_t H = E_mu_h.size();
  
  // Precalculate the terms that don't depend on the voxel- and IC-specific info
  double cumulative_sum_term = 0; // accumulates log 1-nu terms
  Eigen::VectorXd leading_terms = Eigen::VectorXd::Zero(H);
  for (int h = 0; h < H; h++){
    leading_terms(h) = E_log_nu_h(h) + cumulative_sum_term - E_log_gamma_h_sq(h);
    cumulative_sum_term += E_log_1_minus_nu_h(h);
  }
  
  // Now calculate the probability of belonging to each cluster for each location/component
  Eigen::VectorXd unnorm_log_probs = Eigen::VectorXd::Zero(H);
  Eigen::VectorXd normed_log_probs = Eigen::VectorXd::Zero(H);
  double max_log_prob = 0.0;
  double normalization_term = 0.0;
  double sum_unnormalized_probs = 0.0;
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      
      unnorm_log_probs *= 0.0;
      
      // Loop over clusters, probabilities are proportional to these terms
      for (int h = 0; h < H; h++){
        unnorm_log_probs(h) = leading_terms(h) - 
          0.5 * E_sigma_sq_inv * E_gamma_h_sq_inv(h) * 
          (E_S_sq(iv, iq) - 2.0 * E_S(iv, iq) * E_mu_h(h) + E_mu_h_sq(h));
      }
      
      // Log sum exp trick to normalize
      max_log_prob = unnorm_log_probs.maxCoeff();
      for (int h = 0; h < H; h++){
        normed_log_probs(h) = exp(unnorm_log_probs(h) - max_log_prob);
      }
      normalization_term = max_log_prob + log(normed_log_probs.sum());
      
      // Store the updated probabilities
      for (int h = 0; h < H; h++){
        pi_Zqv_eq_h(h, iq*V + iv) = exp(unnorm_log_probs(h) - normalization_term);
      }
      
      // Check and fix underflow
      for (int h = 0; h < H; h++){
        if (pi_Zqv_eq_h(h, iq*V + iv) < 0.0000000001){
          pi_Zqv_eq_h(h, iq*V + iv) = 0.0000000001;
        }
      }
      normalization_term = pi_Zqv_eq_h.col(iq*V + iv).sum();
      for (int h = 0; h < H; h++){
        pi_Zqv_eq_h(h, iq*V + iv) /= normalization_term;
      }

    } // v
  } // q
  
}




// pi is h by QV
// [[Rcpp::export]]
void calibrate_spatial_map_cluster_parameters(Eigen::Map<Eigen::VectorXd> & VarHyper_M_h,
                                              Eigen::Map<Eigen::VectorXd> & VarHyper_L_h,
                                              Eigen::Map<Eigen::VectorXd> & VarHyper_A_h,
                                              Eigen::Map<Eigen::VectorXd> & VarHyper_B_h,
                                              const Eigen::Map<Eigen::MatrixXd> & pi_Zqv_eq_h,
                                              const Eigen::Map<Eigen::MatrixXd> & E_S,
                                              const Eigen::Map<Eigen::MatrixXd> & E_S_sq, 
                                              double E_sigma_sq_inv, 
                                              double base_measure_lambda,
                                              double base_measure_mu, 
                                              double base_measure_alpha, 
                                              double base_measure_beta){
  
  // Number of clusters in truncated DPM
  size_t H  = VarHyper_M_h.size();
  size_t V  = E_S.rows();
  size_t Q  = E_S.cols();
  
  // Intermediate quantities
  Eigen::VectorXd sum_pi_h           = Eigen::VectorXd::Zero(H);
  Eigen::VectorXd sum_pi_h_times_S   = Eigen::VectorXd::Zero(H);
  Eigen::VectorXd sum_pi_h_times_Ssq = Eigen::VectorXd::Zero(H);
  
  // Calculate summary statistics for each cluster aggregated over all spatial
  // locations
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      for (int h = 0; h < H; h++){
        
        sum_pi_h(h)           += pi_Zqv_eq_h(h, iq*V + iv);
        sum_pi_h_times_S(h)   += pi_Zqv_eq_h(h, iq*V + iv) * E_S(iv, iq);
        sum_pi_h_times_Ssq(h) += pi_Zqv_eq_h(h, iq*V + iv) * E_S_sq(iv, iq);
        
      } // h
    } // v
  } // q
  
  
  // Update the variational free parameters
  for (int h = 0; h < H; h++){
    
    // Mean
    VarHyper_M_h(h) = (E_sigma_sq_inv*sum_pi_h_times_S(h) + 
      base_measure_lambda*base_measure_mu) / 
      (base_measure_lambda + E_sigma_sq_inv*sum_pi_h(h));
    
    VarHyper_L_h(h) = base_measure_lambda + E_sigma_sq_inv * sum_pi_h(h);
    VarHyper_A_h(h) = base_measure_alpha  + E_sigma_sq_inv * sum_pi_h(h)/2.0;
    
    // VarHyper_B_h(h) = 0.5 * (E_sigma_sq_inv*sum_pi_h(h) * base_measure_lambda) /
    //   VarHyper_L_h(h) * pow(sum_pi_h_times_S(h) / sum_pi_h(h) - base_measure_mu, 2) +
    //   base_measure_beta + E_sigma_sq_inv * sum_pi_h_times_Ssq(h)/2 -
    //   (1/VarHyper_L_h(h)) * E_sigma_sq_inv * (pow(sum_pi_h_times_S(h), 2) /sum_pi_h(h))/2 -
    //   (1/VarHyper_L_h(h)) * E_sigma_sq_inv * pow(sum_pi_h_times_S(h), 2) /2;
    
    // std::cout << "explicit check:" << (sum_pi_h_times_Ssq(h)/sum_pi_h(h)  - pow(sum_pi_h_times_S(h)/sum_pi_h(h),2)) << std::endl;
    
    VarHyper_B_h(h) = base_measure_beta + 
      0.5 * E_sigma_sq_inv * 
      (
          sum_pi_h(h) *(sum_pi_h_times_Ssq(h)/sum_pi_h(h)  - pow(sum_pi_h_times_S(h)/sum_pi_h(h),2)) +
            base_measure_lambda * sum_pi_h(h) / (base_measure_lambda + sum_pi_h(h)) * pow(sum_pi_h_times_S(h)/sum_pi_h(h), 2)
      );

    // VarHyper_B_h(h) = 0.5 * (E_sigma_sq_inv*sum_pi_h(h) * base_measure_lambda) / 
    //   VarHyper_L_h(h) * pow(sum_pi_h_times_S(h) / sum_pi_h(h) - base_measure_mu, 2) + 
    //   base_measure_beta + E_sigma_sq_inv * sum_pi_h_times_Ssq(h)/2 -
    //   (1/VarHyper_L_h(h)) * E_sigma_sq_inv * (pow(sum_pi_h_times_S(h), 2) /sum_pi_h(h))/2 - 
    //   (base_measure_lambda * E_sigma_sq_inv * sum_pi_h(h)/VarHyper_L_h(h)) * E_sigma_sq_inv * sum_pi_h_times_S(h) /2;
    // 
    // current_calc    <- 0.5 * (E_sigma_sq_inv*sum_pi_h * base_measure_lambda) / 
    //   VarHyper_L_h[h] * (sum_pi_s_qv/sum_pi_h - base_measure_mu)^2 + 
    //   base_measure_beta 
    // current_calc    <- current_calc + E_sigma_sq_inv*sum(pi_ssq_qv)/2 - 
    //   (1/VarHyper_L_h[h]) * (sum_pi_s_qv^2 /sum_pi_h)/2 - 
    //   (1/VarHyper_L_h[h]) * E_sigma_sq_inv * sum_pi_s_qv^2 /2

  }

  
}

