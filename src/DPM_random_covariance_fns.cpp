// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// E_tr_BiBit_qv is Prand x QV. is not actually trace, is involved in that term
// [[Rcpp::export]]
void calibrate_Rqv_cluster_probabilities(Eigen::Map<Eigen::MatrixXd> & pi_Ztildeqv_eq_h,
                                         const Eigen::Map<Eigen::MatrixXd> & E_tr_BiBit_qv,
                                         const Eigen::Map<Eigen::MatrixXd> & E_log_det_Rqv_inv,
                                         const Eigen::Map<Eigen::MatrixXd> & E_Rqv_inv,
                                         const Eigen::Map<Eigen::VectorXd> & E_Rqv_log_nu_h,
                                         const Eigen::Map<Eigen::VectorXd> & E_Rqv_log_1_minus_nu_h,
                                         double E_sigma_sq_inv,
                                         int V,
                                         int N){
  
  // Dimensions
  size_t Prand    = E_tr_BiBit_qv.rows();
  size_t Hmax_Rqv = E_Rqv_log_nu_h.size();
  size_t Q        = E_tr_BiBit_qv.cols() / Prand / V;
  
  // Precalculate the terms that don't depend on the voxel- and IC-specific info
  double cumulative_sum_term = 0; // accumulates log 1-nu terms
  Eigen::VectorXd leading_terms = Eigen::VectorXd::Zero(Hmax_Rqv);
  for (int h = 0; h < Hmax_Rqv; h++){
    leading_terms(h)     = E_Rqv_log_nu_h(h) + cumulative_sum_term + 0.5 * N * E_log_det_Rqv_inv(h);
    cumulative_sum_term += E_Rqv_log_1_minus_nu_h(h);
  }
  
  // Now calculate the probability of belonging to each cluster for each location/component
  Eigen::VectorXd unnorm_log_probs = Eigen::VectorXd::Zero(Hmax_Rqv);
  Eigen::VectorXd normed_log_probs = Eigen::VectorXd::Zero(Hmax_Rqv);
  double max_log_prob = 0.0;
  double normalization_term = 0.0;
  double sum_unnormalized_probs = 0.0;
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      
      unnorm_log_probs *= 0.0;
      
      // Loop over clusters, probabilities are proportional to these terms
      for (int h = 0; h < Hmax_Rqv; h++){
        unnorm_log_probs(h) = leading_terms(h) - 
          0.5 * (E_tr_BiBit_qv.block(0, Prand*V*iq + Prand*iv, Prand, Prand) * 
          E_Rqv_inv.block(0, Prand*h, Prand, Prand)).trace();
      }
      
      // if (iq == 0){
      //   if (iv == 0){
      //     std::cout << "leading terms:" << std::endl;
      //     std::cout << leading_terms.transpose() << std::endl;
      //     std::cout << "unnormalized probabilities:" << std::endl;
      //     std::cout << unnorm_log_probs.transpose() << std::endl;
      //   }
      // }
      
      // Log sum exp trick to normalize
      max_log_prob = unnorm_log_probs.maxCoeff();
      for (int h = 0; h < Hmax_Rqv; h++){
        normed_log_probs(h) = exp(unnorm_log_probs(h) - max_log_prob);
      }
      normalization_term = max_log_prob + log(normed_log_probs.sum());
      
      // Store the updated probabilities
      for (int h = 0; h < Hmax_Rqv; h++){
        pi_Ztildeqv_eq_h(h, iq*V + iv) = exp(unnorm_log_probs(h) - normalization_term);
      }
      
      // Check and fix underflow
      for (int h = 0; h < Hmax_Rqv; h++){
        if (pi_Ztildeqv_eq_h(h, iq*V + iv) < 0.0000000001){
          pi_Ztildeqv_eq_h(h, iq*V + iv) = 0.0000000001;
        }
      }
      normalization_term = pi_Ztildeqv_eq_h.col(iq*V + iv).sum();
      for (int h = 0; h < Hmax_Rqv; h++){
        pi_Ztildeqv_eq_h(h, iq*V + iv) /= normalization_term;
      }
      
      
    } // v
  } // q
  
}




// pi is h by QV
// [[Rcpp::export]]
void calibrate_Rqv_cluster_parameters(Eigen::Map<Eigen::VectorXd> & VarHyper_Rqv_df,
                                      Eigen::Map<Eigen::MatrixXd> & VarHyper_Rqv_scale,
                                      const Eigen::Map<Eigen::MatrixXd> & pi_Ztildeqv_eq_h,
                                      const Eigen::Map<Eigen::MatrixXd> & E_tr_BiBit_qv,
                                      const Eigen::Map<Eigen::MatrixXd> & base_measure_scale,
                                      double E_sigma_sq_inv,
                                      double base_measure_df,
                                      int V,
                                      int N){
  
  // Dimensions
  size_t Prand    = E_tr_BiBit_qv.rows();
  size_t Hmax_Rqv = pi_Ztildeqv_eq_h.rows();
  size_t Q        = E_tr_BiBit_qv.cols() / Prand / V;
  
  // Intermediate quantities
  Eigen::VectorXd sum_pi_h           = Eigen::VectorXd::Zero(Hmax_Rqv);
  Eigen::MatrixXd sum_pi_h_times_BBt = Eigen::MatrixXd::Zero(Hmax_Rqv, Hmax_Rqv * Prand);

  // Calculate summary statistics for each cluster aggregated over all spatial
  // locations
  for (int iq = 0; iq < Q; iq++){
    for (int iv = 0; iv < V; iv++){
      for (int h = 0; h < Hmax_Rqv; h++){
        
        sum_pi_h(h)           += pi_Ztildeqv_eq_h(h, iq*V + iv);
        sum_pi_h_times_BBt.block(0, Prand*h, Prand, Prand) += 
          pi_Ztildeqv_eq_h(h, iq*V + iv) * 
          E_tr_BiBit_qv.block(0, Prand*V*iq + Prand*iv, Prand, Prand);

      } // h
    } // v
  } // q
  
  
  // Update the variational free parameters
  for (int h = 0; h < Hmax_Rqv; h++){
    
    VarHyper_Rqv_df(h) = base_measure_df + N * sum_pi_h(h);
    VarHyper_Rqv_scale.block(0, h*Prand, Prand, Prand) = base_measure_scale + 
      sum_pi_h_times_BBt.block(0, Prand*h, Prand, Prand);
    
  }
  
  
}

