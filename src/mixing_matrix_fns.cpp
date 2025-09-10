// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// Y is N x QV
// Base version
// [[Rcpp::export]]
void update_all_mixing_matrices(Eigen::Map<Eigen::MatrixXd> & A,
                                const Eigen::Map<Eigen::MatrixXd> & Y_mixed_QxVN,
                                const Eigen::Map<Eigen::MatrixXd> & E_Si,
                                double E_sigma_sq_inv){
  
  
  size_t Q = A.rows();
  size_t N = A.cols() / Q;
  size_t V = Y_mixed_QxVN.cols() / N;

  // Initialize variables
  Eigen::MatrixXd mean_direction = Eigen::MatrixXd::Zero(Q, Q);
  Eigen::MatrixXd E_Si_VxQ       = Eigen::MatrixXd::Zero(V, Q);

  for (int i = 0; i < N; i++){
    
    for (int iq = 0; iq < Q; iq++){
      E_Si_VxQ.col(iq) = E_Si.block(i, iq*V, 1, V).transpose();
    }
    
    // Calculate the direction parameter
    mean_direction = E_sigma_sq_inv * (E_Si_VxQ.transpose() * Y_mixed_QxVN.block(0, i*V, Q, V).transpose()).transpose();
    
    // Perform a singular value decomposition
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mean_direction, Eigen::ComputeFullU | Eigen::ComputeFullV);

    // Get the mode (should also be expected value for rows = cols case)
    // and update stored version of mixing matrix
    A.block(0, i*Q, Q, Q) = svd.matrixU() * svd.matrixV().transpose();
    
  } // end of loop over scans
  
}

// Y_unmixed is N x QV
// [[Rcpp::export]]
void apply_unmixing_matrices(Eigen::Map<Eigen::MatrixXd> & Y_unmixed,
                             const Eigen::Map<Eigen::MatrixXd> & Y_mixed_QxNV,
                             const Eigen::Map<Eigen::MatrixXd> & A){
  
  size_t N_time_courses = Y_unmixed.rows();
  size_t Q = A.rows();
  size_t V = Y_unmixed.cols() / Q;
  
  // Storage space
  Eigen::MatrixXd VxQ_time_course = Eigen::MatrixXd::Zero(V, Q);
  
  for (int i = 0; i < N_time_courses; i++){
    
    // Turn time series i into a V x Q matrix for unmixing
    VxQ_time_course = Y_mixed_QxNV.block(0, V*i, Q, V).transpose();
    
    // Apply the unmixing
    // Note - normally A' Y, when Y is Q x V, so in this case Y' A
    VxQ_time_course = VxQ_time_course * A.block(0, Q*i, Q, Q);
    
    // Store in Y_unmixed
    for (int iq = 0; iq < Q; iq++){
      Y_unmixed.block(i, iq*V, 1, V) = VxQ_time_course.col(iq).transpose();
    }
    
  } // end of loop over subjects
  
}