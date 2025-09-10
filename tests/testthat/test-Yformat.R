test_that("Longitudinal List of Lists Works", {
  
  true_Ji <- c(2, 3)
  subj1_scan1 <- matrix(c(1,2,3,4,5,6),     nrow = 2)
  subj1_scan2 <- matrix(c(1,2,3,4,5,6) + 1, nrow = 2)
  subj2_scan1 <- matrix(c(1,2,3,4,5,6) + 2, nrow = 2)
  subj2_scan2 <- matrix(c(1,2,3,4,5,6) + 3, nrow = 2)
  subj2_scan3 <- matrix(c(1,2,3,4,5,6) + 4, nrow = 2)
  subj1 <- list()
  subj1[[1]] <- subj1_scan1
  subj1[[2]] <- subj1_scan2
  subj2 <- list()
  subj2[[1]] <- subj2_scan1
  subj2[[2]] <- subj2_scan2
  subj2[[3]] <- subj2_scan3
  
  Y <- list()
  Y[[1]] <- subj1
  Y[[2]] <- subj2
  
  id_variable <- c("1", "1", "2", "2", "2")
  
  # Run the function
  result <- format_time_series_input(Y, longit_adjusted = TRUE, debug = FALSE,
                                     id_variable = id_variable)
  
  # Checking the result is correct
  check_subj1_scan1 <- all(result$Y_mixed_QxVNstar[1:2, 1:3]   == subj1_scan1)
  check_subj1_scan2 <- all(result$Y_mixed_QxVNstar[1:2, 4:6]   == subj1_scan2)
  check_subj2_scan1 <- all(result$Y_mixed_QxVNstar[1:2, 7:9]   == subj2_scan1)
  check_subj2_scan2 <- all(result$Y_mixed_QxVNstar[1:2, 10:12] == subj2_scan2)
  check_subj2_scan3 <- all(result$Y_mixed_QxVNstar[1:2, 13:15] == subj2_scan3)
  
  testthat::expect_true(all(check_subj1_scan1, check_subj1_scan2, check_subj2_scan1,
                            check_subj2_scan2, check_subj2_scan3))
  
})






test_that("Longitudinal List of Arrays Works", {
  
  true_Ji <- c(2, 3)
  subj1_scan1 <- matrix(c(1,2,3,4,5,6),     nrow = 2)
  subj1_scan2 <- matrix(c(1,2,3,4,5,6) + 1, nrow = 2)
  subj2_scan1 <- matrix(c(1,2,3,4,5,6) + 2, nrow = 2)
  subj2_scan2 <- matrix(c(1,2,3,4,5,6) + 3, nrow = 2)
  subj2_scan3 <- matrix(c(1,2,3,4,5,6) + 4, nrow = 2)
  subj1 <- array(0, dim = c(2, 3, 2))
  subj1[,,1] <- subj1_scan1
  subj1[,,2] <- subj1_scan2
  subj2 <- array(0, dim = c(2, 3, 3))
  subj2[,,1] <- subj2_scan1
  subj2[,,2] <- subj2_scan2
  subj2[,,3] <- subj2_scan3
  
  Y <- list()
  Y[[1]] <- subj1
  Y[[2]] <- subj2
  
  id_variable <- c("1", "1", "2", "2", "2")
  
  # Run the function
  result <- format_time_series_input(Y, longit_adjusted = TRUE, debug = FALSE,
                                     id_variable = id_variable)
  
  # Checking the result is correct
  check_subj1_scan1 <- all(result$Y_mixed_QxVNstar[1:2, 1:3]   == subj1_scan1)
  check_subj1_scan2 <- all(result$Y_mixed_QxVNstar[1:2, 4:6]   == subj1_scan2)
  check_subj2_scan1 <- all(result$Y_mixed_QxVNstar[1:2, 7:9]   == subj2_scan1)
  check_subj2_scan2 <- all(result$Y_mixed_QxVNstar[1:2, 10:12] == subj2_scan2)
  check_subj2_scan3 <- all(result$Y_mixed_QxVNstar[1:2, 13:15] == subj2_scan3)
  
  testthat::expect_true(all(check_subj1_scan1, check_subj1_scan2, check_subj2_scan1,
                            check_subj2_scan2, check_subj2_scan3))
  
})








test_that("3D Input Array Works", {
  
  true_Ji <- c(2, 3)
  subj1_scan1 <- matrix(c(1,2,3,4,5,6),     nrow = 2)
  subj1_scan2 <- matrix(c(1,2,3,4,5,6) + 1, nrow = 2)
  subj2_scan1 <- matrix(c(1,2,3,4,5,6) + 2, nrow = 2)
  subj2_scan2 <- matrix(c(1,2,3,4,5,6) + 3, nrow = 2)
  subj2_scan3 <- matrix(c(1,2,3,4,5,6) + 4, nrow = 2)
  
  Y <- array(0, dim = c(2, 3, 5))
  Y[,,1] <- subj1_scan1
  Y[,,2] <- subj1_scan2
  Y[,,3] <- subj2_scan1
  Y[,,4] <- subj2_scan2
  Y[,,5] <- subj2_scan3
  
  id_variable <- c("1", "1", "2", "2", "2")
  
  # Run the function
  result <- format_time_series_input(Y, longit_adjusted = TRUE, debug = FALSE,
                                     id_variable = id_variable)
  
  # Checking the result is correct
  check_subj1_scan1 <- all(result$Y_mixed_QxVNstar[1:2, 1:3]   == subj1_scan1)
  check_subj1_scan2 <- all(result$Y_mixed_QxVNstar[1:2, 4:6]   == subj1_scan2)
  check_subj2_scan1 <- all(result$Y_mixed_QxVNstar[1:2, 7:9]   == subj2_scan1)
  check_subj2_scan2 <- all(result$Y_mixed_QxVNstar[1:2, 10:12] == subj2_scan2)
  check_subj2_scan3 <- all(result$Y_mixed_QxVNstar[1:2, 13:15] == subj2_scan3)
  
  testthat::expect_true(all(check_subj1_scan1, check_subj1_scan2, check_subj2_scan1,
                            check_subj2_scan2, check_subj2_scan3))
  
})



test_that("4D Input Array Works", {
  
  true_Ji <- c(2, 3)
  subj1_scan1 <- matrix(c(1,2,3,4,5,6),     nrow = 2)
  subj1_scan2 <- matrix(c(1,2,3,4,5,6) + 1, nrow = 2)
  subj2_scan1 <- matrix(c(1,2,3,4,5,6) + 2, nrow = 2)
  subj2_scan2 <- matrix(c(1,2,3,4,5,6) + 3, nrow = 2)
  subj2_scan3 <- matrix(c(1,2,3,4,5,6) + 4, nrow = 2)
  
  # Q, V, N, Jmax
  Y <- array(0, dim = c(2, 3, 2, 3))
  Y[,,1,1] <- subj1_scan1
  Y[,,1,2] <- subj1_scan2
  Y[,,2,1] <- subj2_scan1
  Y[,,2,2] <- subj2_scan2
  Y[,,2,3] <- subj2_scan3
  
  id_variable <- c("1", "1", "2", "2", "2")
  
  # Run the function
  result <- format_time_series_input(Y, longit_adjusted = TRUE, debug = FALSE,
                                     id_variable = id_variable)
  
  # Checking the result is correct
  check_subj1_scan1 <- all(result$Y_mixed_QxVNstar[1:2, 1:3]   == subj1_scan1)
  check_subj1_scan2 <- all(result$Y_mixed_QxVNstar[1:2, 4:6]   == subj1_scan2)
  check_subj2_scan1 <- all(result$Y_mixed_QxVNstar[1:2, 7:9]   == subj2_scan1)
  check_subj2_scan2 <- all(result$Y_mixed_QxVNstar[1:2, 10:12] == subj2_scan2)
  check_subj2_scan3 <- all(result$Y_mixed_QxVNstar[1:2, 13:15] == subj2_scan3)
  
  testthat::expect_true(all(check_subj1_scan1, check_subj1_scan2, check_subj2_scan1,
                            check_subj2_scan2, check_subj2_scan3))
  
})