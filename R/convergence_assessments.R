#' Plot relative changes at each CAVI iteration for key model parameters
#' 
#' @param obj an object of \link{class} "BBSSModel"
#' @param plot_type Either "Variance" or "DPM", depending on which parameters you wish to view
#' 
#' @details
#' This function generates a ggplot object showing relative change over each iteration for either the variance parameters or the DPM parameters, depending on user selection.
#' 
#' @return A ggplot object
#' 
#' @export
plot_relative_changes <- function(obj, plot_type = "Variance"){
  
  # Extract the relative change tracking information
  RCT <- obj$relative_change_tracking
  
  # Start for storing plot information
  plot_tibble <- NULL
  plt <- NULL
  
  if (plot_type == "Variance"){
    plot_tibble <- tibble(
      error_variance   = RCT$relchange_sigma_sq_error,
      site_variance    = RCT$relchange_sigma_sq_site,
      global_precision = RCT$relchange_tau_sq
    ) %>%
      mutate(Iteration = row_number()) 
    
    plt <- plot_tibble %>%
      group_by(Iteration) %>%
      pivot_longer(-c("Iteration")) %>%
      ggplot(aes(x = Iteration, y = value)) +
      theme_minimal(base_size = 10) + 
      geom_point() + 
      geom_line() +
      ylab("Relative Change") + 
      facet_wrap(vars(name), scales = "free")
    
  }
  
  if (plot_type == "DPM"){
    
    # DPM cluster means
    mu_h <- RCT$relchange_mu_h %>% 
      reshape2::melt() %>% 
      rename("Group" = "Var1", "Iteration" = "Var2") %>%
      as_tibble() %>%
      mutate(label = "DPM Cluster Mean") %>%
      mutate(Group = factor(Group))
    
    sigma_h_sq <- RCT$relchange_sigma_h_sq %>% 
      reshape2::melt() %>% 
      rename("Group" = "Var1", "Iteration" = "Var2") %>%
      as_tibble() %>%
      mutate(label = "DPM Cluster Variance")  %>%
      mutate(Group = factor(Group))
    
    plot_data <- bind_rows(mu_h, sigma_h_sq)
    
    plt <- plot_data %>%
      ggplot(aes(x = Iteration, y = value, color = Group)) +
      theme_minimal(base_size = 10) + 
      geom_point() + 
      geom_line() +
      ylab("Relative Change") + 
      facet_wrap(vars(label), scales = "free")
    
  }
  
  return(plt)
  
  
}