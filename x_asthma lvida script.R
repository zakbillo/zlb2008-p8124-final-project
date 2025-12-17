library(tidyverse)
library(glmnet)
library(pcalg)
library(foreach)
library(doParallel)
library(reshape2)
library(graph) 

X_asthma_stable <- read_csv("data/X_asthma_stable.csv")

source("lv-ida/lvida.R") 
source("lv-ida/iscyclic.R")

estimate_downstream_effects_safe <- function(data, cause_var, n_boots = 100, file_name = NULL) {
  
  # 1. Setup
  alphas <- c(0.05)
  cause_idx <- which(colnames(data) == cause_var)
  
  if (length(cause_idx) == 0) stop("Cause variable not found in data.")
  
  potential_targets <- colnames(data)[-cause_idx]
  n_targets <- length(potential_targets)
  n_obs <- nrow(data)
  
  final_summary <- data.frame()
  
  # 2. File Name Setup
  clean_cause <- gsub("[^[:alnum:]]", "_", cause_var)
  if (is.null(file_name)) {
    csv_name <- paste0("Downstream_Effects_", clean_cause, ".csv")
  } else {
    csv_name <- ifelse(grepl("\\.csv$", file_name), file_name, paste0(file_name, ".csv"))
  }
  
  dir_path <- dirname(csv_name)
  if (dir_path != "." && !dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
  
  # 3. Main Loop
  for (alpha in alphas) {
    message(paste("\n=== Processing Alpha =", alpha, "==="))
    
    boot_min_effects <- matrix(NA, nrow = n_targets, ncol = n_boots)
    boot_max_effects <- matrix(NA, nrow = n_targets, ncol = n_boots)
    
    skipped_count <- 0
    
    for (b in 1:n_boots) {
      if (b %% 10 == 0) message(paste("   Bootstrap", b, "/", n_boots))
      
      # Resample
      boot_indices <- sample(1:n_obs, n_obs, replace = TRUE)
      boot_data <- data[boot_indices, ]
      
      # Stats
      suffStat <- list(C = cor(boot_data, use = "pairwise.complete.obs"), n = n_obs)
      mcov <- cov(boot_data, use = "pairwise.complete.obs")
      
      # Run RFCI
      pag.est <- tryCatch({
        rfci(suffStat, indepTest = gaussCItest, alpha = alpha, 
            labels = colnames(data), verbose = FALSE, m.max = 4)
      }, error = function(e) return(NULL))
      
      if (is.null(pag.est)) {
        skipped_count <- skipped_count + 1
        next
      }
      
      amat <- pag.est@amat
      
      # --- CRITICAL FIX: CHECK FOR CYCLES ---
      # If the graph has a cycle, lv.ida will crash with Stack Overflow.
      # We must skip this iteration.
      if (exists("is.cyclic") && is.cyclic(amat)) {
        # message("Skipping cyclic graph...") # Optional: uncomment to see how often this happens
        skipped_count <- skipped_count + 1
        next
      }
      
      # Check Targets
      for (i in seq_along(potential_targets)) {
        target_name <- potential_targets[i]
        target_idx <- which(colnames(data) == target_name)
        
        # Optimization
        is_possible_descendant <- target_idx %in% pcalg::possibleDe(amat, cause_idx)
        
        if (!is_possible_descendant) {
          boot_min_effects[i, b] <- 0
          boot_max_effects[i, b] <- 0
        } else {
          # Run LV-IDA
          # Note: Checking if function is 'lv.ida' (dot) or 'lvida' (no dot) based on your source file
          effs <- tryCatch({
            if (exists("lv.ida")) {
              lv.ida(cause_idx, target_idx, mcov, amat, method = "local")
            } else {
              lvida(cause_idx, target_idx, mcov, amat) 
            }
          }, error = function(e) return(NA))
          
          if (!all(is.na(effs))) {
            boot_min_effects[i, b] <- min(abs(effs))
            boot_max_effects[i, b] <- max(abs(effs))
          } else {
            boot_min_effects[i, b] <- 0
            boot_max_effects[i, b] <- 0
          }
        }
      }
    } 
    
    message(paste("   Skipped", skipped_count, "cyclic/failed graphs out of", n_boots))
    
    # 4. Summarize
    alpha_summary <- data.frame(
      Target_Gene = potential_targets,
      Alpha = alpha,
      Stability_Pct = rowMeans(boot_max_effects > 0, na.rm = TRUE) * 100,
      Median_Min_Effect = apply(boot_min_effects, 1, function(x) median(x[x > 0], na.rm = TRUE)),
      Median_Max_Effect = apply(boot_max_effects, 1, function(x) median(x[x > 0], na.rm = TRUE))
    )
    
    alpha_summary[is.na(alpha_summary)] <- 0
    final_summary <- rbind(final_summary, alpha_summary)
  } 
  
  write.csv(final_summary, csv_name, row.names = FALSE)
  message(paste("Saved downstream effects summary to:", csv_name))
  return(final_summary)
}

# Run
results_asthma <- estimate_downstream_effects_safe(X_asthma_stable, "MAP3K5-AS2", n_boots = 100, file_name = "Asthma_MAP3K5_AS2.csv")