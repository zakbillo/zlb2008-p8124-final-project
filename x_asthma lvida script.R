library(tidyverse)
library(glmnet)
library(pcalg)
library(foreach)
library(doParallel)
library(reshape2)
library(graph)

# Load Data
X_asthma_stable <- read_csv("data/X_asthma_stable.csv", show_col_types = FALSE)
X_asthma_stable <- X_asthma_stable[, -1]

# --- PARALLEL FUNCTION ---
estimate_downstream_effects_safe_parallel <- function(data, cause_var, n_boots = 100, file_name = NULL) {
  
  # 1. Setup
  alphas <- c(0.01)
  cause_idx <- which(colnames(data) == cause_var)
  
  if (length(cause_idx) == 0) stop("Cause variable not found in data.")
  
  potential_targets <- colnames(data)[-cause_idx]
  n_targets <- length(potential_targets)
  n_obs <- nrow(data)
  node_names <- colnames(data)
  
  # 2. File Name Setup
  clean_cause <- gsub("[^[:alnum:]]", "_", cause_var)
  if (is.null(file_name)) {
    csv_name <- paste0("Downstream_Effects_", clean_cause, ".csv")
  } else {
    csv_name <- ifelse(grepl("\\.csv$", file_name), file_name, paste0(file_name, ".csv"))
  }
  
  final_summary <- data.frame()
  
  # --- START PARALLEL CLUSTER ---
  n_cores <- parallel::detectCores() - 1 
  cl <- makeCluster(n_cores) 
  registerDoParallel(cl)
  message(paste("Running on", n_cores, "cores..."))
  
  # CRITICAL FIX: Source the scripts on EVERY worker node
  # This ensures all hidden helper functions are available
  clusterEvalQ(cl, {
    if(file.exists("lv-ida/lvida.R")) source("lv-ida/lvida.R") 
    if(file.exists("lv-ida/iscyclic.R")) source("lv-ida/iscyclic.R")
    library(pcalg)
    library(graph)
  })
  
  # 3. Main Loop (Alphas)
  for (alpha in alphas) {
    message(paste("\n=== Processing Alpha =", alpha, "==="))
    
    # 4. PARALLEL BOOTSTRAP LOOP
    # We remove .export since we sourced everything above
    boot_results <- foreach(b = 1:n_boots, 
                            .packages = c("pcalg", "graph", "methods")) %dopar% {
                              
                              # --- INSIDE WORKER ---
                              # Resample
                              boot_indices <- sample(1:n_obs, n_obs, replace = TRUE)
                              boot_data <- data[boot_indices, ]
                              
                              # Stats
                              suffStat <- list(C = cor(boot_data, use = "pairwise.complete.obs"), n = n_obs)
                              mcov <- cov(boot_data, use = "pairwise.complete.obs")
                              
                              # Run RFCI
                              pag.est <- tryCatch({
                                rfci(suffStat, indepTest = gaussCItest, alpha = alpha, 
                                     labels = node_names, verbose = FALSE, m.max = 4)
                              }, error = function(e) { return(NULL) })
                              
                              if (is.null(pag.est)) return(NULL)
                              
                              amat <- pag.est@amat
                              
                              # Check Cycle
                              if (exists("is.cyclic") && is.cyclic(amat)) return(NULL)
                              
                              row_min <- rep(0, n_targets)
                              row_max <- rep(0, n_targets)
                              
                              descendants <- pcalg::possibleDe(amat, cause_idx)
                              
                              for (i in 1:n_targets) {
                                target_name <- potential_targets[i]
                                target_idx <- which(node_names == target_name)
                                
                                if (target_idx %in% descendants) {
                                  effs <- tryCatch({
                                    # Ensure we use the correct function name (lv.ida vs lvida)
                                    if (exists("lv.ida")) {
                                      lv.ida(cause_idx, target_idx, mcov, amat, method = "local")
                                    } else {
                                      lvida(cause_idx, target_idx, mcov, amat)
                                    }
                                  }, error = function(e) {
                                    # DEBUG: Return a special code or print if needed
                                    # For now, return NA but we know the environment is fixed
                                    return(NA) 
                                  })
                                  
                                  if (!all(is.na(effs))) {
                                    row_min[i] <- min(abs(effs))
                                    row_max[i] <- max(abs(effs))
                                  }
                                }
                              }
                              
                              return(list(min = row_min, max = row_max))
                            }
    
    # 5. Process Results
    valid_results <- boot_results[!sapply(boot_results, is.null)]
    n_valid <- length(valid_results)
    message(paste("Successful Bootstraps:", n_valid, "/", n_boots))
    
    if (n_valid == 0) {
      message("Warning: All bootstraps failed. Check 'lv.ida' compatibility.")
      next
    }
    
    boot_min_effects <- do.call(cbind, lapply(valid_results, function(x) x$min))
    boot_max_effects <- do.call(cbind, lapply(valid_results, function(x) x$max))
    
    alpha_summary <- data.frame(
      Target_Gene = potential_targets,
      Alpha = alpha,
      Stability_Pct = rowMeans(boot_max_effects > 0, na.rm = TRUE) * 100,
      Median_Min_Effect = apply(boot_min_effects, 1, function(x) median(x[x > 0], na.rm = TRUE)),
      Median_Max_Effect = apply(boot_max_effects, 1, function(x) median(x[x > 0], na.rm = TRUE)),
      CI_Lower          = apply(boot_min_mat, 1, function(x) quantile(x, 0.025, na.rm=T)),
      CI_Upper          = apply(boot_min_mat, 1, function(x) quantile(x, 0.975, na.rm=T))
    )
    
    alpha_summary[is.na(alpha_summary)] <- 0
    final_summary <- rbind(final_summary, alpha_summary)
  }
  
  stopCluster(cl)
  
  write.csv(final_summary, csv_name, row.names = FALSE)
  return(final_summary)
}

# Run
results_asthmaa <- estimate_downstream_effects_safe_parallel(X_asthma_stable, "MAP3K5-AS2", n_boots = 100, file_name = "Asthmaa_MAP3K5_AS2.csv")