library(tidyverse)
library(glmnet)
library(pcalg)
library(reshape2)
library(graph)

# --- Setup ---
if(!file.exists("lv-ida/lvida.R")) stop("lv-ida/lvida.R not found.")
source("lv-ida/lvida.R")
if(file.exists("lv-ida/iscyclic.R")) source("lv-ida/iscyclic.R")

X_asthma_stable <- read_csv("data/X_asthma_stable.csv", show_col_types = FALSE)
X_asthma_stable <- X_asthma_stable[, -1]

# --- CORE FUNCTION (SERIAL) ---
estimate_downstream_effects_serial <- function(data, cause_var, n_boots = 100, file_name = NULL) {
  
  # 1. Setup and Pre-indexing
  alpha_val <- 0.01  
  node_names <- colnames(data)
  cause_idx <- which(node_names == cause_var)
  
  if (length(cause_idx) == 0) stop("Cause variable not found in data.")
  
  potential_targets <- node_names[-cause_idx]
  target_indices <- setdiff(seq_along(node_names), cause_idx)
  n_targets <- length(target_indices)
  n_obs <- nrow(data)
  
  # 2. File Name Setup
  clean_cause <- gsub("[^[:alnum:]]", "_", cause_var)
  csv_name <- if(is.null(file_name)) paste0("Causal_Stability_", clean_cause, ".csv") else file_name
  
  message(paste("Processing", n_boots, "bootstraps serially..."))
  
  # 3. Serial Bootstrap Loop
  boot_results <- vector("list", n_boots)
  pb <- txtProgressBar(min = 0, max = n_boots, style = 3)
  
  for (b in 1:n_boots) {
    
    # Resample Data
    boot_indices <- sample.int(n_obs, n_obs, replace = TRUE)
    boot_data <- data[boot_indices, ]
    suffStat <- list(C = cor(boot_data, use = "pairwise.complete.obs"), n = n_obs)
    mcov <- cov(boot_data, use = "pairwise.complete.obs")
    
    # Estimate PAG via RFCI
    pag_est <- tryCatch({
      rfci(suffStat, indepTest = gaussCItest, alpha = alpha_val, 
           labels = node_names, verbose = FALSE, m.max = 4)
    }, error = function(e) {
      message(paste("\nBootstrap", b, "failed at RFCI"))
      return(NULL)
    })
    
    if (is.null(pag_est)) {
      boot_results[[b]] <- NULL
      setTxtProgressBar(pb, b)
      next
    }
    
    amat <- pag_est@amat
    
    # Proper cycle checking with tryCatch
    is_cyclic <- tryCatch({
      if (exists("is.cyclic", mode = "function")) {
        is.cyclic(amat)
      } else {
        FALSE  # If function doesn't exist, proceed without check
      }
    }, error = function(e) {
      FALSE  # If error occurs, proceed without check
    })
    
    if (is_cyclic) {
      message(paste("\nBootstrap", b, "skipped - cyclic graph detected"))
      boot_results[[b]] <- NULL
      setTxtProgressBar(pb, b)
      next
    }
    
    # Calculate Causal Effects
    row_eff <- numeric(n_targets)
    descendants <- pcalg::possibleDe(amat, cause_idx)
    
    # Optimization: Only check nodes reachable from the cause
    active_targets <- intersect(target_indices, descendants)
    
    if (length(active_targets) > 0) {
      for (target_idx in active_targets) {
        # Positional mapping in the output vector (1 to n_targets)
        out_pos <- if(target_idx < cause_idx) target_idx else target_idx - 1
        
        effs <- tryCatch({
          if (exists("lv.ida", mode = "function")) {
            lv.ida(cause_idx, target_idx, mcov, amat, method = "local")
          } else {
            lvida(cause_idx, target_idx, mcov, amat)
          }
        }, error = function(e) {
          message(paste("\nBootstrap", b, "- Effect calculation failed for target", target_idx))
          return(NA)
        })
        
        # Take the MEDIAN of the multiset for this specific bootstrap
        if (!all(is.na(effs)) && length(effs) > 0) {
          row_eff[out_pos] <- median(effs, na.rm = TRUE)
        }
      }
    }
    
    boot_results[[b]] <- row_eff
    setTxtProgressBar(pb, b)
  }
  
  close(pb)
  
  # 4. Result Aggregation and CI Calculation
  # Filter out failed bootstraps
  valid_results <- boot_results[!sapply(boot_results, is.null)]
  n_valid <- length(valid_results)
  message(paste("\nSuccessful Bootstraps:", n_valid, "/", n_boots))
  
  if (n_valid < 10) stop("Insufficient successful bootstraps for statistics.")
  
  # Create effect matrix [genes x bootstraps]
  eff_matrix <- do.call(cbind, valid_results)
  
  # Helper to compute stats for each gene across bootstraps
  # Only considers instances where a causal path was found (non-zero)
  compute_stats <- function(row) {
    causal_vals <- row[row != 0]
    if (length(causal_vals) == 0) return(c(0, 0, 0))
    
    c(median(causal_vals), 
      quantile(causal_vals, 0.025), 
      quantile(causal_vals, 0.975))
  }
  
  stats_data <- t(apply(eff_matrix, 1, compute_stats))
  
  final_summary <- data.frame(
    Target_Gene   = potential_targets,
    Alpha         = alpha_val,
    Stability_Pct = (rowSums(eff_matrix != 0) / n_valid) * 100,
    Median_Effect = stats_data[, 1],
    CI_Lower      = stats_data[, 2],
    CI_Upper      = stats_data[, 3]
  )
  
  # Save and Return
  write.csv(final_summary, csv_name, row.names = FALSE)
  return(final_summary)
}

# --- EXAMPLE USAGE ---
results_asthma <- estimate_downstream_effects_serial(X_asthma_stable, "MAP3K5-AS2", n_boots = 100, "asthma lvida.csv")