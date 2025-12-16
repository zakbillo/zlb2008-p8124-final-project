library(tidyverse)
library(glmnet)
library(pcalg)
library(foreach)
library(doParallel)
library(reshape2)
library(graph) 

X_normal_stable <- read_csv("data/X_normal_stable.csv")

source("lv-ida/lvida.R") 
source("lv-ida/iscyclic.R")

estimate_downstream_effects <- function(data, cause_var, n_boots = 100, title = NULL, file_name = NULL) {
  
  # Setup
  alphas <- c(0.01, 0.05, 0.10)
  cause_idx <- which(colnames(data) == cause_var)
  
  if (length(cause_idx) == 0) stop("Cause variable not found in data.")
  
  # "Predictors" here are actually potential TARGETS (Effects)
  potential_targets <- colnames(data)[-cause_idx]
  n_targets <- length(potential_targets)
  n_obs <- nrow(data)
  
  final_summary1 <- data.frame()
  
  # Title/File setup
  if (is.null(title)) title <- paste("Downstream Analysis | Cause:", cause_var)
  if (is.null(file_name)) {
    clean_cause <- gsub("[^[:alnum:]]", "_", cause_var)
    file_name <- paste0("normal_Downstream_Effects_", clean_cause, ".pdf")
    csv_name <- paste0("normal_Downstream_Effects_", clean_cause, ".csv")
  } else {
    # derive csv name from pdf name
    csv_name <- sub(".pdf$", ".csv", file_name)
  }
  
  pdf(file_name, width = 12, height = 10)
  
  for (alpha in alphas) {
    message(paste("\n=== Processing Alpha =", alpha, "==="))
    
    # Rows = Targets, Cols = Bootstrap Iterations
    boot_min_effects <- matrix(NA, nrow = n_targets, ncol = n_boots)
    boot_max_effects <- matrix(NA, nrow = n_targets, ncol = n_boots)
    
    # --- BOOTSTRAP LOOP ---
    for (b in 1:n_boots) {
      if (b %% 50 == 0) message(paste("  Bootstrap", b, "/", n_boots))
      
      # Resample
      boot_indices <- sample(1:n_obs, n_obs, replace = TRUE)
      boot_data <- data[boot_indices, ]
      
      # Re-calculate stats
      suffStat <- list(C = cor(boot_data, use = "complete.obs"), n = n_obs)
      mcov <- cov(boot_data, use = "complete.obs")
      
      # Run FCI
      pag.est <- tryCatch({
        fci(suffStat, indepTest = gaussCItest, alpha = alpha, labels = colnames(data), verbose = FALSE, m.max = 4)
      }, error = function(e) return(NULL))
      
      if (is.null(pag.est)) next
      
      amat <- pag.est@amat
      
      # Loop through potential DOWNSTREAM targets
      for (i in seq_along(potential_targets)) {
        target_name <- potential_targets[i]
        target_idx <- which(colnames(data) == target_name)
        
        # --- LOGIC FLIP: Check if Cause -> Target is possible ---
        # We check if 'target_idx' is in the descendants of 'cause_idx'
        is_possible_descendant <- target_idx %in% pcalg::possibleDe(amat, cause_idx)
        
        if (!is_possible_descendant) {
          boot_min_effects[i, b] <- 0
          boot_max_effects[i, b] <- 0
        } else {
          
          # --- DYNAMIC FUNCTION CALL ---
          tryCatch({
            if (exists("lvida", mode = "function")) {
              # args: x (cause), y (effect), cov, graph
              effs <- lvida(cause_idx, target_idx, mcov, amat)
            } else {
              effs <- lvIda(cause_idx, target_idx, mcov, amat)
            }
            
            if (!all(is.na(effs))) {
              boot_min_effects[i, b] <- min(abs(effs))
              boot_max_effects[i, b] <- max(abs(effs))
            } else {
              boot_min_effects[i, b] <- 0
              boot_max_effects[i, b] <- 0
            }
          }, error = function(e) {
            boot_min_effects[i, b] <- NA
            boot_max_effects[i, b] <- NA
          })
        }
      }
    } 
    
    # Summarize Results for this alpha
    alpha_summary <- data.frame(
      Target_Gene = potential_targets,
      Alpha = alpha,
      Stability_Pct = rowMeans(boot_max_effects > 0, na.rm = TRUE) * 100,
      Median_Min_Effect = apply(boot_min_effects, 1, function(x) median(x[x > 0], na.rm = TRUE)),
      Median_Max_Effect = apply(boot_max_effects, 1, function(x) median(x[x > 0], na.rm = TRUE))
    )
    
    alpha_summary[is.na(alpha_summary)] <- 0
    final_summary1 <- rbind(final_summary1, alpha_summary)
    
    # Reference Plot
    suffStat_full <- list(C = cor(data), n = nrow(data))
    pag_full <- fci(suffStat_full, indepTest = gaussCItest, alpha = alpha, labels = colnames(data))
    
    # Label top downstream targets
    top_targets <- alpha_summary[order(-alpha_summary$Stability_Pct), ][1:3, ]
    caption <- paste("Top Downstream Targets:", paste(top_targets$Target_Gene, "(", round(top_targets$Stability_Pct,0), "%)", collapse=", "))
    
    plot(pag_full, main = paste0(title, "\nalpha=", alpha, "\n", caption))
    
  } 
  
  dev.off()
  
  # Save to CSV 
  write.csv(final_summary1, csv_name, row.names = FALSE)
  message(paste("Saved edge list to:", csv_name))
  
  return(final_summary1)
}

results_normal <- estimate_downstream_effects(X_normal_stable, cause_var = "HLA-DQB1", n_boots = 500)
