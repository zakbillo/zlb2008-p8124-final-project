library(tidyverse)
library(glmnet)
library(pcalg)
library(foreach)
library(doParallel)
library(reshape2)
library(graph) 

analyze_pag_sensitivity <- function(data, title = NULL, file_name = NULL) {
  
  # Title and filename setup
  data_arg_name <- deparse(substitute(data)) 
  if (is.null(title)) {
    clean_name <- if(data_arg_name == ".") "Data" else data_arg_name
    title <- paste("Estimated PAG:", clean_name)
  }
  
  if (is.null(file_name)) {
    clean_title <- gsub("[^[:alnum:]]", "_", title)
    file_name <- paste0(clean_title, ".pdf")
  }
  
  # Define alphas
  alphas <- c(0.01, 0.05, 0.10)
  
  # Initialize results df
  results_df1 <- data.frame(Alpha = numeric(), 
                            Edge_Count = integer(), 
                            Time_Sec = numeric())
  
  # suffStat
  suffStat <- list(C = cor(data, use = "pairwise.complete.obs"), n = nrow(data))
  
  pdf(file_name, width = 10, height = 10)
  
  for (alpha in alphas) {
    
    message(paste("Running FCI for Alpha:", alpha, "..."))
    
    # system.time
    pag.est <- NULL
    
    timer <- system.time({
      pag.est <- fci(suffStat, indepTest = gaussCItest, alpha = alpha, 
                     labels = colnames(data), verbose = FALSE)
    })
    
    # Extract elapsed time
    duration <- timer["elapsed"]
    
    # Count Edges
    amat <- pag.est@amat
    current_edges <- sum(amat[upper.tri(amat)] != 0)
    
    # Store Results
    results_df1 <- rbind(results_df1, data.frame(Alpha = alpha, 
                                                 Edge_Count = current_edges, 
                                                 Time_Sec = as.numeric(duration)))
    
    # Plot PAG
    plot(pag.est)
    
    # Title
    title_text <- paste0(title, 
                         "\nALPHA: ", alpha, 
                         "   |   Edges: ", current_edges, 
                         "   |   Time: ", round(duration, 3), "s")
    
    title(main = title_text, col.main = "blue", cex.main = 1.5)
  }
  
  dev.off()
  
  return(results_df1)
}

results_normal <- analyze_pag_sensitivity(X_normal_top, 
                                          title = "normal PAG Analysis", 
                                          file_name = "normal_PAG_Analysis.pdf")

print(results_normal)

# Visualization of the runtime
plot(results_normal$Alpha, results_normal$Time_Sec, 
     type = "b", pch = 19, col = "red",
     main = "Computation Time by Alpha",
     xlab = "Alpha", ylab = "Time (seconds)")