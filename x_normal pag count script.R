library(tidyverse)
library(glmnet)
library(pcalg)
library(foreach)
library(doParallel)
library(reshape2)

X_normal_top25 <- read_csv("data/X_normal_top25.csv")

pag_plot <- function(data, title = NULL, file_name = NULL) {
  
  # Title setup
  data_arg_name <- deparse(substitute(data)) 
  if (is.null(title)) {
    clean_name <- if(data_arg_name == ".") "Data" else data_arg_name
    title <- paste("Estimated PAG:", clean_name)
  }
  
  alphas <- c(0.01, 0.05, 0.10)
  
  # Store results in a df
  edge_counts <- data.frame(Alpha = numeric(), Edge_Count = integer())
  
  # Filename setup
  if (is.null(file_name)) {
    clean_title <- gsub("[^[:alnum:]]", "_", title)
    file_name <- paste0(clean_title, ".pdf")
  }
  
  suffStat <- list(C = cor(data), n = nrow(data))
  
  pdf(file_name, width = 10, height = 10)
  
  for (alpha in alphas) {
    
    # Run FCI
    pag.est <- fci(suffStat, indepTest = gaussCItest, alpha = alpha, labels = colnames(data))
    
    # Extract amat
    amat <- pag.est@amat
    
    # Count non-zero entries in the upper triangle
    current_edges <- sum(amat[upper.tri(amat)] != 0)
    
    # Append to results
    edge_counts <- rbind(edge_counts, data.frame(Alpha = alpha, Edge_Count = current_edges))
    
    # Plot
    page_title <- paste0(title, " (alpha = ", alpha, " | Edges = ", current_edges, ")")
    plot(pag.est, main = page_title)
  }
  
  dev.off()
  
  # Return the table of counts
  return(edge_counts)
}

counts_normal <- pag_plot(X_normal_top25)
print(counts_normal) # 12 at 0.01, 12 at 0.05, 14 at 0.10