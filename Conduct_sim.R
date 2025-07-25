################################################################################
# 1. Setup: Clear Workspace, Set Working Directory, and Load Functions
################################################################################
rm(list = ls())
mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

# Load required functions
source("Function.R")
source("LCCA_function.R")

################################################################################
# 2. Simulation Parameters and Implementation
################################################################################
# Simulation parameters
rho      <- 0.8      # Correlation
comp     <- 4        # Number of components
time_num <- 8        # maximum of time points
p_ture   <- 10       # Number of important variables

sample_list <- c(50, 100, 150)  # Sample sizes
ratio       <- c(0.5, 1, 1.5)    # p/n ratios, where p is the dimension

# Run simulations for each combination of sample size and ratio
for(n in sample_list){
  for(r in ratio){
    p_X <- p_Y <- n * r
    source("Simulation.R")
  }
}

################################################################################
# 3. Correlation Error Analysis and LaTeX Table Output
################################################################################
latex_res <- NULL
for(i in 1:3){
  cor_error <- NULL
  for(k in 1:3){
    n    <- sample_list[i]
    r    <- ratio[k]
    p_X  <- p_Y <- n * r
    load(paste0("Result/Result_", n, "_", p_X, ".rda"))
    
    # Combine correlation errors over all simulation replicates
    result <- sapply(1:length(Result), function(h) {
      unlist(Result[[h]]$cor_error)
    })
    result   <- apply(result, 1, function(x) mean(x, na.rm = TRUE))
    cor_error <- rbind(cor_error, result)
  }
  
  # Reorganize correlation errors into a summary matrix
  f <- rbind(
    c(c(cor_error[, 3]), c(cor_error[, 8])),
    c(c(cor_error[, 2]), c(cor_error[, 7])),
    c(c(cor_error[, 4]), c(cor_error[, 9])),
    c(c(cor_error[, 5]), c(cor_error[, 10])),
    c(c(cor_error[, 1]), c(cor_error[, 6]))
  )
  latex_res <- rbind(latex_res, f)
}
latex_res <- data.frame(rep(c("LCCA", "SCCA", "MFPCA", "FSCCA", "CCVA"), 3), latex_res)
xtable::xtable(latex_res)

################################################################################
# 4. Variable Selection Error Analysis
################################################################################
var_error <- NULL
for(i in 1:3){
  for(k in 1:3){
    n    <- sample_list[i]
    r    <- ratio[k]
    p_X  <- p_Y <- n * r
    load(paste0("Result/Result_", n, "_", p_X, ".rda"))
    
    result <- sapply(1:length(Result), function(h) {
      # --- CCVA ---
      FNR_X_SLCCA <- sum(sapply(1:p_ture, function(j) {
        sum(j == Result[[h]]$Result_SLCCA$var_X[[1]])
      }) == 0) / p_ture
      
      FPR_X_SLCCA <- sum(sapply((p_ture + 1):p_X, function(j) {
        sum(j == Result[[h]]$Result_SLCCA$var_X[[1]])
      }) == 1) / (p_X - p_ture)
      
      Ove_X_SLCCA <- (sum(sapply(1:p_ture, function(j) {
        sum(j == Result[[h]]$Result_SLCCA$var_X[[1]])
      }) == 0) +
        sum(sapply((p_ture + 1):p_X, function(j) {
          sum(j == Result[[h]]$Result_SLCCA$var_X[[1]])
        }) == 1))
      Ove_X_SLCCA <- 1 - Ove_X_SLCCA / p_X
      
      FNR_Y_SLCCA <- sum(sapply(1:p_ture, function(j) {
        sum(j == Result[[h]]$Result_SLCCA$var_Y[[1]])
      }) == 0) / p_ture
      
      FPR_Y_SLCCA <- sum(sapply((p_ture + 1):p_Y, function(j) {
        sum(j == Result[[h]]$Result_SLCCA$var_Y[[1]])
      }) == 1) / (p_Y - p_ture)
      
      Ove_Y_SLCCA <- (sum(sapply(1:p_ture, function(j) {
        sum(j == Result[[h]]$Result_SLCCA$var_Y[[1]])
      }) == 0) +
        sum(sapply((p_ture + 1):p_Y, function(j) {
          sum(j == Result[[h]]$Result_SLCCA$var_Y[[1]])
        }) == 1))
      Ove_Y_SLCCA <- 1 - Ove_Y_SLCCA / p_Y
      
      # --- SCCA ---
      FNR_X_SCCA <- sum(sapply(1:p_ture, function(j) {
        sum(j == Result[[h]]$Result_GSCCA$var_X)
      }) == 0) / p_ture
      
      FPR_X_SCCA <- sum(sapply((p_ture + 1):p_X, function(j) {
        sum(j == Result[[h]]$Result_GSCCA$var_X)
      }) == 1) / (p_X - p_ture)
      
      Ove_X_SCCA <- (sum(sapply(1:p_ture, function(j) {
        sum(j == Result[[h]]$Result_GSCCA$var_X)
      }) == 0) +
        sum(sapply((p_ture + 1):p_X, function(j) {
          sum(j == Result[[h]]$Result_GSCCA$var_X)
        }) == 1))
      Ove_X_SCCA <- 1 - Ove_X_SCCA / p_X
      
      FNR_Y_SCCA <- sum(sapply(1:p_ture, function(j) {
        sum(j == Result[[h]]$Result_GSCCA$var_Y)
      }) == 0) / p_ture
      
      FPR_Y_SCCA <- sum(sapply((p_ture + 1):p_Y, function(j) {
        sum(j == Result[[h]]$Result_GSCCA$var_Y)
      }) == 1) / (p_Y - p_ture)
      
      Ove_Y_SCCA <- (sum(sapply(1:p_ture, function(j) {
        sum(j == Result[[h]]$Result_GSCCA$var_Y)
      }) == 0) +
        sum(sapply((p_ture + 1):p_Y, function(j) {
          sum(j == Result[[h]]$Result_GSCCA$var_Y)
        }) == 1))
      Ove_Y_SCCA <- 1 - Ove_Y_SCCA / p_Y
      
      # --- FSCCA ---
      FNR_X_FSCCA <- sum(sapply(1:p_ture, function(j) {
        sum(j == Result[[h]]$Result_FSCCA$var_X)
      }) == 0) / p_ture
      
      FPR_X_FSCCA <- sum(sapply((p_ture + 1):p_X, function(j) {
        sum(j == Result[[h]]$Result_FSCCA$var_X)
      }) == 1) / (p_X - p_ture)
      
      Ove_X_FSCCA <- (sum(sapply(1:p_ture, function(j) {
        sum(j == Result[[h]]$Result_FSCCA$var_X)
      }) == 0) +
        sum(sapply((p_ture + 1):p_X, function(j) {
          sum(j == Result[[h]]$Result_FSCCA$var_X)
        }) == 1))
      Ove_X_FSCCA <- 1 - Ove_X_FSCCA / p_X
      
      FNR_Y_FSCCA <- sum(sapply(1:p_ture, function(j) {
        sum(j == Result[[h]]$Result_FSCCA$var_Y)
      }) == 0) / p_ture
      
      FPR_Y_FSCCA <- sum(sapply((p_ture + 1):p_Y, function(j) {
        sum(j == Result[[h]]$Result_FSCCA$var_Y)
      }) == 1) / (p_Y - p_ture)
      
      Ove_Y_FSCCA <- (sum(sapply(1:p_ture, function(j) {
        sum(j == Result[[h]]$Result_FSCCA$var_Y)
      }) == 0) +
        sum(sapply((p_ture + 1):p_Y, function(j) {
          sum(j == Result[[h]]$Result_FSCCA$var_Y)
        }) == 1))
      Ove_Y_FSCCA <- 1 - Ove_Y_FSCCA / p_Y
      
      c(FNR_X_SLCCA, FPR_X_SLCCA, Ove_X_SLCCA,
        FNR_Y_SLCCA, FPR_Y_SLCCA, Ove_Y_SLCCA,
        FNR_X_SCCA, FPR_X_SCCA, Ove_X_SCCA,
        FNR_Y_SCCA, FPR_Y_SCCA, Ove_Y_SCCA,
        FNR_X_FSCCA, FPR_X_FSCCA, Ove_X_FSCCA,
        FNR_Y_FSCCA, FPR_Y_FSCCA, Ove_Y_FSCCA)
    })
    
    # Summarize results by computing 2.5%, mean, and 97.5% quantiles
    result <- cbind(
      apply(result, 1, function(x) quantile(x, 0.025)),
      apply(result, 1, mean),
      apply(result, 1, function(x) quantile(x, 0.975))
    )
    
    var_error <- rbind(var_error, 
                       data.frame(result * 300,
                                  Method = rep(c(rep("CCVA", 6), rep("SCCA", 6), rep("FSCCA", 6)), 1),
                                  Label  = rep(c(rep("The first dataset", 3), rep("The second dataset", 3)), 3),
                                  Type   = rep(c("FNR", "FPR", "OER"), 6),
                                  n      = rep(paste0("n=", n), 18),
                                  p      = rep(r, 18))
    )
  }
}

# Format factors and column names for plotting
var_error$n <- factor(var_error$n, levels = paste0("n=", sample_list))
colnames(var_error)[1:3] <- c("low", "mean", "upper")
var_error$Method <- factor(var_error$Method, levels = c("SCCA", "FSCCA", "CCVA"))

################################################################################
# 5. Plotting: Variable Selection Errors (FNR and FPR)
################################################################################
library(ggplot2)
library(patchwork)

# Plot for False Negative Rate (FNR)
plot_dat <- var_error[var_error$Type == "FNR", ]
p1 <- ggplot(plot_dat) +
  geom_bar(aes(x = p, y = mean, group = interaction(Method, Label), fill = Method),
           stat = "identity", position = "dodge", width = 0.3) +
  facet_wrap(n ~ Label, nrow = 3, scales = "free_y") +
  labs(x = "p / n", y = "",
       title = "(B) False Negative Rate (%)") +
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2)) +
  # scale_y_log10(limits = c(1, 100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("blue", "orange", "red"))
p1

# Plot for False Positive Rate (FPR)
plot_dat <- var_error[var_error$Type == "FPR", ]
p2 <- ggplot(plot_dat) +
  geom_bar(aes(x = p, y = mean, group = interaction(Method, Label), fill = Method),
           stat = "identity", position = "dodge", width = 0.3) +
  facet_wrap(n ~ Label, nrow = 3) +
  labs(x = "p / n", y = "",
       title = "(A) False Positive Rate (%)") +
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2)) +
  # scale_y_log10(
  #   # limits = c(0.01, NA),         # <-- start at y = 10; NA means “auto” upper limit
  #   expand = c(0, 0),           # optional: remove padding
  # ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("blue", "orange", "red"))
p2

combined_var_plot <- p1 + p2 & theme(legend.position = "bottom")
combined_var_plot + plot_layout(guides = "collect", ncol = 2)
ggsave(paste0("Figure/", "CCA_var_sim", ".pdf"), width = 9, height = 5, dpi = 300)

################################################################################
# 6. CCA Loading Error Analysis and Plotting
################################################################################
loa_error <- NULL
for(i in 1:3){
  for(k in 1:3){
    n    <- sample_list[i]
    r    <- ratio[k]
    p_X  <- p_Y <- n * r
    load(paste0("Result/Result_", n, "_", p_X, ".rda"))
    
    result <- sapply(1:length(Result), function(h) {
      unlist(Result[[h]]$loa_error) / p_X
    })
    
    loa_error <- rbind(loa_error, 
                       data.frame(Error  = c(t(result)),
                                  Method = rep(c(rep("CCVA", 300), rep("LCCA", 300), rep("FSCCA", 300), rep("MFPCA", 300)), 2),
                                  Label  = c(rep("The first dataset", 4 * 300), rep("The second dataset", 4 * 300)),
                                  n      = rep(paste0("n=", n), length(c(t(result)))),
                                  p      = rep(r, length(c(t(result))))
                       )
    )
  }
}

loa_error$n <- factor(loa_error$n, levels = paste0("n=", sample_list))
loa_error$Method <- factor(loa_error$Method, levels = c("LCCA", "FSCCA", "MFPCA", "CCVA"))
loa_error$p <- factor(loa_error$p)

# Plot loading errors for the first dataset
plot_dat <- loa_error[loa_error$Label == "The first dataset", ]
p1_load <- ggplot(plot_dat) +
  geom_boxplot(aes(x = p, y = Error, fill = Method), position = "dodge", size = 0.2) +
  facet_wrap(~ n, ncol = 1, scales = "free_y") +
  labs(x = "p / n", y = "Errors of functional loading", title = "(A)") +
  scale_y_log10() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("#31a354", "orange", "purple", "red"))

# Plot loading errors for the second dataset
plot_dat <- loa_error[loa_error$Label == "The second dataset", ]
p2_load <- ggplot(plot_dat) +
  geom_boxplot(aes(x = p, y = Error, fill = Method), position = "dodge", size = 0.2) +
  facet_wrap(~ n, ncol = 1, scales = "free_y") +
  labs(x = "p / n", y = "Errors of functional loading", title = "(B)") +
  scale_y_log10() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("#31a354", "orange", "purple", "red"))

combined_load_plot <- p1_load + p2_load & theme(legend.position = "bottom")
combined_load_plot + plot_layout(guides = "collect", ncol = 2)

ggsave(paste0("Figure/", "CCA_loading_sim", ".pdf"), width = 8.5, height = 5.5, dpi = 300)

