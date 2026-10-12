# Set working directory
rm(list=ls())
# set working directory to be where the current script is located
mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

# Input functions and packages
source("Function.R")
source("LCCA_function.R")

# Input data
load("Data/Protein.rda")
load("Data/Metabolomic.rda")

dat_X <- Protein
dat_Y <- Metabolomic
rm(Protein, Metabolomic)

n <- dat_X$n # Sample size
p_X <- dat_X$p # Dimension for the first dataset
p_Y <- dat_Y$p # Dimension for the second dataset

grid_num <- 201 # Number of evaluated time points for functions
time_grid <- seq(0, 1, length.out = grid_num)
max_var <- floor(c(min(p_X / 2, n), min(p_Y / 2, n)))

################################################################################
# CV correlation
cv_num <- 6
cv_mark <- lapply(1:cv_num, function(k){
  floor(seq(k, n, cv_num))
}) # CV mark

## FACD
error_FACD <- sapply(1:cv_num, function(k){
  
  dat_X_train <- dat_X 
  dat_Y_train <- dat_Y
  sam_mark <- (1:n)[-cv_mark[[k]]]
  
  dat_X_train$dat <- lapply(1:p_X, function(j){
    lapply(1:length(sam_mark), function(i){
      dat_X_train$dat[[j]][[sam_mark[[i]]]]
    })
  })
  dat_X_train$n <- length(sam_mark)
  dat_X_train$sub_label <- dat_X_train$sub_label[sam_mark]
  
  dat_Y_train$dat <- lapply(1:p_Y, function(j){
    lapply(1:length(sam_mark), function(i){
      dat_Y_train$dat[[j]][[sam_mark[[i]]]]
    })
  })
  dat_Y_train$n <- length(sam_mark)
  dat_Y_train$sub_label <- dat_Y_train$sub_label[sam_mark]
  
  Result_FACD_train <- FACD(grid_num = grid_num, 
                              dat_X, dat_Y,
                              max_var = max_var, # Maximum number of selected features
                              length_tun = 30, # Number of evaluated tuning parameters for two datasets
                              print_tuning = F,
                              comp = 1,
                              FVE = 0.95 # Fraction of variance explained
  )
  
  load_X <- Result_FACD_train$load_X[[1]]
  X_cca <- sapply(cv_mark[[k]], function(i){
    sum(sapply(1:p_X, function(l){
      mark <- sapply(1:length(dat_X$dat[[l]][[i]]$time), function(k) which.min(abs(dat_X$dat[[l]][[i]]$time[k] - time_grid)))
      return(mean(dat_X$dat[[l]][[i]]$val * load_X[l,mark]))
    }))
  })
  
  load_Y <- Result_FACD_train$load_Y[[1]]
  Y_cca <- sapply(cv_mark[[k]], function(i){
    sum(sapply(1:p_Y, function(l){
      mark <- sapply(1:length(dat_Y$dat[[l]][[i]]$time), function(k) which.min(abs(dat_Y$dat[[l]][[i]]$time[k] - time_grid)))
      return(mean(dat_Y$dat[[l]][[i]]$val * load_Y[l,mark]))
    }))
  })
  
  cor <- cor(X_cca, Y_cca)
  
  if(is.finite(cor) == F){
    cor <- 0
  }
  
  return(abs(cor))
})

### CV correlation
mean(error_FACD) 

################################################################################
## LCCA
error_LCCA <- sapply(1:cv_num, function(k){
  
  sam_mark <- (1:n)[-cv_mark[[k]]]
  dat_x <- NULL
  for(i in 1:length(sam_mark)){
    time_i <- unique(unlist(lapply(1:p_X, function(j){
      dat_X$dat[[j]][[sam_mark[[i]]]]$time
    })))
    
    mark <- sapply(1:length(time_i), function(k){
      prod(sapply(1:p_X, function(j){
        sum(dat_X$dat[[j]][[sam_mark[[i]]]]$time == time_i[k])
      }))
    })
    time_i <- time_i[mark == 1]
    
    dat_x <- rbind(dat_x, cbind(rep(sam_mark[[i]], length(time_i)),
                                time_i,
                                sapply(1:p_X, function(j){
                                  sapply(1:length(time_i), function(k){
                                    time_ob <- time_i[k]
                                    mark <- which(time_ob == dat_X$dat[[j]][[sam_mark[[i]]]]$time)
                                    dat_X$dat[[j]][[sam_mark[[i]]]]$val[mark]
                                  })
                                })
    ))
  }
  
  dat_y <- NULL
  for(i in 1:length(sam_mark)){
    time_i <- unique(unlist(lapply(1:p_Y, function(j){
      dat_Y$dat[[j]][[sam_mark[[i]]]]$time
    })))
    mark <- sapply(1:length(time_i), function(k){
      prod(sapply(1:p_Y, function(j){
        sum(dat_Y$dat[[j]][[sam_mark[[i]]]]$time == time_i[k])
      }))
    })
    time_i <- time_i[mark == 1]
    
    dat_y <- rbind(dat_y, cbind(rep(sam_mark[[i]], length(time_i)),
                                time_i,
                                sapply(1:p_Y, function(j){
                                  sapply(1:length(time_i), function(k){
                                    time_ob <- time_i[k]
                                    mark <- which(time_ob == dat_Y$dat[[j]][[sam_mark[[i]]]]$time)
                                    dat_Y$dat[[j]][[sam_mark[[i]]]]$val[mark]
                                  })
                                })
    ))
  }
  
  x <- list(X = as.matrix(t(dat_x[,-c(1:2)])),
            time = dat_x[,2],
            I = length(sam_mark),
            J = nrow(dat_x),
            visit = sapply(1:length(sam_mark), function(i) length(which(dat_x[,1] == sam_mark[[i]]))))
  
  y <- list(X = as.matrix(t(dat_y[,-c(1:2)])),
            time = dat_y[,2],
            I = length(sam_mark),
            J = nrow(dat_y),
            visit = sapply(1:length(sam_mark), function(i) length(which(dat_y[,1] == sam_mark[[i]]))))
  
  Result_LCCA <- lcca.linear_c(x = x, y = y, varthresh = 0.9)
  
  load_X <- matrix(Result_LCCA$re.X$phix0, nrow = p_X)[,1] +  matrix(Result_LCCA$re.X$phix1, nrow = p_X)[,1] %*% t(seq(0, 1, 1 / 200))
  X_cca <- sapply(cv_mark[[k]], function(i){
    sum(sapply(1:p_X, function(l){
      mark <- sapply(1:length(dat_X$dat[[l]][[i]]$time), function(k) which.min(abs(dat_X$dat[[l]][[i]]$time[k] - time_grid)))
      return(mean(dat_X$dat[[l]][[i]]$val * load_X[l,mark]))
    }))
  })
  
  load_Y <- matrix(Result_LCCA$re.Y$phix0, nrow = p_Y)[,1] +  matrix(Result_LCCA$re.Y$phix1, nrow = p_Y)[,1] %*% t(seq(0, 1, 1 / 200))
  Y_cca <- sapply(cv_mark[[k]], function(i){
    sum(sapply(1:p_Y, function(l){
      mark <- sapply(1:length(dat_Y$dat[[l]][[i]]$time), function(k) which.min(abs(dat_Y$dat[[l]][[i]]$time[k] - time_grid)))
      return(mean(dat_Y$dat[[l]][[i]]$val * load_Y[l,mark]))
    }))
  })
  
  cor <- cor(X_cca, Y_cca)
  
  if(is.finite(cor) == F){
    cor <- 0
  }
  return(abs(cor))
})

### CV correlation
mean(error_LCCA) 

################################################################################
## FSCCA
error_FSCCA <- sapply(1:cv_num, function(k){
  
  sam_mark <- (1:n)[-cv_mark[[k]]]
  
  norm_ploy <- sqrt(colMeans(poly(seq(0, 1, length.out = 201), degree = 2) ^ 2))
  dat_x_g <- t(sapply(sam_mark, function(i){
    c(sapply(1:p_X, function(j){
      (t(poly(dat_X$dat[[j]][[i]]$time, degree = 2)) / norm_ploy) %*% dat_X$dat[[j]][[i]]$val / length(dat_X$dat[[j]][[i]]$time)
    }))
  }))
  group_x_g <- rep(1:p_X, each = 2)
  
  dat_y_g <- t(sapply(sam_mark, function(i){
    c(sapply(1:p_Y, function(j){
      (t(poly(dat_Y$dat[[j]][[i]]$time, degree = 2)) / norm_ploy) %*% dat_Y$dat[[j]][[i]]$val / length(dat_Y$dat[[j]][[i]]$time)
    }))
  }))
  group_y_g <- rep(1:p_Y, each = 2)
  
  Result_FSCCA <- GSCCA(dat_x_g, dat_y_g, group_x_g, group_y_g, print_tuning = F, length_tun = 30, max_var = max_var)
  
  load_X <- matrix(Result_FSCCA$fit_opt$fit_cca$X_coef, ncol = 2, byrow = T) %*% (t(poly(seq(0, 1, 1 / 200), degree = 2)) / norm_ploy)
  X_cca <- sapply(cv_mark[[k]], function(i){
    sum(sapply(1:p_X, function(l){
      mark <- sapply(1:length(dat_X$dat[[l]][[i]]$time), function(k) which.min(abs(dat_X$dat[[l]][[i]]$time[k] - time_grid)))
      return(mean(dat_X$dat[[l]][[i]]$val * load_X[l,mark]))
    }))
  })
  
  load_Y <- matrix(Result_FSCCA$fit_opt$fit_cca$Y_coef, ncol = 2, byrow = T) %*% (t(poly(seq(0, 1, 1 / 200), degree = 2)) / norm_ploy)
  Y_cca <- sapply(cv_mark[[k]], function(i){
    sum(sapply(1:p_Y, function(l){
      mark <- sapply(1:length(dat_Y$dat[[l]][[i]]$time), function(k) which.min(abs(dat_Y$dat[[l]][[i]]$time[k] - time_grid)))
      return(mean(dat_Y$dat[[l]][[i]]$val * load_Y[l,mark]))
    }))
  })
  
  cor <- cor(X_cca, Y_cca)
  
  if(is.finite(cor) == F){
    cor <- 0
  }
  return(abs(cor))
})

### CV correlation
mean(error_FSCCA) 

################################################################################
error_MFPCA <- sapply(1:cv_num, function(k){
  
  sam_mark <- (1:n)[-cv_mark[[k]]]
  
  X <- multiFunData(lapply(1:p_X, function(j){
    funData(argvals = seq(0, 1, length.out = 201),
            X = t(sapply(sam_mark, function(i){
              obs <- rep(NA, 201)
              mark <- sapply(1:length(dat_X$dat[[j]][[i]]$time), function(k) which.min(abs(dat_X$dat[[j]][[i]]$time[k] - time_grid)))
              obs[mark] <- dat_X$dat[[j]][[i]]$val
              return(obs)
            })))
  }))
  MFPCA_X <- MFPCA(X, M = 1, uniExpansions = lapply(1:p_X, function(j) list(type = "uFPCA", pve = 0.9, nbasis = 5)))
  
  X <- multiFunData(lapply(1:p_Y, function(j){
    funData(argvals = seq(0, 1, length.out = 201),
            X = t(sapply(sam_mark, function(i){
              obs <- rep(NA, 201)
              mark <- sapply(1:length(dat_Y$dat[[j]][[i]]$time), function(k) which.min(abs(dat_Y$dat[[j]][[i]]$time[k] - time_grid)))
              obs[mark] <- dat_Y$dat[[j]][[i]]$val
              return(obs)
            })))
  }))
  MFPCA_Y <- MFPCA(X, M = 1, uniExpansions = lapply(1:p_Y, function(j) list(type = "uFPCA", pve = 0.9, nbasis = 5)))
  
  load_X <- t(sapply(1:p_X, function(j){
    MFPCA_X$functions[[j]]@X
  }))
  X_cca <- sapply(cv_mark[[k]], function(i){
    sum(sapply(1:p_X, function(l){
      mark <- sapply(1:length(dat_X$dat[[l]][[i]]$time), function(k) which.min(abs(dat_X$dat[[l]][[i]]$time[k] - time_grid)))
      return(mean(dat_X$dat[[l]][[i]]$val * load_X[l,mark]))
    }))
  })
  
  load_Y <- t(sapply(1:p_Y, function(j){
    MFPCA_Y$functions[[j]]@X
  }))
  Y_cca <- sapply(cv_mark[[k]], function(i){
    sum(sapply(1:p_Y, function(l){
      mark <- sapply(1:length(dat_Y$dat[[l]][[i]]$time), function(k) which.min(abs(dat_Y$dat[[l]][[i]]$time[k] - time_grid)))
      return(mean(dat_Y$dat[[l]][[i]]$val * load_Y[l,mark]))
    }))
  })
  
  cor <- cor(X_cca, Y_cca)
  
  if(is.finite(cor) == F){
    cor <- 0
  }
  return(abs(cor))
})

### CV correlation
mean(error_MFPCA) 

################################################################################
# Plot
dat_plot <- rbind(data.frame(Cor = error_FACD, Method = rep("FACD", length(error_FACD))),
                  data.frame(Cor = error_LCCA, Method = rep("LCCA", length(error_LCCA))),
                  data.frame(Cor = error_FSCCA, Method = rep("FSCCA", length(error_FACD))),
                  data.frame(Cor = error_MFPCA, Method = rep("MFPCA", length(error_FACD)))
)
dat_plot$Method <- factor(dat_plot$Method, levels = c("LCCA", "FSCCA", "MFPCA", "FACD"))

save(dat_plot, file = "Result/CV_result.rda", version = 2)

p_1 <- ggplot(dat_plot) + 
  geom_boxplot(aes(x = Method, y = Cor), outliers = F, width = 0.65) +
  theme_bw(base_family = "Times") +
  # scale_fill_manual(values = c("blue", "orange", "purple", "red")) +
  xlab("Method") + ylab("Testing correlation") + 
  theme(legend.position = "none",
        text = element_text(size = 15),
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  labs(title = "(B)", fill = "")

# ggsave(paste0("Figure/", "correlation", ".pdf"), width = 4.5, height = 3, dpi = 300)

################################################################################
# Correlation of scores
file_name <- c("fit_Protein_Metabolomic.rda", "fit_Protein_Lipid.rda",
               "fit_Metabolomic_Lipid.rda")

var_name <- cbind(var1 = c("Protein", "Protein", "Metabolite"),
                  var2 = c("Metabolite", "Lipid", "Lipid"))

dat_plot <- data_frame()
cor <- vector()
for(i in 1:3){
  load(paste0("Result/", file_name[i]))
  X_cca <- Result$Z_X[[1]]
  Y_cca <- Result$Z_Y[[1]]
  cor[i] <- abs(crossprod(X_cca, Y_cca) / sqrt(crossprod(X_cca, X_cca) * crossprod(Y_cca, Y_cca)))
  dat_plot <- rbind(dat_plot, 
                    data.frame(
                      X = X_cca, Y = Y_cca,
                      cor = c(cor[i], rep(NA, n - 1)), var1 = rep(var_name[i,1], n), var2 = rep(var_name[i,2], n),
                      text_x = c(min(X_cca), rep(NA, n - 1)), text_y = c(max(Y_cca) + 0.2, rep(NA, n - 1))
                    ))
  
}

dat_plot$var1 <- factor(dat_plot$var1, levels = c("Protein", "Metabolite", "Lipid"))
dat_plot$var2 <- factor(dat_plot$var2, levels = c("Protein", "Metabolite", "Lipid"))

p_2 <- ggplot(dat_plot) + 
  geom_point(aes(x = X, y = Y)) +  
  ggrepel::geom_text_repel(aes(x = text_x, y = text_y, label = paste0("Correlation ", round(cor, 3))), size = 4, max.overlaps = 100) +
  facet_wrap(~ var1 + var2, scales = "free") +
  theme_bw(base_family = "Times") +
  xlab("The first dataset") + ylab("The second dataset") + 
  theme(legend.position = "none",
        text = element_text(size = 15),
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  labs(title = "(A)", fill = "")

# Combined figure
gridExtra::grid.arrange(p_2, p_1, ncol = 2, layout_matrix = t(c(1, 1, 1, 1, 1, 2, 2, 2)))
p <- gridExtra::grid.arrange(p_2, p_1, ncol = 2, layout_matrix = t(c(1, 1, 1, 1, 1, 2, 2, 2)))

ggsave(paste0("Figure/", "correlation", ".pdf"), p, width = 15, height = 4, dpi = 300)
