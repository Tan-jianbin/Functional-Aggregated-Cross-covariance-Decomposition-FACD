# Package
library(fdapace)
library(TruncatedNormal)
library(mvtnorm)
library(fields)
library(filling)
# library(knitr)
library(ggplot2)
library(mgcv)
library(sparsegl)
library(dplyr)
library(expm)
library(fda)
library(CCA)
library(gss)
library(MFPCA)
library(funData)
library(face)

###############################################################################################
# Wrapping function of pre-processing
wrap.dat <- function(dat, dat_count, time_label){
  
  feature_label <- colnames(dat)
  sub_label <- unique(dat_count$SubjectID)
  
  n <- length(sub_label) # Sample size
  p <- length(feature_label) # Dimension
  
  fda_object <- lapply(1:p, function(i){
    lapply(1:n, function(j){
      dat_ij <- data.frame(time = dat_count$timepoint[which(dat_count$SubjectID == sub_label[j])],
                           val = dat[which(dat_count$SubjectID == sub_label[j]), which(colnames(dat) == feature_label[i])])
      time <- sort(unique(dat_ij$time))
      
      val <- sapply(1:length(time), function(k) mean(dat_ij$val[which(dat_ij$time == time[k])]))
      time <- sapply(1:length(time), function(k) time_label[time[k]])
      dat_ij <- data.frame(time = time, val= val)[is.finite(val),]
      dat_ij$val <- dat_ij$val
      return(dat_ij)
    })
  })
  
  ### Scale and time transformation
  x_grid <- sort(unique(unlist(lapply(1:p, function(i) unlist(sapply(1:n, function(j)  fda_object[[i]][[j]]$time))))))
  
  min_x <- min(x_grid)
  max_x <- max(x_grid)
  x_grid <- (x_grid - min_x) / (max_x - min_x)
  
  for(i in 1:p){
    mean_fac <- mean(unlist(sapply(1:n, function(j)  fda_object[[i]][[j]]$val)))
    scale_fac <- sqrt(mean((unlist(sapply(1:n, function(j)  fda_object[[i]][[j]]$val))) ^ 2))
    for(j in 1:n){
      if(length(fda_object[[i]][[j]]$time) != 0){
        fda_object[[i]][[j]]$time <- (fda_object[[i]][[j]]$time - min_x) /  (max_x - min_x)
        fda_object[[i]][[j]]$val <- (fda_object[[i]][[j]]$val - mean_fac) / scale_fac
      }
    }
  }
  
  return(list(dat = fda_object, n = n, p = p,
              sub_label = sub_label, feature_label = feature_label,
              time_grid = x_grid, time_min = min_x, time_max = max_x))
}

# Data generation
data_gen <- function(n, rho, time_num, p_X, p_Y, comp, p_ture){
  
  var_comp <- sapply(1:comp, function(k) 0.5 * k ^ (-2))
  mat <- lapply(1:comp, function(k){
    matrix(c(1, rho / k ^ 2, rho / k ^ 2, 1), 2, 2)
  })
  val <- lapply(1:comp, function(k){
    t(mvtnorm::rmvnorm(n, mean = rep(0, 2), sigma = mat[[k]] * var_comp[k]))
  }) 
  
  alpha_X <- sapply(1:comp, function(k){
    A <- sapply(1:3, function(g) c(runif(p_ture, 1, 2 * (5 - g)) * sample(c(1,-1), p_ture, replace = T), rep(0, p_X - p_ture)))
    A <- A / sqrt(sum(A ^ 2))
    return(A)
  })
  alpha_X <- qr.Q(qr(alpha_X))
  alpha_X <- lapply(1:comp, function(k) matrix(alpha_X[,k], ncol = 3))
  
  alpha_Y <- sapply(1:comp, function(k){
    A <- sapply(1:3, function(g) c(runif(p_ture, 1, 2 * (5 - g)) * sample(c(1,-1), p_ture, replace = T), rep(0, p_X - p_ture)))
    A <- A / sqrt(sum(A ^ 2))
    return(A)
  })
  alpha_Y <- qr.Q(qr(alpha_Y))
  alpha_Y <- lapply(1:comp, function(k) matrix(alpha_Y[,k], ncol = 3))
  
  # The first dataset
  time_tol <- seq(0, 1, length.out = 201)
  L <- lapply(1:n, function(i){
    time <- sample((time_num - 3):(time_num), 1)
    time_point <- sort(sample(time_tol, time, replace = F))
    basis <- fourier(time_point, nbasis =  4, period = 1)[,c(2, 3, 4)]
    Ly <- t(sapply(1:p_X, function(j)  rowSums(sapply(1:comp, function(k) val[[k]][1,i] * basis %*% alpha_X[[k]][j,]))
    ))
    return(list(Lt = time_point, Ly = Ly))
  })
  
  load_X <- lapply(1:comp, function(k){
    t(sapply(1:p_X, function(j) fourier(seq(0, 1, 1 / (201 - 1)), nbasis = 4, period = 1)[,c(2, 3, 4)] %*% alpha_X[[k]][j,]))
  })
  load_X <- lapply(1:comp, function(k){
    load_X[[k]] / sqrt(sum(rowMeans(load_X[[k]] ^ 2)))
  })
  
  sigma <- sqrt(sum(sapply(1:p_X, function(j) mean(sapply(1:n, function(i) mean(L[[i]]$Ly[j,] ^ 2))))) / p_ture)
  dist_X <- rdist(1:p_X, 1:p_X) 
  matrix_X <- exp(- (dist_X) ^ 2 / 10) * (sigma) ^ 2  * 0.5 / 3
  
  L_noise <- lapply(1:n, function(i){
    mvtnorm::rmvnorm(3, mean = rep(0, p_X), sigma = matrix_X)
  })
  
  dat_X <- lapply(1:p_X, function(j){
    lapply(1:n, function(i){
      list(
        time = L[[i]]$Lt,
        val = L[[i]]$Ly[j,] + fourier(L[[i]]$Lt, nbasis =  4, period = 1)[,c(2, 3, 4)] %*% L_noise[[i]][,j] + rnorm(length(L[[i]]$Ly[j,]), 0,  sqrt(sigma) * 0.1),
        sigma = sigma
      )
    })
  })
  
  ## The second dataset
  time_tol <- seq(0, 1, length.out = 201)
  L <- lapply(1:n, function(i){
    time <- sample((time_num - 3):(time_num), 1)
    time_point <- sort(sample(time_tol, time, replace = F))
    # time_point <- sort(sample(seq(0, 1, length.out = 20), time, replace = F))
    basis <- fourier(time_point, nbasis =  10, period = 1)[,c(2,3,4)]
    Ly <- t(sapply(1:p_Y, function(j)  rowSums(sapply(1:comp, function(k) val[[k]][2,i] * basis %*% alpha_Y[[k]][j,]))))
    return(list(Lt = time_point, Ly = Ly))
  })
  
  load_Y <- lapply(1:comp, function(k){
    t(sapply(1:p_Y, function(j) fourier(seq(0, 1, 1 / (201 - 1)), nbasis = 10, period = 1)[,c(2,3,4)] %*% alpha_Y[[k]][j,]))
  })
  load_Y <- lapply(1:comp, function(k){
    load_Y[[k]] / sqrt(sum(rowMeans(load_Y[[k]] ^ 2)))
  })
  
  sigma <- sqrt(sum(sapply(1:p_Y, function(j) mean(unlist(sapply(1:n, function(i) L[[i]]$Ly[j,] ^ 2)))) / p_ture))
  dist_Y <- rdist(1:p_Y, 1:p_Y) 
  matrix_Y <- exp(- (dist_Y) ^ 2 / 10) * (sigma) ^ 2  * 0.5 / 3
  
  L_noise <- lapply(1:n, function(i){
    mvtnorm::rmvnorm(3, mean = rep(0, p_Y), sigma = matrix_Y)
  })
  
  dat_Y <- lapply(1:p_Y, function(j){
    lapply(1:n, function(i){
      list(
        time = L[[i]]$Lt,
        val = L[[i]]$Ly[j,] + fourier(L[[i]]$Lt, nbasis =  4, period = 1)[,c(2, 3, 4)] %*% L_noise[[i]][,j] + rnorm(length(L[[i]]$Ly[j,]), 0,  sqrt(sigma) * 0.1),
        sigma = sigma
      )
    })
  })
  
  return(list(dat_X = dat_X, dat_Y = dat_Y, n = n, 
              p_X = p_X, p_Y = p_Y, val = val, alpha_X = alpha_X, alpha_Y = alpha_Y,
              load_X = load_X, load_Y = load_Y))
}

# Function
fit_cov_surface <- function(M, time_mark_X, time_mark_Y, time){
  
  n <- length(M)
  matrix_sample <- vector()
  for(i in 1:n){
    A <- sapply(1:length(time_mark_X[[i]]), function(k){
      cbind(rep(time[time_mark_X[[i]]][k], length(time_mark_Y[[i]])),
            time[time_mark_Y[[i]]],
            M[[i]][k,]
      )
    }, simplify = "array")
    
    matrix_sample <- rbind(matrix_sample, data.frame(
      x = c(A[,1,]),
      y = c(A[,2,]),
      z = c(A[,3,]),
      w = sapply(1:length(A[,1,]), function(k) 1 / (length(time_mark_X[[i]]) * length(time_mark_Y[[i]])))
    ))
  }
  
  uni_time_x <- sort(unique(matrix_sample$x))
  uni_time_y <- sort(unique(matrix_sample$y))
  num_k <- min(c(length(uni_time_x), length(uni_time_y), 15))
  
  x <- rep(time, length(time))
  y <- c(sapply(time, function(k) rep(k, length(time))))
  
  if(length(uni_time_x) * length(uni_time_x) <= 25){
    dat_tran <- sapply(1:length(uni_time_x), function(k){
      sapply(1:length(uni_time_y), function(j){
        obs <- matrix_sample$z[intersect(which(matrix_sample$x == uni_time_x[k]), which(matrix_sample$y == uni_time_y[j]))]
        return(c(uni_time_x[k], uni_time_y[j], mean(obs), length(obs)))
      }, simplify = "array")
    }, simplify = "array")
    
    matrix_sample <- data.frame(
      x = c(dat_tran[1,,]),
      y = c(dat_tran[2,,]),
      z = c(dat_tran[3,,])
    )
    rm(dat_tran)
    
    fit <- Tps(x = cbind(matrix_sample$x, matrix_sample$y), Y = matrix_sample$z, df = nrow(matrix_sample) - 2, method = "GCV.model")
    
    Cov_fit <- c(predictSurface(fit, grid.list = list(x = time, y = time),
                                extrap = T)$z)
    
  }else{
    
    fit <- gam(z ~ te(x, y, bs = "cr", k = c(num_k, num_k)), data = matrix_sample, weights = matrix_sample$w, family = "gaussian", method="GCV.Cp")
    Cov_fit <- predict(fit, newdata = list(x = x, y = y))
    
  }
  rm(fit)
  
  return(matrix(Cov_fit, nrow = length(time)))
}

fit_basis <- function(M, time_mark, time){
  
  fit_cov <- fit_cov_surface(M, time_mark, time_mark, time)
  fit <- eigen(fit_cov, symmetric = T)
  eigenfunction <- fit$vectors[,1] * sqrt(nrow(fit_cov))
  
  return(
    list(
      eigenvalue = c(t(eigenfunction) %*% fit_cov %*% eigenfunction / nrow(fit_cov) ^ 2),
      eigenfunction = eigenfunction,
      eigendocomposition = fit,
      tol_nor = sum(fit[[1]]) / nrow(fit_cov)
    )
  )
}

FSVD_r <- function(M, time, time_mark, R_max, FVE){
  
  fit <- fit_basis(M, time_mark, time)
  tol_norm <- fit$tol_nor
  op_basis <- list()
  
  op_basis$phi <- fit$eigendocomposition$vectors[,1:R_max] * sqrt(length(time))
  op_basis$workGrid <- time
  op_basis$eigenvalue <- fit$eigendocomposition$values[1:R_max] / length(time)
  
  op_basis$comp_num <- max(sum(cumsum(op_basis$eigenvalue) / tol_norm <= FVE)  + 1, 2)
  
  return(op_basis)
}

Optimal_basis_est <- function(dat_X, dat_Y, grid_num, max_comp_num = 5, mean_zero = F, FVE){
  
  n <- dat_X$n
  p_X <- dat_X$p
  p_Y <- dat_Y$p
  
  fda_object_X <- dat_X$dat
  fda_object_Y <- dat_Y$dat
  
  # Estimation of mean functions
  if(mean_zero == T){
    Fit_mean_X <- lapply(1:p_X, function(i){
      list(mu = rep(0, grid_num), workGrid = seq(0, 1, length.out = grid_num))
    })
    
    Fit_mean_Y <- lapply(1:p_Y, function(i){
      list(mu = rep(0, grid_num), workGrid = seq(0, 1, length.out = grid_num))
    })
  }else{
    
    Lt <- unique(unlist((lapply(1:n, function(j) fda_object_X[[1]][[j]]$time))))
    
    if(length(Lt) <= 5){
      Fit_mean_X <- lapply(1:p_X, function(i, fda_object, grid_num, n){
        Lt <- unique(unlist((lapply(1:n, function(j) fda_object[[1]][[j]]$time))))
        Ly <- sapply(1:n, function(j) sapply(1:length(Lt), function(k){
          mark <- which(fda_object[[i]][[j]]$time == Lt[k])
          if(length(mark) > 0){
            fda_object[[i]][[j]]$val[mark]
          }else{
            NA
          }
        }))
        
        dis_mu <- sapply(1:length(Lt), function(t) mean(Ly[t,], na.rm = T))
        mu <- spline(x = Lt, y = dis_mu, xout = seq(0, 1, length.out = grid_num))$y
        
        return(list(mu = mu, workGrid = seq(0, 1, length.out = grid_num)))
      }, fda_object = fda_object_X, grid_num = grid_num, n = n)
      
    }else{
      
      Fit_mean_X <- lapply(1:p_X, function(i, fda_object, grid_num, n){
        Ly <- sapply(1:n, function(j) fda_object[[i]][[j]]$val)
        Lt <- sapply(1:n, function(j) fda_object[[i]][[j]]$time)
        
        mu <- GetMeanCurve(Ly, Lt, optns = list(error = T, nRegGrid = grid_num, 
                                                methodBwMu = "GCV",
                                                methodMuCovEst = "smooth"
        ))$mu
        
        return(list(mu = mu, workGrid = seq(0, 1, length.out = grid_num)))
      }, fda_object = fda_object_X, grid_num = grid_num, n = n)
    }
    
    Lt <- unique(unlist((lapply(1:n, function(j) fda_object_Y[[1]][[j]]$time))))
    
    if(length(Lt) <= 5){
      Fit_mean_Y <- lapply(1:p_Y, function(i, fda_object, grid_num, n){
        Lt <- unique(unlist((lapply(1:n, function(j) fda_object[[i]][[j]]$time))))
        Ly <- sapply(1:n, function(j) sapply(1:length(Lt), function(k){
          mark <- which(fda_object[[i]][[j]]$time == Lt[k])
          if(length(mark) > 0){
            fda_object[[i]][[j]]$val[mark]
          }else{
            NA
          }
        }))
        
        dis_mu <- sapply(1:length(Lt), function(t) mean(Ly[t,], na.rm = T))
        mu <- spline(x = Lt, y = dis_mu, xout = seq(0, 1, length.out = grid_num))$y
        
        return(list(mu = mu, workGrid = seq(0, 1, length.out = grid_num)))
      }, fda_object = fda_object_Y, grid_num = grid_num, n = n)
      
    }else{
      Fit_mean_Y <- lapply(1:p_Y, function(i, fda_object, grid_num, n){
        Ly <- sapply(1:n, function(j) fda_object[[i]][[j]]$val)
        Lt <- sapply(1:n, function(j) fda_object[[i]][[j]]$time)
        
        mu <- GetMeanCurve(Ly, Lt, optns = list(error = T, nRegGrid = grid_num, 
                                                methodBwMu = "GCV",
                                                methodMuCovEst = "smooth"
        ))$mu
        
        return(list(mu = mu, workGrid = seq(0, 1, length.out = grid_num)))
      }, fda_object = fda_object_Y, grid_num = grid_num, n = n)
    }
  }
  
  time_mark_X <- lapply(1:n, function(i){
    sapply(1:length(fda_object_X[[1]][[i]]$time), function(h){
      which.min(abs(fda_object_X[[1]][[i]]$time[h] - Fit_mean_X[[1]]$workGrid))
    })
  }) 
  time_mark_Y <- lapply(1:n, function(i){
    sapply(1:length(fda_object_Y[[1]][[i]]$time), function(h){
      which.min(abs(fda_object_Y[[1]][[i]]$time[h] - Fit_mean_Y[[1]]$workGrid))
    })
  }) 
  
  # Optimal basis estimation
  G_X <- list()
  G_Y <- list()
  
  time <- seq(0, 1, length.out = grid_num)
  
  for(i in 1:n){
    G <- sapply(1:p_X, function(j_1){
      sapply(1:p_Y, function(j_2){
        tcrossprod(fda_object_X[[j_1]][[i]]$val - Fit_mean_X[[j_1]]$mu[time_mark_X[[i]]],
                   fda_object_Y[[j_2]][[i]]$val - Fit_mean_Y[[j_2]]$mu[time_mark_Y[[i]]])
      }, simplify = "array")
    }, simplify = "array")
    
    G_X[[i]] <- sapply(1:length(time_mark_X[[i]]), function(k){
      sapply(1:length(time_mark_X[[i]]), function(g){
        mean(G[k,,,] * G[g,,,])
      })
    }, simplify = "array")
    
    G_Y[[i]] <- sapply(1:length(time_mark_Y[[i]]), function(k){
      sapply(1:length(time_mark_Y[[i]]), function(g){
        mean(G[,k,,] * G[,g,,])
      })
    }, simplify = "array")
    
    rm(G)
    gc()
  }
  
  ## X data
  op_basis_X <- FSVD_r(M = G_X, time, time_mark = time_mark_X, R_max = max_comp_num, FVE)
  
  ## Y data
  op_basis_Y <- FSVD_r(M = G_Y, time, time_mark = time_mark_Y, R_max = max_comp_num, FVE)
  
  return(list(
    Mean = list(X = Fit_mean_X, Y = Fit_mean_Y),
    Time_mark = list(X = time_mark_X, Y = time_mark_Y),
    Optimal_bais = list(X = op_basis_X, Y = op_basis_Y))
  )
}


cca_AX <- function(lambda, abs, Initial_cca, n, group_x, group_y){
  
  a_X <- Initial_cca$X_coef
  a_Y <- Initial_cca$Y_coef
  
  error <- 1
  ite <- 0
  
  num_X <- length(unique(group_x))
  num_Y <- length(unique(group_y))
  
  while((error > abs) & (ite < 200)){
    ite <- ite + 1
    
    a_X_t <- a_X
    a_Y_t <- a_Y
    
    a_X <- c(Initial_cca$Cov %*% a_Y) 
    a_X <- c(sapply(1:num_X, function(j){
      norm <- sqrt(sum(a_X[group_x == j] ^ 2))
      if(norm == 0){
        rep(0, sum(group_x == j))
      }else{
        max((1 - lambda[1] / norm), 0) * a_X[group_x == j]
      }
    })) 
    x_norm <- c(sqrt(crossprod(a_X)))
    if(x_norm != 0){
      a_X <- a_X / x_norm
    }
    
    a_Y <- c(t(a_X) %*% Initial_cca$Cov) 
    a_Y <- c(sapply(1:num_Y, function(j){
      norm <- sqrt(sum(a_Y[group_y == j] ^ 2))
      if(norm == 0){
        rep(0, sum(group_y == j))
      }else{
        max((1 - lambda[2] / norm), 0) * a_Y[group_y == j]
      }
    })) 
    y_norm <- c(sqrt(crossprod(a_Y)))
    if(y_norm != 0){
      a_Y <- a_Y / y_norm
    }
    
    if(x_norm != 0){
      error_x <- (sum(abs(a_X_t - a_X)) / sum(abs(a_X_t))) 
    }else{
      error_x <- 0
    }
    
    if(y_norm != 0){
      error_y <- (sum(abs(a_Y_t - a_Y)) / sum(abs(a_Y_t)))
    }else{
      error_y <- 0
    }
    
    error <- error_x + error_y
  }
  
  
  return(list(X_coef = as.numeric(a_X),
              Y_coef = as.numeric(a_Y)
  ))
}

cv_cca <- function(lambda, Initial_cca, n, zx, zy, group_x, group_y, max_var, cv_num = min(5, n)){
  
  cv_mark <- lapply(1:cv_num, function(k){
    floor(seq(k, n, cv_num))
  })
  
  fit_cca <- cca_AX(lambda, abs = 0.001, Initial_cca, n, group_x, group_y)
  num_X <- length(which(sapply(unique(group_x), function(j) sum(fit_cca$X_coef[group_x == j] ^ 2)) != 0))
  num_Y <- length(which(sapply(unique(group_y), function(j) sum(fit_cca$Y_coef[group_x == j] ^ 2)) != 0))
  
  if((num_X > max_var[1]) | (num_X <= 0) |
     (num_Y > max_var[2]) | (num_Y <= 0)
  ){
    error <- 0
  }else{
    error <- mean(sapply(1:cv_num, function(k){
      cca_cv_fit <- cca_AX(lambda, abs = 0.001, 
                           Initial_cca = list(X_coef = fit_cca$X_coef, Y_coef = fit_cca$Y_coef, 
                                              Cov = t(zx[-cv_mark[[k]],]) %*% zy[-cv_mark[[k]],] / (n - length(cv_mark[[k]])),
                                              var_X = rep(1, length(fit_cca$X_coef)), var_Y = rep(1, length(fit_cca$Y_coef))),
                           n, group_x, group_y)
      
      X_cca <- zx[cv_mark[[k]],] %*% cca_cv_fit$X_coef
      Y_cca <- zy[cv_mark[[k]],] %*% cca_cv_fit$Y_coef
      
      cor <- as.numeric(crossprod(X_cca, Y_cca) / sqrt(crossprod(X_cca, X_cca) * crossprod(Y_cca, Y_cca)))
      
      if(is.finite(cor) == F){
        error <- 0
      }else{
        error <- abs(cor)
      }
      return(error)
    }))
  }
  
  if(is.finite(error) == F){
    error <- 0
  }
  
  return(list(error = error, 
              fit_cca = fit_cca,
              lambda = lambda
  ))
}

cv_slcca <- function(lambda, Initial_cca, n, group_x, group_y, max_var){
  
  fit_cca <- cca_AX(lambda, abs = 0.001, Initial_cca, n, group_x, group_y)
  num_X <- length(which(sapply(unique(group_x), function(j) sum(fit_cca$X_coef[group_x == j] ^ 2)) != 0))
  num_Y <- length(which(sapply(unique(group_y), function(j) sum(fit_cca$Y_coef[group_y == j] ^ 2)) != 0))
  
  if((num_X > max_var[1]) | (num_Y > max_var[2])){
    error <- - 10 ^ 11
  }else if((num_X == 0) | (num_Y == 0)){
    error <- - 10 ^ 10
  }else{
    error <- mean(sapply(1:Initial_cca$cv_num, function(k){
      cca_cv_fit <- cca_AX(lambda, abs = 0.001, 
                           Initial_cca = list(X_coef = fit_cca$X_coef, Y_coef = fit_cca$Y_coef, 
                                              Cov = Initial_cca$CV[[k]]$train), n, 
                           group_x, group_y)
      
      cor <- c(t(cca_cv_fit$X_coef) %*% Initial_cca$CV[[k]]$test %*% cca_cv_fit$Y_coef)
      
      return(cor)
    }))
  }
  
  return(list(error = error, 
              fit_cca = fit_cca,
              lambda = lambda
  ))
}

GSCCA <- function(dat_x, dat_y, group_x, group_y, print_tuning = T, length_tun = 30, max_var,
                  dmean = T){
  
  n <- nrow(dat_x)
  px <- ncol(dat_x)
  py <- ncol(dat_y)
  
  if(dmean == T){
    zx <- t(t(dat_x) - colMeans(dat_x))
    zy <- t(t(dat_y) - colMeans(dat_y))
  }else{
    zx <- dat_x
    zy <- dat_y
  }
  
  Cov <- t(zx) %*% (zy) / n
  fit <- svd(Cov)
  Initial_cca <- list(X_coef = fit$u[,1], Y_coef = fit$v[,1], Cov = Cov, var_X = rep(1, length(fit$u[,1])), var_Y = length(fit$v[,1]))
  rm(Cov)
  
  # Use grid search to determine the initial tuning parameters
  max_t <- fit$d[1] * max(sapply(unique(group_x), function(j) sqrt(sum(fit$u[which(group_x == j),1] ^ 2)))) 
  
  fit <- cv_cca(lambda = c(max_t, max_t), Initial_cca, n, zx, zy, group_x, group_y, max_var)
  var_num_X <- length(which(sapply(unique(group_x), function(j) sum(fit$fit_cca$X_coef[group_x == j] ^ 2)) != 0))
  var_num_Y <- length(which(sapply(unique(group_y), function(j) sum(fit$fit_cca$Y_coef[group_y == j] ^ 2)) != 0))
  
  while(var_num_X + var_num_Y != 0){
    max_t <- max_t * 5
    fit <- cv_cca(lambda = c(max_t, max_t), Initial_cca, n, zx, zy, group_x, group_y, max_var)
    var_num_X <- length(which(sapply(unique(group_x), function(j) sum(fit$fit_cca$X_coef[group_x == j] ^ 2)) != 0))
    var_num_Y <- length(which(sapply(unique(group_y), function(j) sum(fit$fit_cca$Y_coef[group_y == j] ^ 2)) != 0))
  }
  
  lambda_list <- c(max_t / 10, max_t / 100, max_t / 1000, max_t / 10000)
  var_num_X <- matrix(NA, length(lambda_list), length(lambda_list))
  var_num_Y <- matrix(NA, length(lambda_list), length(lambda_list))
  
  for(a in 1:length(lambda_list)){
    for(b in 1:length(lambda_list)){
      fit <- cv_cca(lambda = c(lambda_list[a], lambda_list[b]), Initial_cca, n, zx, zy, group_x, group_y, max_var)
      var_num_X[b,a] <- length(which(sapply(unique(group_x), function(j) sum(fit$fit_cca$X_coef[group_x == j] ^ 2)) != 0))
      var_num_Y[b,a] <- length(which(sapply(unique(group_y), function(j) sum(fit$fit_cca$Y_coef[group_y == j] ^ 2)) != 0))
    }
  }
  
  A <- var_num_X > max_var[1]
  if(length(which(colMeans(A) > 0)) == 0){
    min_X <- max_t / 100000
    fit <- cv_cca(lambda = c(min_X, max_t), Initial_cca, n, zx, zy, group_x, group_y, max_var)
    while(length(which(sapply(unique(group_x), function(j) sum(fit$fit_cca$X_coef[group_x == j] ^ 2)) != 0)) <= max_var[1]){
      min_X <- min_X / 10
      fit <- cv_cca(lambda = c(min_X, max_t), Initial_cca, n, zx, zy, group_x, group_y, max_var)
    }
  }else{
    min_X <- lambda_list[min(which(colMeans(A) > 0))]
  }
  
  A <- var_num_Y > max_var[2]
  if(length(which(rowMeans(A) > 0)) == 0){
    min_Y <- max_t / 100000
    fit <- cv_cca(lambda = c(min_X, min_Y), Initial_cca, n, zx, zy, group_x, group_y, max_var)
    while(length(which(sapply(unique(group_y), function(j) sum(fit$fit_cca$Y_coef[group_y == j] ^ 2)) != 0)) <= max_var[2]){
      min_Y <- min_Y / 10
      fit <- cv_cca(lambda = c(min_X, min_Y), Initial_cca, n, zx, zy, group_x, group_y, max_var)
    }
  }else{
    min_Y <- lambda_list[min(which(rowMeans(A) > 0))]
  }
  
  lambda_list_X <- exp(seq(log(max_t), log(min_X), length.out = length_tun + 1))[-1] 
  lambda_list_Y <- matrix(NA, length_tun, length_tun)
  CV_error <- matrix(NA, length(lambda_list_X), length(lambda_list_X))
  
  for(a in 1:length(lambda_list_X)){
    a_X <- c(Initial_cca$Cov %*% Initial_cca$Y_coef) 
    a_X <- c(sapply(1:p_X, function(j){
      norm <- sqrt(sum(a_X[group_x == j] ^ 2))
      if(norm == 0){
        rep(0, sum(group_x == j))
      }else{
        max((1 - lambda_list_X[a] / norm), 0) * a_X[group_x == j]
      }
    }))
    x_norm <- c(sqrt(crossprod(a_X)))
    if(x_norm != 0){
      a_X <- a_X / x_norm
    }
    
    a_Y <- c(t(a_X) %*% Initial_cca$Cov) 
    max_Y <- max(c(sapply(1:p_Y, function(j){
      sqrt(sum(a_Y[group_y == j] ^ 2))
    })))
    if(max_Y == 0){
      max_Y <- max_t
    }
    
    lambda_list_Y[,a] <- exp(seq(log(max_Y), log(min_Y), length.out = length_tun + 1))[-1] 
    
    for(b in 1:length(lambda_list_Y[,a])){
      fit <- cv_cca(lambda = c(lambda_list_X[a], lambda_list_Y[b,a]), Initial_cca, n, zx, zy, group_x, group_y, max_var)
      CV_error[b,a] <- fit$error
      if(CV_error[b,a] == max(CV_error, na.rm = T)){
        lambda <- c(lambda_list_X[a], lambda_list_Y[b,a])
      }
      if(print_tuning == T){
        print(CV_error)
      }
    }
  }
  
  fit <- optim(par = lambda, fn = function(lambda){
    fit <- cv_cca(lambda, Initial_cca, n, zx, zy, group_x, group_y, max_var)
    return(- fit$error)
  }, method = "L-BFGS-B", lower = c(min(lambda_list_X), min(lambda_list_Y)), 
  upper = c(max(lambda_list_X), max(lambda_list_Y)))
  
  fit_opt <- cv_cca(lambda = fit$par, Initial_cca, n, zx, zy, group_x, group_y, max_var)
  cor <- fit_opt$error
  
  return(list(fit_opt = fit_opt,
              error = cor,
              CV_error = CV_error,
              X_cca = zx %*% fit_opt$fit_cca$X_coef,
              Y_cca = zy %*% fit_opt$fit_cca$Y_coef,
              var_X = which(sapply(unique(group_x), function(j) sum(fit_opt$fit_cca$X_coef[group_x == j] ^ 2)) != 0),
              var_Y = which(sapply(unique(group_y), function(j) sum(fit_opt$fit_cca$Y_coef[group_y == j] ^ 2)) != 0)
  ))
}

SLCCA_tr <- function(group_x, group_y, print_tuning = F, length_tun = 30, max_var, Initial_cca, n, p_X, p_Y){
  
  # Use grid search to determine the initial tuning parameters
  vec <- c(Initial_cca$Cov %*% Initial_cca$Y_coef) 
  max_t <- max(sapply(unique(group_x), function(j) sqrt(sum(vec[which(group_x == j)] ^ 2)))) 
  
  fit <- cca_AX(lambda = c(max_t, max_t), abs = 0.001, Initial_cca, n, group_x, group_y)
  var_num_X <- length(which(sapply(unique(group_x), function(j) sum(fit$X_coef[group_x == j] ^ 2)) != 0))
  var_num_Y <- length(which(sapply(unique(group_y), function(j) sum(fit$Y_coef[group_y == j] ^ 2)) != 0))
  
  while(var_num_X + var_num_Y != 0){
    max_t <- max_t * 5
    fit <- cca_AX(lambda = c(max_t, max_t), abs = 0.001, Initial_cca, n, group_x, group_y)
    var_num_X <- length(which(sapply(unique(group_x), function(j) sum(fit$fit_cca$X_coef[group_x == j] ^ 2)) != 0))
    var_num_Y <- length(which(sapply(unique(group_y), function(j) sum(fit$fit_cca$Y_coef[group_y == j] ^ 2)) != 0))
  }
  
  lambda_list <- c(max_t / 2, max_t / 10, max_t / 100, max_t / 1000, max_t / 10000)
  var_num_X <- matrix(NA, length(lambda_list), length(lambda_list))
  var_num_Y <- matrix(NA, length(lambda_list), length(lambda_list))
  
  for(a in 1:length(lambda_list)){
    for(b in 1:length(lambda_list)){
      fit <- cca_AX(lambda = c(lambda_list[a], lambda_list[b]), abs = 0.001, Initial_cca, n, group_x, group_y)
      var_num_X[b,a] <- length(which(sapply(unique(group_x), function(j) sum(fit$X_coef[group_x == j] ^ 2)) != 0))
      var_num_Y[b,a] <- length(which(sapply(unique(group_y), function(j) sum(fit$Y_coef[group_y == j] ^ 2)) != 0))
    }
  }
  
  A <- var_num_X >= max_var[1]
  if(length(which(colMeans(A) > 0)) == 0){
    min_X <- max_t / 100000
    fit <- cca_AX(lambda = c(min_X, max_t), abs = 0.001, Initial_cca, n, group_x, group_y)
    while(length(which(sapply(unique(group_x), function(j) sum(fit$X_coef[group_x == j] ^ 2)) != 0)) <= max_var[1]){
      min_X <- min_X / 10
      fit <- cca_AX(lambda = c(min_X, max_t), abs = 0.001, Initial_cca, n, group_x, group_y)
    }
  }else{
    min_X <- lambda_list[min(which(colMeans(A) > 0))]
  }
  
  A <- var_num_Y >= max_var[2]
  if(length(which(rowMeans(A) > 0)) == 0){
    min_Y <- max_t / 100000
    fit <- cca_AX(lambda = c(min_X, min_Y), abs = 0.001, Initial_cca, n, group_x, group_y)
    while(length(which(sapply(unique(group_y), function(j) sum(fit$Y_coef[group_y == j] ^ 2)) != 0)) <= max_var[2]){
      min_Y <- min_Y / 10
      fit <- cca_AX(lambda = c(min_X, min_Y), abs = 0.001, Initial_cca, n, group_x, group_y)
    }
  }else{
    min_Y <- lambda_list[min(which(rowMeans(A) > 0))]
  }
  
  lambda_list_X <- exp(seq(log(max_t), log(min_X), length.out = length_tun + 1))[-1] 
  lambda_list_Y <- matrix(NA, length_tun, length_tun)
  CV_error <- matrix(NA, length(lambda_list_X), length(lambda_list_X))
  
  for(a in 1:length(lambda_list_X)){
    a_X <- c(Initial_cca$Cov %*% Initial_cca$Y_coef) 
    a_X <- c(sapply(1:p_X, function(j){
      norm <- sqrt(sum(a_X[group_x == j] ^ 2))
      if(norm == 0){
        rep(0, sum(group_x == j))
      }else{
        max((1 - lambda_list_X[a] / norm), 0) * a_X[group_x == j]
      }
    })) 
    x_norm <- c(sqrt(crossprod(a_X)))
    if(x_norm != 0){
      a_X <- a_X / x_norm
    }
    
    a_Y <- c(t(a_X) %*% Initial_cca$Cov)
    max_Y <- max(c(sapply(1:p_Y, function(j){
      sqrt(sum(a_Y[group_y == j] ^ 2))
    })))
    if(max_Y == 0){
      max_Y <- max_t
    }
    
    lambda_list_Y[,a] <- exp(seq(log(max_Y), log(min_Y), length.out = length_tun + 1))[-1] 
    for(b in 1:length(lambda_list_Y[,a])){
      if(b > 1){
        if(CV_error[b-1,a] != -10 ^ (11)){
          fit <- cv_slcca(lambda = c(lambda_list_X[a], lambda_list_Y[b,a]), Initial_cca, n, group_x, group_y, max_var)
          CV_error[b,a] <- fit$error
          if(CV_error[b,a] == max(CV_error, na.rm = T)){
            lambda <- c(lambda_list_X[a], lambda_list_Y[b,a])
          }
        }else{
          CV_error[b,a] <- -10 ^ (11)
          if(CV_error[b,a] == max(CV_error, na.rm = T)){
            lambda <- c(lambda_list_X[a], lambda_list_Y[b,a])
          }
        }
        
      }else{
        fit <- cv_slcca(lambda = c(lambda_list_X[a], lambda_list_Y[b,a]), Initial_cca, n, group_x, group_y, max_var)
        CV_error[b,a] <- fit$error
        if(CV_error[b,a] == max(CV_error, na.rm = T)){
          lambda <- c(lambda_list_X[a], lambda_list_Y[b,a])
        }
      }
      
      if(print_tuning == T){
        print(CV_error)
      }
      # gc()
    }
  }
  
  CV_error[CV_error == - 10 ^ 10] <- NA
  CV_error[CV_error == - 10 ^ 11] <- NA
  sd_CV <- 0.1 * sd(CV_error[is.finite(CV_error)])
  
  mark <- (CV_error >= max(CV_error, na.rm = T) - sd_CV)
  
  row_mark <- min(which(apply(mark, c(2), function(x) mean(x, na.rm = T)) != 0))
  lambda_1 <- c(lambda_list_X[row_mark], lambda_list_Y[min(which(mark[,row_mark] == T)), row_mark])
  
  col_mark <- min(which(apply(mark, c(1), function(x) mean(x, na.rm = T)) != 0))
  lambda_2 <- c(lambda_list_X[min(which(mark[col_mark,] == T))], lambda_list_Y[col_mark, min(which(mark[col_mark,] == T))])
  
  fit_opt_1 <- cv_slcca(lambda = lambda_1, Initial_cca, n, group_x, group_y, max_var)
  fit_opt_2 <- cv_slcca(lambda = lambda_2, Initial_cca, n, group_x, group_y, max_var)
  
  if(sum(sapply(unique(group_x), function(j) sum(fit_opt_1$fit_cca$X_coef[group_x == j] ^ 2)) != 0) + 
     sum(sapply(unique(group_y), function(j) sum(fit_opt_1$fit_cca$Y_coef[group_y == j] ^ 2)) != 0) <=
     sum(sapply(unique(group_x), function(j) sum(fit_opt_2$fit_cca$X_coef[group_x == j] ^ 2)) != 0) + 
     sum(sapply(unique(group_y), function(j) sum(fit_opt_2$fit_cca$Y_coef[group_y == j] ^ 2)) != 0)
  ){
    fit_opt <- fit_opt_1
    print(c(lambda_1))
  }else{
    fit_opt <- fit_opt_2
    print(c(lambda_2))
  }
  cor <- fit_opt$error
  
  return(list(fit_opt = fit_opt,
              error = cor,
              CV_error = CV_error,
              var_X = which(sapply(unique(group_x), function(j) sum(fit_opt$fit_cca$X_coef[group_x == j] ^ 2)) != 0),
              var_Y = which(sapply(unique(group_y), function(j) sum(fit_opt$fit_cca$Y_coef[group_y == j] ^ 2)) != 0)
  ))
}

CCVA <- function(grid_num = 201, dat_X, dat_Y, mean_zero = F, print_tuning = T, length_tun = 30, 
                 max_var = c(min(n, p_X), min(n, p_Y)), Sparse = T, comp, FVE = 0.95){
  
  ## First-step estimation
  time_grid <- seq(0, 1, length.out = grid_num)
  n <- dat_X$n
  p_X <- dat_X$p
  p_Y <- dat_Y$p
  
  Optimal_basis <- Optimal_basis_est(dat_X, dat_Y, grid_num, max_comp_num = 5, mean_zero, FVE)
  
  print(c(Optimal_basis$Optimal_bais$X$comp_num, Optimal_basis$Optimal_bais$Y$comp_num))
  group_x <- rep(1:p_X, each = Optimal_basis$Optimal_bais$X$comp_num)
  group_y <- rep(1:p_Y, each = Optimal_basis$Optimal_bais$Y$comp_num)
  
  ## SVD on cross-covariance
  dat_x_g <- t(sapply(1:n, function(i){
    c(sapply(1:p_X, function(j){
      t(Optimal_basis$Optimal_bais$X$phi[Optimal_basis$Time_mark$X[[i]],1:Optimal_basis$Optimal_bais$X$comp_num]) %*% c(dat_X$dat[[j]][[i]]$val - Optimal_basis$Mean$X[[j]]$mu[Optimal_basis$Time_mark$X[[i]]]) / length(dat_X$dat[[j]][[i]]$time)
    }))
  }))
  
  dat_y_g <- t(sapply(1:n, function(i){
    c(sapply(1:p_Y, function(j){
      t(Optimal_basis$Optimal_bais$Y$phi[Optimal_basis$Time_mark$Y[[i]],1:Optimal_basis$Optimal_bais$Y$comp_num]) %*% c(dat_Y$dat[[j]][[i]]$val - Optimal_basis$Mean$Y[[j]]$mu[Optimal_basis$Time_mark$Y[[i]]]) / length(dat_Y$dat[[j]][[i]]$time)
    }))
  }))
  
  Optimal_basis$G_matrix <- sapply(1:n, function(i){
    tcrossprod(dat_x_g[i,], dat_y_g[i,])
  }, simplify = 'array')
  
  Cov <- apply(Optimal_basis$G_matrix, c(1, 2), mean)
  fit <- svd(Cov)
  fit_cca <- list()
  fit_cca$X_coef <- fit$u
  fit_cca$Y_coef <- fit$v
  fit_cca$cov <- fit$d
  tol_cov <- sqrt(sum(fit$d ^ 2))
  rm(fit)
  gc()
  
  if(Sparse == T){
    cv_num <- 5
    cv_mark <- lapply(1:cv_num, function(k){
      floor(seq(k, n, cv_num))
    })
    
    Initial_cca <- list(max_sin = fit_cca$cov[1], X_coef = fit_cca$X_coef[,1], Y_coef = fit_cca$Y_coef[,1],
                        Cov = Cov,
                        CV = lapply(1:cv_num, function(k){
                          list(train = rowMeans(Optimal_basis$G_matrix[,,seq_len(dim(Optimal_basis$G_matrix)[3])[-cv_mark[[k]]]], dims = 2),
                               test = rowMeans(Optimal_basis$G_matrix[,,seq_len(dim(Optimal_basis$G_matrix)[3])[cv_mark[[k]]]], dims = 2))
                        }), cv_num = cv_num
    )
    Optimal_basis$G_matrix <- NULL
    gc()
    
    fit <- SLCCA_tr(group_x = group_x, group_y = group_y, print_tuning, length_tun, max_var, Initial_cca, n, p_X, p_Y)
    gc()
    
    var_X <- list()
    var_Y <- list()
    
    var_X[[1]] <- fit$var_X
    var_Y[[1]] <- fit$var_Y
    
    fit_cca$X_coef[,1] <- fit$fit_opt$fit_cca$X_coef
    
    fit_cca$Y_coef[,1] <- fit$fit_opt$fit_cca$Y_coef
    fit_cca$cov[1] <- t(fit_cca$X_coef[,1]) %*% Cov %*% fit_cca$Y_coef[,1]
    
    if(comp > 1){
      for(r in 2:comp){
        fit_matrix <- fit_cca$cov[r-1] * tcrossprod(fit_cca$X_coef[,r-1], fit_cca$Y_coef[,r-1])
        Cov <- Cov - fit_matrix
        fit <- svd(Cov)
        Initial_cca <- list(max_sin = fit$d[1], X_coef = fit$u[,1], Y_coef = fit$v[,1],
                            Cov = Cov,
                            CV = lapply(1:cv_num, function(k){
                              list(train = Initial_cca$CV[[k]]$train - fit_matrix,
                                   test = Initial_cca$CV[[k]]$test - fit_matrix)
                            }), cv_num = cv_num
        )
        gc()
        
        fit <- SLCCA_tr(group_x = group_x, group_y = group_y, print_tuning, length_tun, max_var, Initial_cca, n, p_X, p_Y)
        gc()
        
        var_X[[r]] <- fit$var_X
        var_Y[[r]] <- fit$var_Y
        
        fit_cca$X_coef[,r] <- fit$fit_opt$fit_cca$X_coef
        
        fit_cca$Y_coef[,r] <- fit$fit_opt$fit_cca$Y_coef
        
        fit_cca$cov[r] <- t(fit_cca$X_coef[,r]) %*% Cov %*% fit_cca$Y_coef[,r]
        
      }
    }
    rm(Initial_cca)
    
  }else{
    var_X <- lapply(1:comp, function(k) 1:p_X)
    var_Y <- lapply(1:comp, function(k) 1:p_Y)
  }
  rm(Cov)
  
  load_X <- lapply(1:comp, function(k){
    matrix(fit_cca$X_coef[,k], p_X, byrow = T) %*% t(Optimal_basis$Optimal_bais$X$phi[,1:Optimal_basis$Optimal_bais$X$comp_num])
  })
  load_Y <- lapply(1:comp, function(k){
    matrix(fit_cca$Y_coef[,k], p_Y, byrow = T) %*% t(Optimal_basis$Optimal_bais$Y$phi[,1:Optimal_basis$Optimal_bais$Y$comp_num])
  })
  
  Z_X <- lapply(1:comp, function(k){
    sapply(1:n, function(i){
      sum(sapply(1:p_X, function(j){
        mean((dat_X$dat[[j]][[i]]$val - Optimal_basis$Mean$X[[j]]$mu[Optimal_basis$Time_mark$X[[i]]]) * load_X[[k]][j,Optimal_basis$Time_mark$X[[i]]])
      }))
    })
  })
  
  Z_Y <- lapply(1:comp, function(k){
    sapply(1:n, function(i){
      sum(sapply(1:p_Y, function(j){
        mean((dat_Y$dat[[j]][[i]]$val - Optimal_basis$Mean$Y[[j]]$mu[Optimal_basis$Time_mark$Y[[i]]]) * load_Y[[k]][j,Optimal_basis$Time_mark$Y[[i]]])
      }))
    })
  })
  
  return(list(Mean = Optimal_basis$Mean,
              Time_mark = Optimal_basis$Time_mark,
              Optimal_basis = Optimal_basis$Optimal_bais,
              var_X = var_X,
              var_Y = var_Y,
              load_X = load_X,
              load_Y = load_Y,
              Z_X = Z_X,
              Z_Y = Z_Y,
              eta = fit_cca$cov[1:comp],
              tol_cov = tol_cov
  ))
}

plot_CCVA <- function(fit, dat_X, dat_Y, hight_feature = 5, k, name_X, name_Y, time_illu, time_label, title_name){
  
  ## Longitudinal canonical vector
  Canonical_X <- fit$load_X[[k]]
  Canonical_X <- Canonical_X * sign(sum(Canonical_X))
  Canonical_Y <- fit$load_Y[[k]]
  Canonical_Y <- Canonical_Y * sign(sum(Canonical_Y))
  
  ### X
  order_X <- order(rowSums(Canonical_X ^ 2), decreasing = T)
  nonzeo_X <- sum(sort(rowSums(Canonical_X ^ 2), decreasing = T) > 0)
  feature_X <- dat_X$feature_label[order_X[1:nonzeo_X]]
  
  mark <- sapply(1:length(time_illu), function(h) which.min(abs(fit$Optimal_basis$X$workGrid - time_illu[h]))) 
  Lcv_X <- matrix(sapply(1:min(hight_feature, nonzeo_X), function(l){
    Canonical_X[order_X[l],mark]
  }), ncol = min(hight_feature, nonzeo_X))
  
  imp_fea <- sapply(1:dim(Lcv_X)[2], function(k) dat_X$feature_label[order_X[k]])
  imp_fea <- factor(c(matrix(imp_fea, nrow = nrow(Lcv_X), ncol = length(imp_fea), byrow = T)),
                    levels = rev(imp_fea))
  dat_plot_X <- data.frame(time = rep(time_illu, ncol(Lcv_X)),
                           val  = c(Lcv_X), 
                           label = imp_fea,
                           data_type = rep(name_X, length(imp_fea))
  )
  
  ### Y
  order_Y <- order(rowSums(Canonical_Y ^ 2), decreasing = T)
  nonzeo_Y <- sum(sort(rowSums(Canonical_Y ^ 2), decreasing = T) > 0)
  feature_Y <- dat_Y$feature_label[order_Y[1:nonzeo_Y]]
  
  Lcv_Y <- matrix(sapply(1:min(hight_feature, nonzeo_Y), function(l){
    Canonical_Y[order_Y[l],mark]
  }), ncol = min(hight_feature, nonzeo_Y))
  
  imp_fea <- sapply(1:dim(Lcv_Y)[2], function(k) dat_Y$feature_label[order_Y[k]])
  imp_fea <- factor(c(matrix(imp_fea, nrow = dim(Lcv_Y)[1], ncol = length(imp_fea), byrow = T)),
                    levels = rev(imp_fea))
  dat_plot_Y <- data.frame(time = rep(time_illu, ncol(Lcv_Y)),
                           val  = c(Lcv_Y), 
                           label = imp_fea,
                           data_type = rep(name_Y, length(imp_fea))
  )
  
  dat_plot <- rbind(dat_plot_X, dat_plot_Y)
  
  pl <- ggplot(dat_plot) + 
    geom_tile(aes(x = time, y = label, fill = val)) +
    scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "orange", midpoint = 0) +
    # geom_line(aes(x = time, y = val, group = label), color = "gray", alpha = 0.3) +
    # geom_line(aes(x = time, y =  high_fea, group = label), color = "black") +
    # ggrepel::geom_text_repel(aes(x = high_x, y = high_y, label = label)) +
    facet_wrap(~data_type, nrow = 1, scales = "free_y") +
    # geom_hline(yintercept = 0, color = "blue", linetype = 2) +
    theme_bw(base_family = "Times") +
    scale_x_continuous(breaks = time_illu, labels = time_label) +
    xlab("") + ylab("") + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = "right",
          text = element_text(size = 12),
          panel.border = element_blank(),
          # text = element_text(family = "STHeiti"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 0)) +
    labs(title = title_name, fill = "")
  
  return(list(plot = pl,
              feature_X = feature_X,
              feature_Y = feature_Y))
}

group_tran <- function(dat, p, n){
  
  Time <- lapply(1:p, function(j){
    unlist(sapply(1:n, function(i){
      dat$dat[[j]][[i]]$time
    }))
  })
  
  mark <- lapply(1:p, function(j){
    which(sapply(dat$time_grid, function(t){
      sum(Time[[j]] == t)
    }) == n)
  })
  
  label <- 1:n
  
  dat <- t(sapply(1:length(label), function(i){
    unlist(lapply(1:p, function(j){
      mean(dat$dat[[j]][[i]]$val)
    }))
  }))
  
  group <- 1:p
  
  return(list(dat = dat, group = group, label = 1:n))
}

simu_func <- function(n, rho, time_num, p_X, p_Y, comp, p_ture, seed){
  
  error <- T
  
  while(error == T){
    tryCatch({
      seed <- seed * 0.8 + 1000
      set.seed(seed)
      dat <- data_gen(n, rho, time_num, p_X, p_Y, comp, p_ture)
      
      dat_X <- list(dat = dat$dat_X, n = dat$n, p = dat$p_X, feature_label = 1:dat$p_X)
      dat_Y <- list(dat = dat$dat_Y, n = dat$n, p = dat$p_Y, feature_label = 1:dat$p_Y)
      
      dat_X$time_grid <- seq(0, 1, length.out = time_num)
      dat_Y$time_grid <- seq(0, 1, length.out = time_num)
      
      # LCCA
      dat_x <- NULL
      for(i in 1:n){
        time_i <- unique(unlist(lapply(1:p_X, function(j){
          dat_X$dat[[j]][[i]]$time
        })))
        
        mark <- sapply(1:length(time_i), function(k){
          prod(sapply(1:p_X, function(j){
            sum(dat_X$dat[[j]][[i]]$time == time_i[k])
          }))
        })
        time_i <- time_i[mark == 1]
        
        dat_x <- rbind(dat_x, cbind(rep(i, length(time_i)), 
                                    time_i,
                                    sapply(1:p_X, function(j){
                                      sapply(1:length(time_i), function(k){
                                        time_ob <- time_i[k]
                                        mark <- which(time_ob == dat_X$dat[[j]][[i]]$time)
                                        dat_X$dat[[j]][[i]]$val[mark]
                                      })
                                    })))
      }
      
      dat_y <- NULL
      for(i in 1:n){
        time_i <- unique(unlist(lapply(1:p_Y, function(j){
          dat_Y$dat[[j]][[i]]$time
        })))
        mark <- sapply(1:length(time_i), function(k){
          prod(sapply(1:p_Y, function(j){
            sum(dat_Y$dat[[j]][[i]]$time == time_i[k])
          }))
        })
        time_i <- time_i[mark == 1]
        
        dat_y <- rbind(dat_y, cbind(rep(i, length(time_i)), 
                                    time_i,
                                    sapply(1:p_Y, function(j){
                                      sapply(1:length(time_i), function(k){
                                        time_ob <- time_i[k]
                                        mark <- which(time_ob == dat_Y$dat[[j]][[i]]$time)
                                        dat_Y$dat[[j]][[i]]$val[mark]
                                      })
                                    })
        ))
      }
      
      x <- list(X = as.matrix(t(dat_x[,-c(1:2)])), 
                time = dat_x[,2], 
                I = n, 
                J = nrow(dat_x), 
                visit = sapply(1:n, function(i) length(which(dat_x[,1] == i))))
      
      y <- list(X = as.matrix(t(dat_y[,-c(1:2)])),
                time = dat_y[,2], 
                I = n,
                J = nrow(dat_y), 
                visit = sapply(1:n, function(i) length(which(dat_y[,1] == i))))
      rm(dat_x, dat_y)
      
      Result_LCCA <- lcca.linear_c(x = x, y = y)
      rm(x, y)
      gc()
      
      # CCVA
      Result_SLCCA <- CCVA(grid_num = 201, dat_X, dat_Y, mean_zero = F, print_tuning = F, length_tun = 30, 
                           max_var = c(2 * p_ture, 2 * p_ture), Sparse = T, comp = 1)
      gc()
      
      # SCCA
      tran_X <- group_tran(dat_X, p_X, n)
      dat_x <- tran_X$dat
      group_x <- tran_X$group
      
      tran_Y <- group_tran(dat_Y, p_Y, n)
      dat_y <- tran_Y$dat
      group_y <- tran_Y$group
      
      label <- intersect(tran_X$label, tran_Y$label)
      
      Result_GSCCA <- GSCCA(dat_x, dat_y, group_x, group_y, print_tuning = F, length_tun = 30, max_var = c(2 * p_ture, 2 * p_ture))
      gc()
      
      # FSCCA-group
      norm_ploy <- sqrt(colMeans(poly(seq(0, 1, length.out = 201), degree = 3) ^ 2))
      dat_x_g <- t(sapply(1:n, function(i){
        c(sapply(1:p_X, function(j){
          (t(poly(dat_X$dat[[j]][[i]]$time, degree = 3)) / norm_ploy) %*% dat_X$dat[[j]][[i]]$val / length(dat_X$dat[[j]][[i]]$time)
        }))
      }))
      group_x_g <- rep(1:p_X, each = 3)
      
      dat_y_g <- t(sapply(1:n, function(i){
        c(sapply(1:p_Y, function(j){
          (t(poly(dat_Y$dat[[j]][[i]]$time, degree = 3)) / norm_ploy) %*% dat_Y$dat[[j]][[i]]$val / length(dat_Y$dat[[j]][[i]]$time)
        }))
      }))
      group_y_g <- rep(1:p_Y, each = 3)
      
      Result_FSCCA <- GSCCA(dat_x_g, dat_y_g, group_x_g, group_y_g, print_tuning = F, length_tun = 30, max_var = c(2 * p_ture, 2 * p_ture))
      gc()
      
      # MFPCA
      X <- multiFunData(lapply(1:p_X, function(j){
        funData(argvals = seq(0, 1, length.out = 201),
                X = t(sapply(1:n, function(i){
                  obs <- rep(NA, 201)
                  obs[Result_SLCCA$Time_mark$X[[i]]] <- dat_X$dat[[j]][[i]]$val
                  return(obs)
                })))
      }))
      MFPCA_X <- MFPCA(X, M = 1, uniExpansions = lapply(1:p_X, function(j) list(type = "uFPCA", pve = 0.9, nbasis = 8)))
      
      X <- multiFunData(lapply(1:p_Y, function(j){
        funData(argvals = seq(0, 1, length.out = 201),
                X = t(sapply(1:n, function(i){
                  obs <- rep(NA, 201)
                  obs[Result_SLCCA$Time_mark$Y[[i]]] <- dat_Y$dat[[j]][[i]]$val
                  return(obs)
                })))
      }))
      MFPCA_Y <- MFPCA(X, M = 1, uniExpansions = lapply(1:p_Y, function(j) list(type = "uFPCA", pve = 0.9, nbasis = 8)))
      rm(X)
      gc()
      
      # Correlation with true scores
      cor_SLCCA_X <- abs(cor(Result_SLCCA$Z_X[[1]], dat$val[[1]][1,]))
      cor_SLCCA_Y <- abs(cor(Result_SLCCA$Z_Y[[1]], dat$val[[1]][2,]))
      
      cor_GSCCA_X <- c(abs(cor(Result_GSCCA$X_cca, dat$val[[1]][1,label])))
      cor_GSCCA_Y <- c(abs(cor(Result_GSCCA$Y_cca, dat$val[[1]][2,label])))
      
      cor_LCCA_X <- abs(cor(matrix(Result_LCCA$scores$x, n)[,1], dat$val[[1]][1,]))
      cor_LCCA_Y <- abs(cor(matrix(Result_LCCA$scores$y, n)[,1], dat$val[[1]][2,]))
      
      cor_FSCCA_X <- c(abs(cor(Result_FSCCA$X_cca, dat$val[[1]][1,label])))
      cor_FSCCA_Y <- c(abs(cor(Result_FSCCA$Y_cca, dat$val[[1]][2,label])))
      
      cor_MFPCA_X <- c(abs(cor(MFPCA_X$scores, dat$val[[1]][1,])))
      cor_MFPCA_Y <- c(abs(cor(MFPCA_Y$scores, dat$val[[1]][2,])))
      
      # Variable selection
      FP_X_SLCCA <- sum(sapply(1:p_ture, function(j){sum(j == Result_SLCCA$var_X[[1]])}) == 0)
      if(length((1:p_X)[-Result_SLCCA$var_X[[1]]]) != 0){
        FN_X_SLCCA <- sum(sapply((1:p_X)[-Result_SLCCA$var_X[[1]]], function(j){sum(j == (1:p_X)[-(1:p_ture)])}))
      }else{
        FN_X_SLCCA <- 0
      }
      
      FP_Y_SLCCA <- sum(sapply(1:p_ture, function(j){sum(j == Result_SLCCA$var_Y[[1]])}) == 0)
      if(length((1:p_Y)[-Result_SLCCA$var_Y[[1]]]) != 0){
        FN_Y_SLCCA <- sum(sapply((1:p_Y)[-Result_SLCCA$var_Y[[1]]], function(j){sum(j == (1:p_Y)[-(1:p_ture)])}))
      }else{
        FN_Y_SLCCA <- 0
      }
      
      FP_X_GSCCA <- sum(sapply(1:p_ture, function(j){sum(j == Result_GSCCA$var_X)}) == 0)
      if(length((1:p_X)[-Result_GSCCA$var_X]) != 0){
        FN_X_GSCCA <- sum(sapply((1:p_X)[-Result_GSCCA$var_X], function(j){sum(j == (1:p_X)[-(1:p_ture)])}))
      }else{
        FN_X_GSCCA <- 0
      }
      
      FP_Y_GSCCA <- sum(sapply(1:p_ture, function(j){sum(j == Result_GSCCA$var_Y)}) == 0)
      if(length((1:p_Y)[-Result_GSCCA$var_Y]) != 0){
        FN_Y_GSCCA <- sum(sapply((1:p_Y)[-Result_GSCCA$var_Y], function(j){sum(j == (1:p_Y)[-(1:p_ture)])}))
      }else{
        FN_Y_GSCCA <- 0
      }
      
      FP_X_FSCCA <- sum(sapply(1:p_ture, function(j){sum(j == Result_FSCCA$var_X)}) == 0)
      if(length((1:p_X)[-Result_FSCCA$var_X]) != 0){
        FN_X_FSCCA <- sum(sapply((1:p_X)[-Result_FSCCA$var_X], function(j){sum(j == (1:p_X)[-(1:p_ture)])}))
      }else{
        FN_X_FSCCA <- 0
      }
      
      FP_Y_FSCCA <- sum(sapply(1:p_ture, function(j){sum(j == Result_FSCCA$var_Y)}) == 0)
      if(length((1:p_Y)[-Result_FSCCA$var_Y]) != 0){
        FN_Y_FSCCA <- sum(sapply((1:p_Y)[-Result_FSCCA$var_Y], function(j){sum(j == (1:p_Y)[-(1:p_ture)])}))
      }else{
        FN_Y_FSCCA <- 0
      }
      
      # CCA loading
      # CCVA
      Load_err_X_SLCCA <- sum(sapply(1:p_X, function(j){
        if(sum(abs(dat$load_X[[1]][j,])) != 0){
          sum((sign(sum(Result_SLCCA$load_X[[1]][j,] * dat$load_X[[1]][j,])) * Result_SLCCA$load_X[[1]][j,] - dat$load_X[[1]][j,]) ^ 2) * 0.005
        }else{
          sum((Result_SLCCA$load_X[[1]][j,] - dat$load_X[[1]][j,]) ^ 2) * 0.005
        }
      })) /  sum(sapply(1:p_X, function(j){
        sum((dat$load_X[[1]][j,]) ^ 2) * 0.005
      }))
      
      Load_err_Y_SLCCA <- sum(sapply(1:p_Y, function(j){
        if(sum(abs(dat$load_Y[[1]][j,])) != 0){
          sum((sign(sum(Result_SLCCA$load_X[[1]][j,] * dat$load_X[[1]][j,])) * Result_SLCCA$load_X[[1]][j,] - dat$load_X[[1]][j,]) ^ 2) * 0.005
        }else{
          sum((Result_SLCCA$load_Y[[1]][j,] - dat$load_Y[[1]][j,]) ^ 2) * 0.005
        }
      })) /  sum(sapply(1:p_Y, function(j){
        sum((dat$load_Y[[1]][j,]) ^ 2) * 0.005
      }))
      
      # LCCA
      load_X <- matrix(Result_LCCA$re.X$phix0, nrow = p_X)[,1] +  matrix(Result_LCCA$re.X$phix1, nrow = p_X)[,1] %*% t(seq(0, 1, 1 / 200))
      norm_X <- sum(sapply(1:p_X, function(j){
        sum((load_X[j,]) ^ 2) * 0.005
      }))
      load_X <- load_X / norm_X
      Load_err_X_LCCA <- sum(sapply(1:p_X, function(j){
        if(sum(abs(dat$load_X[[1]][j,])) != 0){
          sum((sign(sum(load_X[j,] * dat$load_X[[1]][j,])) * load_X[j,] - dat$load_X[[1]][j,]) ^ 2) * 0.005
        }else{
          sum((load_X[j,] - dat$load_X[[1]][j,]) ^ 2) * 0.005
        }
        
      })) /  sum(sapply(1:p_X, function(j){
        sum((dat$load_X[[1]][j,]) ^ 2) * 0.005
      }))
      
      load_Y <- matrix(Result_LCCA$re.Y$phix0, nrow = p_Y)[,1] +  matrix(Result_LCCA$re.Y$phix1, nrow = p_Y)[,1] %*% t(seq(0, 1, 1 / 200))
      norm_Y <- sum(sapply(1:p_Y, function(j){
        sum((load_Y[j,]) ^ 2) * 0.005
      }))
      load_Y <- load_Y / norm_Y
      Load_err_Y_LCCA <- sum(sapply(1:p_Y, function(j){
        if(sum(abs(dat$load_Y[[1]][j,])) != 0){
          sum((sign(sum(load_Y[j,] * dat$load_Y[[1]][j,])) * load_Y[j,] - dat$load_Y[[1]][j,]) ^ 2) * 0.005
        }else{
          sum((load_Y[j,] - dat$load_Y[[1]][j,]) ^ 2) * 0.005
        }
      })) /  sum(sapply(1:p_Y, function(j){
        sum((dat$load_Y[[1]][j,]) ^ 2) * 0.005
      }))
      
      # FSCCA
      load_X <- matrix(Result_FSCCA$fit_opt$fit_cca$X_coef, ncol = 3, byrow = T) %*% (t(poly(seq(0, 1, 1 / 200), degree = 3)) / norm_ploy)
      norm_X <- sum(sapply(1:p_X, function(j){
        sum((load_X[j,]) ^ 2) * 0.005
      }))
      load_X <- load_X / norm_X
      Load_err_X_FSCCA <- sum(sapply(1:p_X, function(j){
        if(sum(abs(dat$load_X[[1]][j,])) != 0){
          sum((sign(sum(load_X[j,] * dat$load_X[[1]][j,])) * load_X[j,] - dat$load_X[[1]][j,]) ^ 2) * 0.005
        }else{
          sum((load_X[j,] - dat$load_X[[1]][j,]) ^ 2) * 0.005
        }
        
      })) /  sum(sapply(1:p_X, function(j){
        sum((dat$load_X[[1]][j,]) ^ 2) * 0.005
      }))
      
      load_Y <- matrix(Result_FSCCA$fit_opt$fit_cca$Y_coef, ncol = 3, byrow = T) %*% (t(poly(seq(0, 1, 1 / 200), degree = 3)) / norm_ploy)
      norm_Y <- sum(sapply(1:p_Y, function(j){
        sum((load_Y[j,]) ^ 2) * 0.005
      }))
      load_Y <- load_Y / norm_Y
      Load_err_Y_FSCCA <- sum(sapply(1:p_Y, function(j){
        if(sum(abs(dat$load_Y[[1]][j,])) != 0){
          sum((sign(sum(load_Y[j,] * dat$load_Y[[1]][j,])) * load_Y[j,] - dat$load_Y[[1]][j,]) ^ 2) * 0.005
        }else{
          sum((load_Y[j,] - dat$load_Y[[1]][j,]) ^ 2) * 0.005
        }
      })) /  sum(sapply(1:p_Y, function(j){
        sum((dat$load_Y[[1]][j,]) ^ 2) * 0.005
      }))
      
      # MFPCA
      load_X <- t(sapply(1:p_X, function(j){
        MFPCA_X$functions[[j]]@X
      }))
      norm_X <- sum(sapply(1:p_X, function(j){
        sum((load_X[j,]) ^ 2) * 0.005
      }))
      load_X <- load_X / norm_X
      Load_err_X_MFPCA <- sum(sapply(1:p_X, function(j){
        if(sum(abs(dat$load_X[[1]][j,])) != 0){
          sum((sign(sum(load_X[j,] * dat$load_X[[1]][j,])) * load_X[j,] - dat$load_X[[1]][j,]) ^ 2) * 0.005
        }else{
          sum((load_X[j,] - dat$load_X[[1]][j,]) ^ 2) * 0.005
        }
      })) /  sum(sapply(1:p_X, function(j){
        sum((dat$load_X[[1]][j,]) ^ 2) * 0.005
      }))
      
      load_Y <- t(sapply(1:p_Y, function(j){
        MFPCA_Y$functions[[j]]@X
      }))
      norm_Y <- sum(sapply(1:p_Y, function(j){
        sum((load_Y[j,]) ^ 2) * 0.005
      }))
      load_Y <- load_Y / norm_Y
      Load_err_Y_MFPCA <- sum(sapply(1:p_Y, function(j){
        if(sum(abs(dat$load_Y[[1]][j,])) != 0){
          sum((sign(sum(load_Y[j,] * dat$load_Y[[1]][j,])) * load_Y[j,] - dat$load_Y[[1]][j,]) ^ 2) * 0.005
        }else{
          sum((load_Y[j,] - dat$load_Y[[1]][j,]) ^ 2) * 0.005
        }
      })) /  sum(sapply(1:p_Y, function(j){
        sum((dat$load_Y[[1]][j,]) ^ 2) * 0.005
      }))
      
      return(list(
        Result_SLCCA = Result_SLCCA,
        Result_GSCCA = Result_GSCCA,
        Result_LCCA = Result_LCCA,
        Result_FSCCA = Result_FSCCA,
        cor_error = list(
          X = list(cor_SLCCA = cor_SLCCA_X,
                   cor_GSCCA = cor_GSCCA_X,
                   cor_LCCA = cor_LCCA_X,
                   cor_MFPCA = cor_MFPCA_X,
                   cor_FSCCA = cor_FSCCA_X
          ),
          Y = list(cor_SLCCA = cor_SLCCA_Y,
                   cor_GSCCA = cor_GSCCA_Y,
                   cor_LCCA = cor_LCCA_Y,
                   cor_MFPCA = cor_MFPCA_Y,
                   cor_FSCCA = cor_FSCCA_Y
          )
        ),
        var_error = list(
          X = list(
            SLCCA = list(FP = FP_X_SLCCA, FN = FN_X_SLCCA),
            GSCCA = list(FP = FP_X_GSCCA, FN = FN_X_GSCCA),
            FSCCA = list(FP = FP_X_FSCCA, FN = FN_X_FSCCA)
          ),
          Y = list(
            SLCCA = list(FP = FP_Y_SLCCA, FN = FN_Y_SLCCA),
            GSCCA = list(FP = FP_Y_GSCCA, FN = FN_Y_GSCCA),
            FSCCA = list(FP = FP_Y_FSCCA, FN = FN_Y_FSCCA)
          )
        ),
        loa_error = list(
          X = list(
            Load_err_X_SLCCA = Load_err_X_SLCCA,
            Load_err_X_LCCA = Load_err_X_LCCA,
            Load_err_X_FSCCA = Load_err_X_FSCCA,
            Load_err_X_MFPCA = Load_err_X_MFPCA
          ),
          Y = list(
            Load_err_Y_SLCCA = Load_err_Y_SLCCA,
            Load_err_Y_LCCA = Load_err_Y_LCCA,
            Load_err_Y_FSCCA = Load_err_Y_FSCCA,
            Load_err_Y_MFPCA = Load_err_Y_MFPCA
          )
        )
      ))
      
      error <- F
    }, error = function(e) e)
  }
  
}