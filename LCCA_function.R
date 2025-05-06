lcca.linear_c <- function (x, y, varthresh = 0.95, projectthresh = 1, method = "Wilks", 
                           verbose = FALSE) 
{
  re.X <- hd_lfpca(x$X, T = x$time, J = sum(x$visit), I = x$I, 
                   visit = x$visit, varthresh = varthresh, projectthresh = projectthresh, 
                   timeadjust = FALSE, verbose = verbose)
  re.Y <- hd_lfpca(Y=y$X, T = y$time, J = y$J, I = y$I, visit = y$visit, 
                   varthresh = varthresh, projectthresh = projectthresh, 
                   timeadjust = FALSE, verbose = verbose)
  cca = cc(t(re.X$xi), t(re.Y$xi))
  # tests = p.asym(rho = cca$cor, N = x$I, p = re.X$Nx, q = re.Y$Nx, 
  #                tstat = method)
  # tests$r = cca$cor
  ccor.dim = length(cca$cor)
  ccor = cca$cor[1:ccor.dim]
  xcv_x0 = re.X$phix0 %*% cca$xcoef[, 1:ccor.dim]
  xcv_x1 = re.X$phix1 %*% cca$xcoef[, 1:ccor.dim]
  xcv_y0 = re.Y$phix0 %*% cca$ycoef[, 1:ccor.dim]
  xcv_y1 = re.Y$phix1 %*% cca$ycoef[, 1:ccor.dim]
  scores = list(x = cca$scores$xscores[, 1:ccor.dim], y = cca$scores$yscores[, 
                                                                             1:ccor.dim])
  out = list(ccor.dim = ccor.dim, ccor = ccor, 
             xcv_x0 = xcv_x0, xcv_x1 = xcv_x1, xcv_y0 = xcv_y0, xcv_y1 = xcv_y1, 
             scores = scores, re.X = re.X, re.Y = re.Y)
  return(out)
}

hd_lfpca <- function (Y, T, J, I, visit, verbose = FALSE, prefix = date(), 
                      Nx = NA, Nw = NA, varthresh = 0.95, projectthresh = 1, timeadjust = FALSE) 
{
  Y = round(t(scale(t(Y), center = TRUE, scale = FALSE)), digit = 5)
  if (dim(Y)[1] >= dim(Y)[2]) {
    svdy = svd(t(Y) %*% Y)
    if (verbose == TRUE) {
      print("Reduce Dimension")
    }
    N = J_projected = sum(cumsum(svdy$d^2)/sum(svdy$d^2) < 
                            projectthresh)
    if (verbose == TRUE) {
      print(J_projected)
    }
    U = svdy$v[, 1:J_projected]
    D = diag(sqrt(svdy$d[1:J_projected]))
    S = diag(svdy$d[1:J_projected])
    V = Y %*% U %*% solve(D)
    Ynew = t(V[, 1:J_projected]) %*% Y
    rm(list = c("svdy"))
  }
  else {
    ymatrix <- t(Y) %*% Y
    svdy <- NULL
    while(length(svdy) == 0){
      tryCatch({
        svdy = svd(ymatrix)
      }, error= function(e) e)
      ymatrix <- ymatrix + matrix(rnorm(length(ymatrix), 0, sqrt(sum(ymatrix ^ 2) * 0.0001)), nrow = nrow(ymatrix))
    }
    rm(ymatrix)
    # svdy <- svd(t(Y) %*% Y)
    
    if (verbose == TRUE) {
      print("Do not Reduce Dimension")
    }
    N = J_projected = nrow(Y) - 1
    if (verbose == TRUE) {
      print(J_projected)
    }
    U = svdy$v[, 1:J_projected]
    D = diag(sqrt(svdy$d[1:J_projected]))
    S = diag(svdy$d[1:J_projected])
    V = Y %*% U %*% solve(D)
    Ynew = t(V[, 1:J_projected]) %*% Y
    rm(list = c("svdy"))
  }
  Yvec = matrix(0, J_projected^2, sum(visit^2))
  X = matrix(0, sum(visit^2), 5)
  k = 0
  J_sq = 0
  for (j in 1:I) {
    s = 1
    Ji = visit[j]
    Ti = T[(1:Ji) + k]
    if (timeadjust == TRUE) {
      Ti = (Ti - mean(Ti))/sqrt(var(Ti))
    }
    Yvec_i = matrix(0, J_projected * J_projected, Ji^2)
    X_i = matrix(0, Ji * Ji, 5)
    Ytmp = Ynew[, (1:Ji) + k]
    for (j1 in 1:Ji) {
      for (j2 in 1:Ji) {
        Yvec_i[, s] = as.vector(Ytmp[, j1] %*% t(Ytmp[, 
                                                      j2]))
        X_i[s, ] = c(1, Ti[j2], Ti[j1], Ti[j1] * Ti[j2], 
                     as.numeric(j1 == j2))
        s = s + 1
      }
    }
    Yvec[, (1:(visit[j]^2)) + J_sq] = Yvec_i
    X[(1:visit[j]^2) + J_sq, ] = X_i
    J_sq = J_sq + Ji^2
    k = k + Ji
  }
  if (verbose == TRUE) {
    print("Estimate Covariance Functions")
  }
  beta <- Yvec %*% X %*% solve(t(X) %*% X)
  K00 = matrix(beta[, 1], J_projected, J_projected)
  K01 = matrix(beta[, 2], J_projected, J_projected)
  K10 = matrix(beta[, 3], J_projected, J_projected)
  K11 = matrix(beta[, 4], J_projected, J_projected)
  Kw = matrix(beta[, 5], J_projected, J_projected)
  Kx = rbind(cbind(K00, K01), cbind(K10, K11))
  Ax = eigen(Kx)
  Aw = eigen(Kw)
  if (is.na(Nx) | is.na(Nw)) {
    lim = 1
    Nw0 = Nx0 = Nw = Nx = max(min(J_projected, sum(Ax$values > 
                                                     0), sum(Aw$values > 0)), 1)
    iter <- 1
    while ((lim > varthresh) & (iter < 1000)) {
      iter <- 1 + iter
      if (min(Aw$values[1:Nw]) < min(Ax$values[1:Nx])) {
        Nw = max(Nw0 - 1, 1)
      }
      else {
        Nx = max(Nx0 - 1, 1)
      }
      lim = (sum(Ax$values[1:Nx]) + sum(Aw$values[1:Nw]))/(sum(Ax$values[1:sum(Ax$values > 
                                                                                 0)]) + sum(Aw$values[1:sum(Aw$values > 0)]))
      Nx0 <- Nx
      Nw0 <- Nw
      if (verbose == TRUE) {
        print(c(Nx, Nw, lim))
      }
    }
  }
  Ax0 = as.matrix(Ax$vectors[1:J_projected, 1:Nx])
  Ax1 = as.matrix(Ax$vectors[1:J_projected + J_projected, 1:Nx])
  S_x = Ax$values[1:Nx]
  S_u = Aw$values[1:Nw]
  Au = Aw$vectors[, 1:Nw]
  phix0 = V %*% Ax0
  phix1 = V %*% Ax1
  phiw = V %*% Aw$vectors[, 1:Nw]
  b = rainbow(200)
  phix = rbind(phix0, phix1)
  total_lambda = sum(Ax$values[1:Nx]) + sum(Aw$values[1:Nw])
  a = (diag(t(phix0) %*% phix0))
  xi_est = matrix(0, Nx, I)
  zeta_est = matrix(0, Nw, J)
  C00 = t(Ax0) %*% Ax0
  C11 = t(Ax1) %*% Ax1
  C01 = t(Ax0) %*% Ax1
  C10 = t(C01)
  C0W = t(Ax0) %*% Au
  C1W = t(Ax1) %*% Au
  if (verbose == TRUE) {
    print("Calcuate Scores")
  }
  Ut = t(U)
  k = 0
  for (i in 1:I) {
    Ji = visit[i]
    Ti = T[k + 1:Ji]
    if (timeadjust == TRUE) {
      Ti = (Ti - mean(Ti))/sqrt(var(Ti))
    }
    Ti_dot = sum(Ti)
    Ti_sq_dot = sum(Ti^2)
    Ui = Ut[, k + 1:Ji]
    ones_Ji = rep(1, Ji)
    BBi_11 = Ji * C00 + Ti_sq_dot * C11 + Ti_dot * (C10 + 
                                                      C01)
    BBi_12 = kronecker(t(ones_Ji), C0W) + kronecker(t(Ti), 
                                                    C1W)
    BBi_21 = t(BBi_12)
    BBi_22 = diag(rep(1, Nw * Ji))
    C1 = ginv(BBi_11 - BBi_12 %*% BBi_22 %*% BBi_21)
    C2 = BBi_22 + BBi_21 %*% ginv(BBi_11 - BBi_12 %*% BBi_21) %*% 
      BBi_12
    C3 = ginv(BBi_11) %*% BBi_12 %*% C2
    BYi_11 = as.vector(t(Ax0) %*% sqrt(S[1:N, ]) %*% Ui %*% 
                         ones_Ji + t(Ax1) %*% sqrt(S[1:N, ]) %*% Ui %*% Ti)
    BYi_12 = as.vector(t(Au) %*% sqrt(S[1:N, ]) %*% Ui)
    BYi = c(BYi_11, BYi_12)
    wi <- rbind(cbind(C1, -C3), cbind(-t(C3), C2)) %*% BYi
    xi_est[, i] = wi[1:Nx]
    zeta_i_vec = wi[Nx + 1:(Ji * Nw)]
    for (j in 1:Ji) {
      zeta_est[, k + j] = zeta_i_vec[((j - 1) * Nw + 1):(j * 
                                                           Nw)]
    }
    k = k + Ji
  }
  if (verbose == TRUE) {
    print("Compute Residual to the Demeaned Data.")
  }
  tmpx = matrix(0, nrow(Ax0), J)
  cumvisit = 0
  for (j in 1:I) {
    for (k in 1:visit[j]) {
      tmpx[, cumvisit + k] = (Ax0 + Ax1 * T[cumvisit + 
                                              k]) %*% xi_est[, j]
    }
    cumvisit = sum(visit[1:j])
  }
  residual = sum((Ynew - tmpx - Aw$vectors[, 1:Nw] %*% zeta_est)^2)
  if (verbose == TRUE) {
    print(paste("Residual of LFPCA Model is:", residual))
  }
  result = list(xi = xi_est, zeta = zeta_est, phix0 = phix0, 
                phix1 = phix1, Nx = Nx, Nw = Nw, phiw = phiw, sx = Ax$values[1:Nx], 
                sw = Aw$values[1:Nw], residual = residual, tij = T)
  return(result)
}

