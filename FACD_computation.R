################################################################################
# FACD convergence and initialization stability analysis
################################################################################

# Set address
rm(list = ls())

mydir <- dirname(
  rstudioapi::getSourceEditorContext()$path
)

setwd(mydir)


################################################################################
# Settings
################################################################################

library(parallel)
library(ggplot2)

source("Function.R")

n_rep <- 100
n_random_init <- 20
n_core <- 10

tol <- 0.001
max_iter <- 200

# Each simulation task is attempted at most twice:
# one original seed and one alternative seed.
max_task_attempts <- 2L

# Seed increment used for the second attempt.
# This value keeps the retry seed different from all original seeds.
retry_seed_offset <- 700000000

# Fixed L2 magnitude of the perturbation added to the SVD initialization
perturbation_level <- 0.1

# Fixed group-thresholding parameters
lambda <- c(0.01, 0.01)

# Common simulation setting
base_setting <- list(
  n = 150,
  rho = 0.8,
  time_num = 8,
  comp = 20,
  p_true = 10,
  grid_num = 201,
  FVE = 0.95
)

# Three dimension settings
dimension_ratio <- c(
  0.5,
  1,
  1.5
)

case_name <- c(
  "0.5n",
  "1.0n",
  "1.5n"
)

case_setting <- lapply(
  seq_along(dimension_ratio),
  function(k) {
    
    list(
      case_id = k,
      case = case_name[k],
      dimension_ratio = dimension_ratio[k],
      n = base_setting$n,
      rho = base_setting$rho,
      time_num = base_setting$time_num,
      
      p_X = as.integer(
        round(
          dimension_ratio[k] *
            base_setting$n
        )
      ),
      
      p_Y = as.integer(
        round(
          dimension_ratio[k] *
            base_setting$n
        )
      ),
      
      comp = base_setting$comp,
      p_true = base_setting$p_true,
      grid_num = base_setting$grid_num,
      FVE = base_setting$FVE
    )
  }
)

task_grid <- expand.grid(
  case_id = seq_along(case_setting),
  replication = seq_len(n_rep),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)


################################################################################
# Iterative matrix factorization with diagnostic outputs
################################################################################

cca_AX_diagnostic <- function(
    lambda,
    tol,
    Cov,
    a_X,
    a_Y,
    group_x,
    group_y,
    max_iter = 200
) {
  
  error <- tol + 1
  ite <- 0L
  
  a_X <- as.numeric(a_X)
  a_Y <- as.numeric(a_Y)
  
  x_norm <- sqrt(
    sum(a_X^2)
  )
  
  y_norm <- sqrt(
    sum(a_Y^2)
  )
  
  degenerate <- FALSE
  
  
  while (
    is.finite(error) &&
    error > tol &&
    ite < max_iter
  ) {
    
    ite <- ite + 1L
    
    a_X_old <- a_X
    a_Y_old <- a_Y
    
    
    ##########################################################################
    # Update a_X
    ##########################################################################
    
    a_X <- as.numeric(
      Cov %*% a_Y
    )
    
    for (j in unique(group_x)) {
      
      index <- which(
        group_x == j
      )
      
      group_norm <- sqrt(
        sum(a_X[index]^2)
      )
      
      if (
        !is.finite(group_norm) ||
        group_norm == 0
      ) {
        
        a_X[index] <- 0
        
      } else {
        
        a_X[index] <-
          max(
            1 -
              lambda[1] /
              group_norm,
            0
          ) *
          a_X[index]
      }
    }
    
    x_norm <- sqrt(
      sum(a_X^2)
    )
    
    if (
      !is.finite(x_norm) ||
      x_norm == 0
    ) {
      
      degenerate <- TRUE
      error <- Inf
      
      break
    }
    
    a_X <- a_X / x_norm
    
    
    ##########################################################################
    # Update a_Y
    ##########################################################################
    
    a_Y <- as.numeric(
      t(a_X) %*% Cov
    )
    
    for (m in unique(group_y)) {
      
      index <- which(
        group_y == m
      )
      
      group_norm <- sqrt(
        sum(a_Y[index]^2)
      )
      
      if (
        !is.finite(group_norm) ||
        group_norm == 0
      ) {
        
        a_Y[index] <- 0
        
      } else {
        
        a_Y[index] <-
          max(
            1 -
              lambda[2] /
              group_norm,
            0
          ) *
          a_Y[index]
      }
    }
    
    y_norm <- sqrt(
      sum(a_Y^2)
    )
    
    if (
      !is.finite(y_norm) ||
      y_norm == 0
    ) {
      
      degenerate <- TRUE
      error <- Inf
      
      break
    }
    
    a_Y <- a_Y / y_norm
    
    
    ##########################################################################
    # Relative L1 error
    ##########################################################################
    
    denominator_X <- sum(
      base::abs(a_X_old)
    )
    
    denominator_Y <- sum(
      base::abs(a_Y_old)
    )
    
    if (
      is.finite(denominator_X) &&
      denominator_X > 0
    ) {
      
      error_X <-
        sum(
          base::abs(
            a_X -
              a_X_old
          )
        ) /
        denominator_X
      
    } else {
      
      error_X <- 0
    }
    
    if (
      is.finite(denominator_Y) &&
      denominator_Y > 0
    ) {
      
      error_Y <-
        sum(
          base::abs(
            a_Y -
              a_Y_old
          )
        ) /
        denominator_Y
      
    } else {
      
      error_Y <- 0
    }
    
    error <- error_X + error_Y
    
    if (!is.finite(error)) {
      
      degenerate <- TRUE
      error <- Inf
      
      break
    }
  }
  
  converged <-
    !degenerate &&
    is.finite(error) &&
    error <= tol
  
  list(
    X_coef = as.numeric(a_X),
    Y_coef = as.numeric(a_Y),
    iteration = ite,
    final_error = error,
    converged = converged,
    degenerate = degenerate
  )
}


################################################################################
# Construct the matrix used in the sparse FACD decomposition
################################################################################

prepare_FACD_matrix <- function(
    dat_X,
    dat_Y,
    grid_num,
    FVE
) {
  
  n <- dat_X$n
  p_X <- dat_X$p
  p_Y <- dat_Y$p
  
  optimal_basis <- Optimal_basis_est(
    dat_X = dat_X,
    dat_Y = dat_Y,
    grid_num = grid_num,
    max_comp_num = 5,
    mean_zero = FALSE,
    FVE = FVE
  )
  
  kappa_X <-
    optimal_basis$Optimal_bais$X$comp_num
  
  kappa_Y <-
    optimal_basis$Optimal_bais$Y$comp_num
  
  group_x <- rep(
    seq_len(p_X),
    each = kappa_X
  )
  
  group_y <- rep(
    seq_len(p_Y),
    each = kappa_Y
  )
  
  
  ##########################################################################
  # Project X data onto the estimated basis
  ##########################################################################
  
  dat_x_g <- t(
    sapply(
      seq_len(n),
      function(i) {
        
        c(
          sapply(
            seq_len(p_X),
            function(j) {
              
              t(
                optimal_basis$Optimal_bais$X$phi[
                  optimal_basis$Time_mark$X[[i]],
                  seq_len(kappa_X),
                  drop = FALSE
                ]
              ) %*%
                (
                  dat_X$dat[[j]][[i]]$val -
                    optimal_basis$Mean$X[[j]]$mu[
                      optimal_basis$Time_mark$X[[i]]
                    ]
                ) /
                length(
                  dat_X$dat[[j]][[i]]$time
                )
            }
          )
        )
      }
    )
  )
  
  
  ##########################################################################
  # Project Y data onto the estimated basis
  ##########################################################################
  
  dat_y_g <- t(
    sapply(
      seq_len(n),
      function(i) {
        
        c(
          sapply(
            seq_len(p_Y),
            function(j) {
              
              t(
                optimal_basis$Optimal_bais$Y$phi[
                  optimal_basis$Time_mark$Y[[i]],
                  seq_len(kappa_Y),
                  drop = FALSE
                ]
              ) %*%
                (
                  dat_Y$dat[[j]][[i]]$val -
                    optimal_basis$Mean$Y[[j]]$mu[
                      optimal_basis$Time_mark$Y[[i]]
                    ]
                ) /
                length(
                  dat_Y$dat[[j]][[i]]$time
                )
            }
          )
        )
      }
    )
  )
  
  
  ##########################################################################
  # Estimated finite-dimensional cross-covariance matrix
  ##########################################################################
  
  Cov <- crossprod(
    dat_x_g,
    dat_y_g
  ) / n
  
  if (any(!is.finite(Cov))) {
    
    stop(
      "The estimated covariance matrix contains non-finite values."
    )
  }
  
  list(
    Cov = Cov,
    group_x = group_x,
    group_y = group_y
  )
}


################################################################################
# Penalized objective
################################################################################

group_penalty <- function(
    coef,
    group
) {
  
  sum(
    sapply(
      unique(group),
      function(j) {
        
        sqrt(
          sum(
            coef[group == j]^2
          )
        )
      }
    )
  )
}


penalized_criterion <- function(
    Cov,
    a_X,
    a_Y,
    lambda,
    group_x,
    group_y
) {
  
  sum(
    (
      Cov -
        tcrossprod(
          a_X,
          a_Y
        )
    )^2
  ) +
    2 *
    lambda[1] *
    group_penalty(
      a_X,
      group_x
    ) +
    2 *
    lambda[2] *
    group_penalty(
      a_Y,
      group_y
    )
}


################################################################################
# Add a fixed-L2 random perturbation to an SVD initial vector
################################################################################

perturb_initial_vector <- function(
    initial_vector,
    perturbation_level
) {
  
  initial_vector <- as.numeric(
    initial_vector
  )
  
  initial_norm <- sqrt(
    sum(initial_vector^2)
  )
  
  if (
    !is.finite(initial_norm) ||
    initial_norm == 0
  ) {
    
    stop(
      "The supplied initial vector has zero or non-finite norm."
    )
  }
  
  initial_vector <-
    initial_vector /
    initial_norm
  
  noise <- rnorm(
    length(initial_vector)
  )
  
  noise_norm <- sqrt(
    sum(noise^2)
  )
  
  if (
    !is.finite(noise_norm) ||
    noise_norm == 0
  ) {
    
    stop(
      "The generated perturbation has zero or non-finite norm."
    )
  }
  
  # Convert Gaussian noise into a unit-L2 random direction
  noise <- noise / noise_norm
  
  # Add a perturbation with pre-normalization L2 magnitude
  # equal to perturbation_level
  perturbed_vector <-
    initial_vector +
    perturbation_level *
    noise
  
  perturbed_norm <- sqrt(
    sum(perturbed_vector^2)
  )
  
  if (
    !is.finite(perturbed_norm) ||
    perturbed_norm == 0
  ) {
    
    stop(
      "The perturbed initial vector has zero or non-finite norm."
    )
  }
  
  perturbed_vector /
    perturbed_norm
}


################################################################################
# Run one attempt of one simulation task
################################################################################

run_one_attempt <- function(
    task_index,
    attempt
) {
  
  case_id <-
    task_grid$case_id[task_index]
  
  replication <-
    task_grid$replication[task_index]
  
  setting <-
    case_setting[[case_id]]
  
  
  ############################################################################
  # Data seed
  #
  # attempt = 1 uses the original seed.
  # attempt = 2 uses a different seed.
  ############################################################################
  
  data_seed <-
    1000000 * case_id +
    replication +
    (attempt - 1L) *
    retry_seed_offset
  
  set.seed(
    data_seed
  )
  
  
  ############################################################################
  # Generate data
  ############################################################################
  
  dat <- data_gen(
    n = setting$n,
    rho = setting$rho,
    time_num = setting$time_num,
    p_X = setting$p_X,
    p_Y = setting$p_Y,
    comp = setting$comp,
    p_ture = setting$p_true
  )
  
  dat_X <- list(
    dat = dat$dat_X,
    n = dat$n,
    p = dat$p_X,
    feature_label = seq_len(dat$p_X),
    
    time_grid = seq(
      0,
      1,
      length.out = setting$time_num
    )
  )
  
  dat_Y <- list(
    dat = dat$dat_Y,
    n = dat$n,
    p = dat$p_Y,
    feature_label = seq_len(dat$p_Y),
    
    time_grid = seq(
      0,
      1,
      length.out = setting$time_num
    )
  )
  
  
  ############################################################################
  # Construct Gamma
  ############################################################################
  
  problem <- prepare_FACD_matrix(
    dat_X = dat_X,
    dat_Y = dat_Y,
    grid_num = setting$grid_num,
    FVE = setting$FVE
  )
  
  Cov <- problem$Cov
  group_x <- problem$group_x
  group_y <- problem$group_y
  
  
  ############################################################################
  # SVD initialization
  ############################################################################
  
  svd_fit <- svd(
    Cov,
    nu = 1,
    nv = 1
  )
  
  initial_X_svd <- as.numeric(
    svd_fit$u[, 1]
  )
  
  initial_Y_svd <- as.numeric(
    svd_fit$v[, 1]
  )
  
  
  ############################################################################
  # Runtime only for the SVD-initialized alternating factorization
  ############################################################################
  
  algorithm_time <- system.time({
    
    fit_svd <- cca_AX_diagnostic(
      lambda = lambda,
      tol = tol,
      Cov = Cov,
      a_X = initial_X_svd,
      a_Y = initial_Y_svd,
      group_x = group_x,
      group_y = group_y,
      max_iter = max_iter
    )
    
  })[["elapsed"]]
  
  
  ############################################################################
  # Rank-one solution under SVD initialization
  ############################################################################
  
  svd_solution_valid <-
    fit_svd$converged &&
    !fit_svd$degenerate &&
    all(is.finite(fit_svd$X_coef)) &&
    all(is.finite(fit_svd$Y_coef))
  
  if (svd_solution_valid) {
    
    rank_one_svd <- tcrossprod(
      fit_svd$X_coef,
      fit_svd$Y_coef
    )
    
    rank_one_svd_norm <- sqrt(
      sum(rank_one_svd^2)
    )
    
  } else {
    
    rank_one_svd <- NULL
    rank_one_svd_norm <- NA_real_
  }
  
  
  ############################################################################
  # Twenty perturbed SVD initializations
  ############################################################################
  
  perturbed_results <- lapply(
    seq_len(n_random_init),
    function(r) {
      
      perturbation_seed <-
        100000000 * case_id +
        100000 * replication +
        r +
        (attempt - 1L) *
        retry_seed_offset
      
      set.seed(
        perturbation_seed
      )
      
      initial_X_perturbed <-
        perturb_initial_vector(
          initial_vector = initial_X_svd,
          perturbation_level = perturbation_level
        )
      
      initial_Y_perturbed <-
        perturb_initial_vector(
          initial_vector = initial_Y_svd,
          perturbation_level = perturbation_level
        )
      
      fit_perturbed <- cca_AX_diagnostic(
        lambda = lambda,
        tol = tol,
        Cov = Cov,
        a_X = initial_X_perturbed,
        a_Y = initial_Y_perturbed,
        group_x = group_x,
        group_y = group_y,
        max_iter = max_iter
      )
      
      objective_perturbed <- penalized_criterion(
        Cov = Cov,
        a_X = fit_perturbed$X_coef,
        a_Y = fit_perturbed$Y_coef,
        lambda = lambda,
        group_x = group_x,
        group_y = group_y
      )
      
      list(
        fit = fit_perturbed,
        objective = objective_perturbed
      )
    }
  )
  
  
  ############################################################################
  # Select the valid perturbed result with the smallest objective
  ############################################################################
  
  valid_perturbed <- vapply(
    perturbed_results,
    function(x) {
      
      isTRUE(x$fit$converged) &&
        !isTRUE(x$fit$degenerate) &&
        is.finite(x$objective) &&
        all(is.finite(x$fit$X_coef)) &&
        all(is.finite(x$fit$Y_coef))
    },
    logical(1)
  )
  
  n_valid_perturbed <- sum(
    valid_perturbed
  )
  
  if (
    svd_solution_valid &&
    n_valid_perturbed > 0L
  ) {
    
    valid_results <-
      perturbed_results[
        valid_perturbed
      ]
    
    perturbed_objective <- vapply(
      valid_results,
      function(x) {
        
        x$objective
      },
      numeric(1)
    )
    
    best_perturbed_index <- which.min(
      perturbed_objective
    )
    
    fit_best_perturbed <-
      valid_results[[
        best_perturbed_index
      ]]$fit
    
    rank_one_best_perturbed <- tcrossprod(
      fit_best_perturbed$X_coef,
      fit_best_perturbed$Y_coef
    )
    
    rank_one_relative_difference <-
      sqrt(
        sum(
          (
            rank_one_svd -
              rank_one_best_perturbed
          )^2
        )
      ) /
      max(
        rank_one_svd_norm,
        .Machine$double.eps
      )
    
  } else {
    
    rank_one_relative_difference <-
      NA_real_
  }
  
  
  ############################################################################
  # Replication-level result
  ############################################################################
  
  data.frame(
    case_id = case_id,
    case = setting$case,
    dimension_ratio = setting$dimension_ratio,
    replication = replication,
    n = setting$n,
    p_X = setting$p_X,
    p_Y = setting$p_Y,
    
    task_success = TRUE,
    attempt_used = attempt,
    data_seed = data_seed,
    
    algorithm_runtime =
      as.numeric(algorithm_time),
    
    svd_iteration =
      fit_svd$iteration,
    
    svd_converged =
      fit_svd$converged,
    
    n_valid_perturbed =
      n_valid_perturbed,
    
    rank_one_relative_difference =
      rank_one_relative_difference,
    
    error_message =
      NA_character_,
    
    stringsAsFactors = FALSE
  )
}


################################################################################
# Run one task, with one automatic retry under a different seed
################################################################################

run_one_task <- function(
    task_index
) {
  
  case_id <-
    task_grid$case_id[task_index]
  
  replication <-
    task_grid$replication[task_index]
  
  setting <-
    case_setting[[case_id]]
  
  attempt_error <- character(
    max_task_attempts
  )
  
  for (attempt in seq_len(max_task_attempts)) {
    
    current_result <- tryCatch(
      {
        
        run_one_attempt(
          task_index = task_index,
          attempt = attempt
        )
      },
      
      error = function(e) {
        
        attempt_error[attempt] <<-
          conditionMessage(e)
        
        NULL
      }
    )
    
    if (!is.null(current_result)) {
      
      return(
        current_result
      )
    }
  }
  
  
  ############################################################################
  # Both attempts failed:
  # return an NA row rather than stopping all parallel workers
  ############################################################################
  
  final_seed <-
    1000000 * case_id +
    replication +
    (max_task_attempts - 1L) *
    retry_seed_offset
  
  data.frame(
    case_id = case_id,
    case = setting$case,
    dimension_ratio = setting$dimension_ratio,
    replication = replication,
    n = setting$n,
    p_X = setting$p_X,
    p_Y = setting$p_Y,
    
    task_success = FALSE,
    attempt_used = max_task_attempts,
    data_seed = final_seed,
    
    algorithm_runtime =
      NA_real_,
    
    svd_iteration =
      NA_integer_,
    
    svd_converged =
      FALSE,
    
    n_valid_perturbed =
      NA_integer_,
    
    rank_one_relative_difference =
      NA_real_,
    
    error_message =
      paste(
        paste0(
          "Attempt ",
          seq_len(max_task_attempts),
          ": ",
          attempt_error
        ),
        collapse = " | "
      ),
    
    stringsAsFactors = FALSE
  )
}


################################################################################
# Optional serial tests
################################################################################

# Test the first replication of the first case:
#
# test_result <- run_one_task(1L)
# print(test_result)


################################################################################
# Parallel computation
################################################################################

function_file <- normalizePath(
  "Function.R",
  mustWork = TRUE
)

cluster <- makeCluster(
  n_core
)

result <- tryCatch(
  {
    
    clusterExport(
      cluster,
      "function_file",
      envir = environment()
    )
    
    clusterEvalQ(
      cluster,
      {
        source(function_file)
        NULL
      }
    )
    
    clusterExport(
      cluster,
      c(
        "case_setting",
        "task_grid",
        "lambda",
        "tol",
        "max_iter",
        "n_random_init",
        "perturbation_level",
        "max_task_attempts",
        "retry_seed_offset",
        "cca_AX_diagnostic",
        "prepare_FACD_matrix",
        "group_penalty",
        "penalized_criterion",
        "perturb_initial_vector",
        "run_one_attempt",
        "run_one_task"
      ),
      envir = environment()
    )
    
    parLapplyLB(
      cluster,
      seq_len(nrow(task_grid)),
      run_one_task
    )
    
  },
  
  finally = {
    
    stopCluster(
      cluster
    )
  }
)


################################################################################
# Combine results
################################################################################

main_result <- do.call(
  rbind,
  result
)


################################################################################
# Summary helper
################################################################################

safe_mean <- function(x) {
  
  if (
    length(x) == 0L ||
    all(is.na(x))
  ) {
    
    return(
      NA_real_
    )
  }
  
  mean(
    x,
    na.rm = TRUE
  )
}


################################################################################
# Summarize the three cases
################################################################################

case_result <- split(
  main_result,
  main_result$case_id
)

summary_result <- do.call(
  rbind,
  lapply(
    case_result,
    function(dat) {
      
      successful_dat <- dat[
        dat$task_success,
        ,
        drop = FALSE
      ]
      
      converged_iteration <-
        successful_dat$svd_iteration[
          successful_dat$svd_converged
        ]
      
      data.frame(
        case =
          dat$case[1],
        
        n =
          dat$n[1],
        
        p_X =
          dat$p_X[1],
        
        p_Y =
          dat$p_Y[1],
        
        successful_replications =
          sum(dat$task_success),
        
        retried_replications =
          sum(
            dat$task_success &
              dat$attempt_used > 1L
          ),
        
        failed_replications =
          sum(!dat$task_success),
        
        runtime_seconds =
          safe_mean(
            successful_dat$algorithm_runtime
          ),
        
        average_iterations =
          safe_mean(
            converged_iteration
          ),
        
        rank_one_relative_difference =
          safe_mean(
            successful_dat$
              rank_one_relative_difference
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(summary_result) <- NULL


################################################################################
# Report tasks that failed under both seeds
################################################################################

failed_result <- main_result[
  !main_result$task_success,
  ,
  drop = FALSE
]

if (nrow(failed_result) > 0L) {
  
  warning(
    nrow(failed_result),
    " simulation task(s) failed under both the original and retry seeds. ",
    "These tasks were excluded from the numerical summaries."
  )
  
  cat(
    "\n",
    "============================================================\n",
    "Tasks That Failed Under Both Seeds\n",
    "============================================================\n",
    sep = ""
  )
  
  print(
    failed_result[
      ,
      c(
        "case",
        "replication",
        "data_seed",
        "error_message"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}


################################################################################
# Print computing environment
################################################################################

cat(
  "\n",
  "============================================================\n",
  "Computing Environment and Algorithm Settings\n",
  "============================================================\n",
  sep = ""
)

cat(
  "R version:",
  R.version.string,
  "\n"
)

cat(
  "Operating system:",
  Sys.info()[["sysname"]],
  "\n"
)

cat(
  "Machine:",
  Sys.info()[["machine"]],
  "\n"
)

cat(
  "Physical cores detected:",
  detectCores(logical = FALSE),
  "\n"
)

cat(
  "Parallel workers used:",
  n_core,
  "\n"
)

cat(
  "Requested replications per case:",
  n_rep,
  "\n"
)

cat(
  "Perturbed SVD initializations per replication:",
  n_random_init,
  "\n"
)

cat(
  "Relative L2 perturbation level:",
  perturbation_level,
  "\n"
)

cat(
  "Fixed lambda_X:",
  lambda[1],
  "\n"
)

cat(
  "Fixed lambda_Y:",
  lambda[2],
  "\n"
)

cat(
  "Stopping tolerance:",
  tol,
  "\n"
)

cat(
  "Maximum iterations:",
  max_iter,
  "\n"
)


################################################################################
# Draw the one-row boxplot figure
################################################################################

FACD_boxplot <- ggplot(
  boxplot_statistics,
  aes(
    x = case_ratio,
    ymin = ymin,
    lower = lower,
    middle = middle,
    upper = upper,
    ymax = ymax,
    fill = measure
  )
) +
  geom_boxplot(
    stat = "identity",
    width = 0.50,
    linewidth = 0.65,
    color = "#28323C",
    staplewidth = 0.42,
    show.legend = FALSE
  ) +
  facet_wrap(
    facets = vars(measure),
    nrow = 1,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "Runtime (seconds)" =
        "#4C78A8",
      
      "Iterations to convergence" =
        "#59A14F",
      
      "Rank-one relative difference (%)" =
        "#E39C37"
    )
  ) +
  
  # Put the y-axis ticks and labels on the left
  scale_y_continuous(
    position = "left",
    expand = expansion(
      mult = c(
        0.05,
        0.07
      )
    )
  ) +
  
  labs(
    x = expression(p / n),
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    panel.background = element_rect(
      fill = "#F5F7F9",
      color = NA
    ),
    
    panel.border = element_rect(
      fill = NA,
      color = "#AEB7C0",
      linewidth = 0.55
    ),
    
    panel.grid.major.x = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major.y = element_line(
      color = "#D8DDE2",
      linewidth = 0.35
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 12,
      margin = margin(
        t = 8
      )
    ),
    
    axis.text.x = element_text(
      color = "#28323C",
      size = 11
    ),
    
    axis.ticks.x = element_line(
      color = "#66717D",
      linewidth = 0.4
    ),
    
    # Display the y-axis on the left
    axis.text.y.left = element_text(
      color = "#28323C",
      size = 10,
      margin = margin(
        r = 5
      )
    ),
    
    axis.ticks.y.left = element_line(
      color = "#66717D",
      linewidth = 0.4
    ),
    
    axis.ticks.length.y.left = grid::unit(
      0.12,
      "cm"
    ),
    
    # Remove the right-side y-axis
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    
    strip.background = element_rect(
      fill = "#E4E9EE",
      color = "#AEB7C0",
      linewidth = 0.55
    ),
    
    strip.text = element_text(
      face = "bold",
      color = "#28323C",
      size = 11,
      margin = margin(
        t = 7,
        b = 7
      )
    ),
    
    panel.spacing.x = grid::unit(
      0.65,
      "lines"
    ),
    
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 12
    )
  )

print(FACD_boxplot)
ggsave(paste0("Figure/", "FACD_boxplot", ".pdf"), width = 10, height = 3.5, dpi = 300)
