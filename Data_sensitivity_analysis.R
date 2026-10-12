################################################################################
# Time-point subsampling stability analysis
#
# Protein versus Metabolite
# First FACD component
################################################################################


################################################################################
# 1. Setup
################################################################################

rm(list = ls())

mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

source("Function.R")
source("LCCA_function.R")

library(parallel)
library(dplyr)
library(ggplot2)

dir.create("Result", showWarnings = FALSE, recursive = TRUE)
dir.create("Figure", showWarnings = FALSE, recursive = TRUE)


################################################################################
# 2. Settings
################################################################################

X_label <- "Protein"
Y_label <- "Metabolite"
output_name <- "Protein_Metabolomic"

component_index <- 1
n_rep <- 500
n_worker <- 40

subsample_probability <- 0.1
keep_proportion <- 0.80
minimum_time_points <- 3
base_seed <- 20260804


################################################################################
# 3. Load data and full-data result
################################################################################

load("Data/Protein.rda")
dat_X <- Protein
rm(Protein)

load("Data/Metabolomic.rda")
dat_Y <- Metabolomic
rm(Metabolomic)

load("Result/fit_Protein_Metabolomic.rda")
full_fit <- Result
rm(Result)

n <- dat_X$n
p_X <- dat_X$p
p_Y <- dat_Y$p

grid_num <- 201
max_var <- floor(c(min(p_X / 2, n), min(p_Y / 2, n)))


################################################################################
# 4. Full-data selected features
################################################################################

get_selected_indices <- function(selected_variables, p) {
  if(is.logical(selected_variables)){
    selected_variables <- which(selected_variables)
  }
  
  selected_variables <- as.integer(selected_variables)
  
  unique(
    selected_variables[
      is.finite(selected_variables) &
        selected_variables >= 1 &
        selected_variables <= p
    ]
  )
}


full_selected_X <- get_selected_indices(
  full_fit$var_X[[component_index]],
  p_X
)

full_selected_Y <- get_selected_indices(
  full_fit$var_Y[[component_index]],
  p_Y
)

full_load_X <- full_fit$load_X[[component_index]]
full_load_Y <- full_fit$load_Y[[component_index]]

feature_lookup <- dplyr::bind_rows(
  data.frame(
    omic = X_label,
    feature_index = full_selected_X,
    feature = as.character(dat_X$feature_label[full_selected_X])
  ),
  data.frame(
    omic = Y_label,
    feature_index = full_selected_Y,
    feature = as.character(dat_Y$feature_label[full_selected_Y])
  )
)


################################################################################
# 5. Common observation times
################################################################################

get_subject_times <- function(dat, subject_index) {
  for(j in seq_len(dat$p)){
    current_times <- dat$dat[[j]][[subject_index]]$time
    
    if(length(current_times) > 0){
      return(sort(unique(round(current_times, 10))))
    }
  }
  
  numeric(0)
}


common_time_list <- lapply(seq_len(n), function(i) {
  intersect(
    get_subject_times(dat_X, i),
    get_subject_times(dat_Y, i)
  )
})

if(any(lengths(common_time_list) < minimum_time_points)){
  stop(
    paste0(
      "At least one subject has fewer than ",
      minimum_time_points,
      " common observation times."
    )
  )
}


################################################################################
# 6. Time-point subsampling
################################################################################

filter_longitudinal_data <- function(dat, retained_times) {
  subsampled_dat <- dat
  
  for(i in seq_len(dat$n)){
    for(j in seq_len(dat$p)){
      current_curve <- dat$dat[[j]][[i]]
      
      retained_index <- round(current_curve$time, 10) %in%
        retained_times[[i]]
      
      subsampled_dat$dat[[j]][[i]] <-
        current_curve[retained_index, , drop = FALSE]
    }
  }
  
  subsampled_dat
}


subsample_time_points <- function(
    dat_X,
    dat_Y,
    common_time_list,
    subsample_probability,
    keep_proportion,
    minimum_time_points
) {
  retained_times <- lapply(common_time_list, function(current_times) {
    
    # Do not subsample this subject.
    if(stats::runif(1) >= subsample_probability){
      return(current_times)
    }
    
    number_retained <- max(
      minimum_time_points,
      floor(keep_proportion * length(current_times))
    )
    
    number_retained <- min(
      number_retained,
      length(current_times)
    )
    
    sort(
      sample(
        current_times,
        size = number_retained,
        replace = FALSE
      )
    )
  })
  
  list(
    dat_X = filter_longitudinal_data(
      dat_X,
      retained_times
    ),
    dat_Y = filter_longitudinal_data(
      dat_Y,
      retained_times
    )
  )
}


################################################################################
# 7. Feature-specific loading alignment
################################################################################

loading_alignment <- function(
    full_loading,
    subsampled_loading,
    eps = 1e-12
) {
  denominator <- sqrt(
    sum(full_loading^2) *
      sum(subsampled_loading^2)
  )
  
  if(!is.finite(denominator) || denominator <= eps){
    return(NA_real_)
  }
  
  min(
    abs(
      sum(full_loading * subsampled_loading) /
        denominator
    ),
    1
  )
}


construct_stability_result <- function(
    replication,
    omic_name,
    full_selected,
    subsample_selected,
    full_loading,
    subsample_loading
) {
  selected <- full_selected %in% subsample_selected
  alignment <- rep(NA_real_, length(full_selected))
  
  for(k in which(selected)){
    feature_index <- full_selected[k]
    
    alignment[k] <- loading_alignment(
      full_loading[feature_index, ],
      subsample_loading[feature_index, ]
    )
  }
  
  data.frame(
    replication = replication,
    omic = omic_name,
    feature_index = full_selected,
    selected = selected,
    alignment = alignment
  )
}


################################################################################
# 8. Run one subsampling repetition
################################################################################

run_one_attempt <- function(replication, seed_offset = 0) {
  set.seed(base_seed + replication + seed_offset)
  
  subsampled_data <- subsample_time_points(
    dat_X = dat_X,
    dat_Y = dat_Y,
    common_time_list = common_time_list,
    subsample_probability = subsample_probability,
    keep_proportion = keep_proportion,
    minimum_time_points = minimum_time_points
  )
  
  subsampled_fit <- FACD(
    grid_num = grid_num,
    dat_X = subsampled_data$dat_X,
    dat_Y = subsampled_data$dat_Y,
    max_var = max_var,
    length_tun = 30,
    print_tuning = FALSE,
    comp = 1,
    FVE = 0.95,
    Sparse = TRUE
  )
  
  subsample_selected_X <- get_selected_indices(
    subsampled_fit$var_X[[component_index]],
    p_X
  )
  
  subsample_selected_Y <- get_selected_indices(
    subsampled_fit$var_Y[[component_index]],
    p_Y
  )
  
  result_X <- construct_stability_result(
    replication = replication,
    omic_name = X_label,
    full_selected = full_selected_X,
    subsample_selected = subsample_selected_X,
    full_loading = full_load_X,
    subsample_loading = subsampled_fit$load_X[[component_index]]
  )
  
  result_Y <- construct_stability_result(
    replication = replication,
    omic_name = Y_label,
    full_selected = full_selected_Y,
    subsample_selected = subsample_selected_Y,
    full_loading = full_load_Y,
    subsample_loading = subsampled_fit$load_Y[[component_index]]
  )
  
  rbind(result_X, result_Y)
}


################################################################################
# 9. Retry once if FACD fails
################################################################################

run_one_replication <- function(replication) {
  result <- tryCatch(
    run_one_attempt(replication),
    error = function(e) NULL
  )
  
  if(!is.null(result)){
    return(result)
  }
  
  result <- tryCatch(
    run_one_attempt(
      replication,
      seed_offset = 500000000
    ),
    error = function(e) NULL
  )
  
  if(!is.null(result)){
    return(result)
  }
  
  rbind(
    data.frame(
      replication = replication,
      omic = X_label,
      feature_index = full_selected_X,
      selected = NA,
      alignment = NA_real_
    ),
    data.frame(
      replication = replication,
      omic = Y_label,
      feature_index = full_selected_Y,
      selected = NA,
      alignment = NA_real_
    )
  )
}


################################################################################
# 10. Parallel computation
################################################################################

cl <- parallel::makeCluster(
  min(n_worker, n_rep)
)

parallel::clusterExport(
  cl,
  "mydir",
  envir = .GlobalEnv
)

parallel::clusterEvalQ(cl, {
  setwd(mydir)
  
  source("Function.R")
  source("LCCA_function.R")
  
  library(fdapace)
  library(mvtnorm)
  library(fields)
  library(mgcv)
  library(fda)
})

parallel::clusterExport(
  cl,
  varlist = c(
    "dat_X",
    "dat_Y",
    "X_label",
    "Y_label",
    "p_X",
    "p_Y",
    "grid_num",
    "max_var",
    "component_index",
    "full_selected_X",
    "full_selected_Y",
    "full_load_X",
    "full_load_Y",
    "common_time_list",
    "subsample_probability",
    "keep_proportion",
    "minimum_time_points",
    "base_seed",
    "get_selected_indices",
    "filter_longitudinal_data",
    "subsample_time_points",
    "loading_alignment",
    "construct_stability_result",
    "run_one_attempt",
    "run_one_replication"
  ),
  envir = .GlobalEnv
)

stability_result <- dplyr::bind_rows(
  parallel::parLapply(
    cl,
    seq_len(n_rep),
    run_one_replication
  )
)

parallel::stopCluster(cl)


################################################################################
# 11. Stability summary
################################################################################

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  
  if(length(x) == 0){
    return(NA_real_)
  }
  
  unname(
    stats::quantile(
      x,
      probs = probability,
      names = FALSE
    )
  )
}


stability_summary <- stability_result %>%
  dplyr::group_by(
    omic,
    feature_index
  ) %>%
  dplyr::summarise(
    successful_repetitions = sum(!is.na(selected)),
    selected_repetitions = sum(selected %in% TRUE),
    number_of_alignments = sum(is.finite(alignment)),
    median_alignment = safe_quantile(alignment, 0.50),
    alignment_q1 = safe_quantile(alignment, 0.25),
    alignment_q3 = safe_quantile(alignment, 0.75),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    selection_probability = ifelse(
      successful_repetitions > 0,
      selected_repetitions / successful_repetitions,
      NA_real_
    )
  ) %>%
  dplyr::left_join(
    feature_lookup,
    by = c("omic", "feature_index")
  ) %>%
  dplyr::select(
    omic,
    feature_index,
    feature,
    successful_repetitions,
    selected_repetitions,
    selection_probability,
    number_of_alignments,
    median_alignment,
    alignment_q1,
    alignment_q3
  ) %>%
  dplyr::arrange(
    omic,
    dplyr::desc(selection_probability)
  )

write.csv(
  stability_summary,
  paste0(
    "Result/stability_summary_",
    output_name,
    ".csv"
  ),
  row.names = FALSE
)

print(stability_summary, n = Inf)


################################################################################
# 12. Prepare plotting data
################################################################################

# Original loading magnitude for each selected feature
loading_magnitude <- dplyr::bind_rows(
  data.frame(
    omic = X_label,
    feature_index = full_selected_X,
    loading_magnitude = sapply(
      full_selected_X,
      function(j) sqrt(sum(full_load_X[j, ]^2))
    )
  ),
  data.frame(
    omic = Y_label,
    feature_index = full_selected_Y,
    loading_magnitude = sapply(
      full_selected_Y,
      function(j) sqrt(sum(full_load_Y[j, ]^2))
    )
  )
)


# Keep the top 10 features within each omic
feature_information <- stability_summary %>%
  dplyr::left_join(
    loading_magnitude,
    by = c("omic", "feature_index")
  ) %>%
  dplyr::group_by(omic) %>%
  dplyr::slice_max(
    order_by = loading_magnitude,
    n = 10,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    feature_label = sprintf(
      "%s (SF = %.2f)",
      feature,
      selection_probability
    ),
    feature_panel = paste(
      omic,
      feature_index,
      sep = "___"
    )
  )


# Ascending factor levels make the largest loading appear at the top
feature_levels <- feature_information %>%
  dplyr::arrange(
    omic,
    loading_magnitude
  ) %>%
  dplyr::pull(feature_panel) %>%
  unique()


feature_axis_labels <- stats::setNames(
  feature_information$feature_label,
  feature_information$feature_panel
)


alignment_plot_data <- stability_result %>%
  dplyr::filter(
    selected %in% TRUE,
    is.finite(alignment)
  ) %>%
  dplyr::inner_join(
    feature_information %>%
      dplyr::select(
        omic,
        feature_index,
        feature_panel
      ),
    by = c(
      "omic",
      "feature_index"
    )
  ) %>%
  dplyr::mutate(
    omic = factor(
      omic,
      levels = c(X_label, Y_label)
    ),
    feature_panel = factor(
      feature_panel,
      levels = feature_levels
    )
  )


plot_backbone <- feature_information %>%
  dplyr::transmute(
    omic = factor(
      omic,
      levels = c(X_label, Y_label)
    ),
    feature_panel = factor(
      feature_panel,
      levels = feature_levels
    ),
    alignment = ifelse(
      omic == X_label,
      0.90,
      0.95
    )
  )


################################################################################
# 13. Feature-specific alignment boxplots
################################################################################

stability_plot <- ggplot2::ggplot() +
  ggplot2::geom_blank(
    data = plot_backbone,
    ggplot2::aes(
      x = alignment,
      y = feature_panel
    )
  ) +
  ggplot2::geom_boxplot(
    data = alignment_plot_data,
    ggplot2::aes(
      x = alignment,
      y = feature_panel,
      fill = omic
    ),
    orientation = "y",
    width = 0.60,
    linewidth = 0.50,
    outlier.size = 0.70
  ) +
  ggplot2::facet_wrap(
    ~ omic,
    nrow = 1,
    scales = "free"
  ) +
  ggh4x::facetted_pos_scales(
    x = list(
      omic == "Protein" ~ ggplot2::scale_x_continuous(
        limits = c(0.90, 1.00),
        breaks = seq(0.90, 1.00, by = 0.02),
        oob = scales::oob_keep,
        expand = ggplot2::expansion(
          mult = c(0.01, 0.02)
        )
      ),
      omic == "Metabolite" ~ ggplot2::scale_x_continuous(
        limits = c(0.95, 1.00),
        breaks = seq(0.95, 1.00, by = 0.01),
        oob = scales::oob_keep,
        expand = ggplot2::expansion(
          mult = c(0.01, 0.02)
        )
      )
    )
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) {
      unname(feature_axis_labels[x])
    }
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Protein" = "#3182BD",
      "Metabolite" = "#D95F0E"
    )
  ) +
  ggplot2::labs(
    x = "Feature loading alignment",
    y = NULL
  ) +
  ggplot2::theme_bw(
    base_family = "Times",
    base_size = 12
  ) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(linewidth = 0.55),
    strip.background = ggplot2::element_rect(
      fill = "#E6E9ED",
      color = "#9EA7B0",
      linewidth = 0.50
    ),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 12
    ),
    axis.title.x = ggplot2::element_text(face = "bold", size = 15),
    axis.text = ggplot2::element_text(color = "black"),
    axis.text.y = ggplot2::element_text(size = 10),
    legend.position = "none",
    panel.spacing.x = grid::unit(1.2, "lines")
  )

print(stability_plot)


################################################################################
# 14. Save figure
################################################################################

ggplot2::ggsave(
  filename = paste0(
    "Figure/stability_",
    output_name,
    ".pdf"
  ),
  plot = stability_plot,
  width = 12,
  height = 4,
  dpi = 300
)

