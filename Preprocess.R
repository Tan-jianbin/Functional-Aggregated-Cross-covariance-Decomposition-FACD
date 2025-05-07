rm(list=ls())
# set working directory to be where the current script is located
mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

###############################################################################################
# Wrapping function
wrap.dat <- function(dat, dat_count, time){
  
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
      time <- sapply(1:length(time), function(k) time[time[k]])
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

###############################################################################################
dat_all <- readRDS("Data/exercise_log10_multiomic.RDS") # Original data

# Protein data
dat <- cbind(dat_all$cytokine, dat_all$protein) # Original data
dat_count <- dat_all$meta_visit # Label of subjects and times
time <- 0:4 / 4 # Observational times

Protein <- wrap.dat(dat = dat, dat_count = dat_count, time = time)
save(Protein, file = "Data/Protein.rda", version = 2)

###############################################################################################
# metabolite
dat <- dat_all$metabolome # Original data

Metabolomic <- wrap.dat(dat = dat, dat_count = dat_count, time = time)

save(Metabolomic, file = "Data/Metabolomic.rda", version = 2)

###############################################################################################
# Lipid
dat <- dat_all$lipid # Original data

Lipid <- wrap.dat(dat = dat, dat_count = dat_count, time = time)

save(Lipid, file = "Data/Lipid.rda", version = 2)

