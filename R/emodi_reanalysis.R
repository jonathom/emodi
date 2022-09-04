.libPaths("/home/j/j_bahl03/R")
library(sf)
library(parallel)
library(dplyr)
library(ggplot2)
library(emdist)
source("~/emodi/R/plot_geodist_oldpipe.R")

discrete_curve <- function(gd, resolution=100) {
  sts <- gd$distances$dist[gd$distances$what == "sample-to-sample"]
  stp <- gd$distances$dist[gd$distances$what == "sample-to-prediction"]
  cv  <- gd$distances$dist[gd$distances$what == "CV-distances"]

  lower <- min(min(sts), min(stp), min(cv))
  upper <- max(max(sts), max(stp), max(cv))

  distance <- upper - lower
  step <- distance / (resolution - 1)
  # step*100 == distance

  l_sts <- length(sts)
  l_stp <- length(stp)
  l_cv  <- length(cv)

  fg <- data.frame(matrix(ncol=4, nrow=resolution))

  dists <- seq(from=lower, to=upper, by=step)
  for (i in 1:resolution) {
    fg[i,1] <- i # assign distance
    fg[i,2] <- length(sts[sts > dists[i] & sts <= dists[i+1]]) / l_sts
    fg[i,3] <- length(stp[stp > dists[i] & stp <= dists[i+1]]) / l_stp
    fg[i,4] <- length(cv[cv > dists[i] & cv <= dists[i+1]]) / l_cv
  }
  # fg <- na.omit(fg)

  names(fg) <- c("dist", "sample-to-sample", "sample-to-prediction", "CV-distances")

  return(fg)
}

sampled_geodist <- function(x, modeldomain, samples, cvfolds = NA, cv_method=TRUE, stat = 'density', showPlot = TRUE) {
  # x and modeldomain should be sf objects (points), folds should be a vector
  row_numbers <- sample(1:nrow(x), samples)
  row_numbers <- sort(row_numbers)
  if (!("ID" %in% names(x))) {x$ID <- 1:nrow(x)}

  sampled_x <- x[row_numbers,]
  sampled_modeldomain <- st_sample(st_as_sf(modeldomain), samples)
  sampled_modeldomain <- st_transform(sampled_modeldomain, st_crs(x))
  if(!is.na(cvfolds)[1]) {
    # this works because row_numbers is ordered
    l <- data.frame(row.names = 1:700)
    l[row_numbers,1] <- 1:10

    if(cv_method == "random" | cv_method == "spatial") {
      sampled_cvfolds <- lapply(cvfolds, function(x) {
        x <- l[x[x %in% row_numbers],]
      })
      gd <- plot_geodist(x = sampled_x, modeldomain = sampled_modeldomain, cvfolds = sampled_cvfolds, stat = stat, showPlot = showPlot)
    }
    if(cv_method == "nndm") {
      index_in <- cvfolds$idx_train[row_numbers]
      index_ex <- cvfolds$idx_test[row_numbers]
      index_in <- lapply(index_in, function(x) {
        x <- l[x[x %in% row_numbers],]
      })
      index_ex <- lapply(index_ex, function(x) {
        x <- l[x[x %in% row_numbers],]
      })
      gd <- plot_geodist(x = sampled_x, modeldomain = sampled_modeldomain, cvfolds = index_in, cvtrain = index_ex, stat = stat, showPlot = showPlot)
    }
  } else {
    gd <- plot_geodist(x = sampled_x, modeldomain = sampled_modeldomain, stat = stat, showPlot = showPlot)
  }
  return(gd)
}

EMD <- function(df, dist1, dist2) {
  emdist::emdw(A=df[,dist1], wA=rep(1,length(df)), B=df[,dist2], wB=rep(1,length(df)), dist = "euclidean")
}

root_path <- "~/emodi"
load(file.path(root_path, "/mask_preds_resp_europe.Rdata"))
# results_root <- "~/iloek_job/wadoux/deBruin_add_nndm/CVresults"
# samples_root <- "~/iloek_job/wadoux/deBruin_add_nndm/samples"
results_root <- "~/deBruin_add_nndm/CVresults"
samples_root <- "~/deBruin_add_nndm/samples"
samples <- list.files(samples_root)

geodist_per_sample <- function(method, smpl, iteration) {
  print(paste(method, smpl, iteration, "starting"))

  filename <- paste0("AGB_", smpl, sprintf("%03d", iteration), ".Rdata")
  coordsname <- paste0(sprintf("%03d", iteration), "_coords.Rdata")
  result_file <- file.path(root_path, "CVresults2", method, filename)
  coords_file <- file.path(samples_root, smpl, coordsname)
  exhaustive_file <- file.path(results_root, "exhaustive", filename)
  sample_file <- file.path(samples_root, smpl, paste0("AGBdata", sprintf("%03d", iteration), ".Rdata"))

  load(exhaustive_file) # this loads "RMSE"
  RMSE_val <- RMSE

  if (!file.exists(result_file)) next

  load(result_file) # this loads and overwrites RMSE
  if(length(RMSE) > 1) {RMSE <- mean(RMSE)}

  # df <- rbind(df, c(method, smpl, iteration, RMSE, RMSE_val, 0, 0, 0))
  # next

  load(coords_file)
  if(class(pts)[1] != "sf") {
    pts_df <- as.data.frame(pts)
    pts_sf <- st_as_sf(pts_df, coords = c("x", "y"))
    st_crs(pts_sf) <- st_crs("EPSG:3035") # set crs EPSG:3035
  } else {
    pts_sf <- pts
  }
  mask <- st_transform(st_as_sf(mask), st_crs(pts_sf))

  s2s_s2p <- numeric(3)
  s2s_cv <- numeric(3)
  s2p_cv <- numeric(3)

  for (i in 1:3) {
    if (method == "spatial") {
      fold_index <- i*2 - 1
      flds_i <- folds[fold_index][[1]]
    }
    if (method == "random") {
      fold_indices <- 1:10 + i*10 - 10
      flds_i <- folds[fold_indices]
    }
    if (method == "nndm") {
      fold_index <- i*7 - 6
      idx_train <- folds[fold_index][[1]]
      idx_test <- folds[fold_index + 1][[1]]
      flds_i <- list("idx_train" = idx_train, "idx_test" = idx_test)
    }
    gd <- sampled_geodist(x = pts_sf, modeldomain = mask, cvfolds = flds_i,
                          stat='density', samples = 200, showPlot = FALSE, cv_method = method)

    # s2s_s2p[i] <- EMD_s2s_s2p(gd)
    # s2s_cv[i] <- EMD_s2s_cv(gd)
    # s2p_cv[i] <- EMD_s2p_cv(gd)

    gd_discrete <- discrete_curve(gd)
    s2s_s2p[i] <- EMD(gd_discrete, "sample-to-sample", "sample-to-prediction")
    s2s_cv[i] <- EMD(gd_discrete, "sample-to-sample", "CV-distances")
    s2p_cv[i] <- EMD(gd_discrete, "sample-to-prediction", "CV-distances")

  }

  s2s_s2p <- mean(s2s_s2p)
  s2s_cv <- mean(s2s_cv)
  s2p_cv <- mean(s2p_cv)

  res <- c(method, smpl, iteration, RMSE, RMSE_val, s2s_s2p, s2s_cv, s2p_cv)
  fname <- paste0(method,"_",smpl,"_",iteration,".Rdata")
  save(res, file=paste0("~/emodi/CVresults2/reanalysis/", fname))
  return(c(method, smpl, iteration, RMSE, RMSE_val, s2s_s2p, s2s_cv, s2p_cv))
}



# df <- data.frame(method=NA, sample=NA, iteration=NA, RMSE=NA, RMSE_val=NA, s2s_s2p=NA, s2s_cv=NA, s2p_cv=NA)

do_list <- list()
x <- 1
for (method in c("nndm", "spatial", "random")) {
  for (smpl in samples) {
    for (iteration in 1:100) {
      do_list[[x]] <- list("method" = method, "smpl" = smpl, "iteration" = iteration)
      x <- x+1
    }
  }
}

all_rows <- lapply(do_list, function(el) {
  geodist_per_sample(el$method, el$smpl, el$iteration)
})

df_palma <- as.data.frame(do.call(rbind,all_rows))
names(df_palma) <- c("method", "sample", "iteration", "RMSE", "RMSE_val", "s2s_s2p", "s2s_cv", "s2p_cv")
df_palma$RMSE <- as.numeric(df_palma$RMSE)
df_palma$RMSE_val <- as.numeric(df_palma$RMSE_val)
df_palma$s2s_s2p <- as.numeric(df_palma$s2s_s2p)
df_palma$s2s_cv <- as.numeric(df_palma$s2s_cv)
df_palma$s2p_cv <- as.numeric(df_palma$s2p_cv)
out_file <- file.path("~/emodi/result_df_palma_200_samples_cast_nndm_reanalysis.Rdata")
save(df_palma, file = out_file)
