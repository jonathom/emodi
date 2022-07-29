# *****************************************************************************
# R Script implementing conventional random cross-validation.  
# Related to the manuscript "Dealing with clustered samples for assessing map 
# accuracy by cross-validation".
# Contact: Sytze de Bruin, Wageningen University, Laboratory of Geo-information
# Science and Remote Sensing, email: sytze.debruin@wur.nl
# May 3, 2022
# *****************************************************************************

# ****** load required libraries *******
library(ranger)
# library(sperrorest, lib.loc="/home/j/j_bahl03/R")
library(parallel)
# library(CAST)
library(caret)
library(sf)
source("/home/j/j_bahl03/R/CAST/R/CreateSpacetimeFolds.R")
source("/home/j/j_bahl03/R/CAST/R/global_validation.R")

# ************ GLOBALS ***************
samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")
infolder <- "~/deBruin_add_nndm/samples"
# infolder <- ("../../deBruin_add_nndm/samples")
outfolder <- "~/emodi/CVresults"
startseed <- 1234567
n_CV      <- 3  # number of cross validation replications
n_samp    <- 100  # # number of sample replicates (for each design)
cores <- 10

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/spatial")))
  dir.create(paste0(outfolder, "/spatial"))


# ************ FUNCTIONS ***************

err_fu <- function(obs, pred){
  rmse <- sqrt(mean((obs-pred)^2))
  muref <- mean(obs)
  SSR <- sum((obs - pred)^2)
  SST <- sum((obs - muref)^2)
  mec <- 1 - SSR/SST
  me <- mean(obs - pred)
  list(me = me, rmse = rmse, mec = mec)
}


predfun <- function(object, newdata){
  pred <- predict(object, newdata)
  pred[[1]]
}


spatialCV <- function(smpl, number, variate, seed){
  
  fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  RMSE <- numeric(n_CV)
  folds <- list()
  
  for(i_CV in 1:n_CV) {
    # fo <- as.formula(paste0("agb~", paste(names(AGBdata)[-1], collapse = "+")))
    # if (!("ID" %in% names(AGBdata))) {AGBdata$ID <- 1:nrow(AGBdata)}
    # flds <- CreateSpacetimeFolds(AGBdata, spacevar = "ID", k = 10)
    # if ("ID" %in% names(AGBdata)) {AGBdata <- AGBdata[,!(names(AGBdata) %in% c("ID"))]}
    
    # find out autocorrelation ######################
    # stack1 <- raster::raster("~/iloek_job/wadoux/investigate_spatial_validation/debruin/data/agb.tif")
    # s_df <- raster::as.data.frame(stack1, xy=T, na.rm=T)
    # sp_s <- s_df[sample(nrow(s_df), 10000), ]
    # coordinates(sp_s) <- ~x+y
    # lzn.vgm = gstat::variogram(agb~1, data = sp_s, cutoff = 700000)
    # plot(lzn.vgm$dist, lzn.vgm$gamma)
    # -> 500000 # pixelsize 500m, -> 1000 pixels autocorrelation
    #################################################
    
    # wadoux spatial CV
    # mdist <- dist(AGBdata[c("xcoord","ycoord")])
    # hc <- hclust(mdist, method="complete")
    # d = 1000  # the maximum distance between pixels within clusters (val.dist m * 1000 m = val.dist km)                   
    # AGBdata$Clust_val.distkm = cutree(hc, h=d) 
    
    # try with coords
    load(file.path(infolder, smpl, paste0(sprintf("%03d", number), "_coords.Rdata")))
    
    if(class(pts)[1] != "sf") {
      pts <- as.data.frame(pts)
      pts <- st_as_sf(pts, coords = c("x", "y"))
    }
    
    mdist2 <- dist(st_coordinates(pts))
    hc2 <- hclust(mdist2, method="complete")
    AGBdata$Clust_val.distkm = cutree(hc2, h=500000) # tune this, leads to many folds!

    flds <- CreateSpacetimeFolds(AGBdata, spacevar = "Clust_val.distkm")
    
    model <- train(AGBdata[,!(names(AGBdata) %in% c("agb", "ID", "Clust_val.distkm"))],
                   AGBdata$agb,
                   method="rf",
                   importance=TRUE,
                   tuneGrid = expand.grid(mtry=2),
                   trControl = trainControl(method="cv",
                                            index=flds$index,
                                            indexOut=flds$indexOut,
                                            savePredictions = "final"))
    
    RMSE[i_CV] <- global_validation(model)["RMSE"][[1]]
    folds <- append(folds, flds)
  }
  
  # if(variate == "AGB"){
  #   fo <- as.formula(paste0("agb~", paste(names(AGBdata)[-1], collapse = "+")))
  #   tst <- sperrorest(fo, data=AGBdata, model_fun=ranger, 
  #                     model_args=list(respect.unordered.factors=TRUE),
  #                     pred_fun=predfun, 
  #                     smp_fun=partition_kmeans, coords=c("xcoord", "ycoord"),
  #                     smp_args=list(balancing_steps = 1, seed1=seed, 
  #                                   repetition=1:n_CV, iter.max = 50), 
  #                     err_fun = err_fu)
  # } else{
  #   fo <- as.formula(paste0("ocs~", paste(names(OCSdata)[-1], collapse = "+")))
  #   tst <- sperrorest(fo, data=OCSdata, model_fun=ranger, 
  #                     model_args=list(respect.unordered.factors=TRUE),
  #                     pred_fun=predfun, 
  #                     smp_fun=partition_kmeans, coords=c("xcoord", "ycoord"),
  #                     smp_args=list(balancing_steps = 1, seed1=seed, 
  #                                   repetition=1:n_CV, iter.max = 50), 
  #                     err_fun = err_fu)
  # }

  #ME   <- tst$error_rep$test_me
  #RMSE <- tst$error_rep$test_rmse
  #MEC  <- tst$error_rep$test_mec
  #rm(tst)
  
  fname <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out <- file.path(outfolder, "spatial", fname)
  save(RMSE, folds, file=f_out)
  
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    spatialCV(smpl, i, "AGB", startseed)
    # spatialCV(smpl, i, "OCS", startseed)
  }
}, mc.cores = cores)
