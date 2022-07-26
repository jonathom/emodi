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
source("/home/j/j_bahl03/R/CAST/R/CreateSpacetimeFolds.R")
source("/home/j/j_bahl03/R/CAST/R/global_validation.R")

# ************ GLOBALS ***************
samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")
infolder <- "~/investigate_spatial_validation/debruin/samples"
outfolder <- "~/investigate_spatial_validation/debruin/CVresults"
startseed <- 1234567
n_CV      <- 5  # number of cross validation replications
n_samp    <- 30  # # number of sample replicates (for each design)
cores <- 15

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
  
  RMSEs <- numeric(n_CV)
  
  for(i_CV in 1:n_CV) {
    # fo <- as.formula(paste0("agb~", paste(names(AGBdata)[-1], collapse = "+")))
    if (!("ID" %in% names(AGBdata))) {AGBdata$ID <- 1:nrow(AGBdata)}
    flds <- CreateSpacetimeFolds(AGBdata, spacevar = "ID", k = 10)
    # if ("ID" %in% names(AGBdata)) {AGBdata <- AGBdata[,!(names(AGBdata) %in% c("ID"))]}
    model <- train(AGBdata[,!(names(AGBdata) %in% c("agb", "ID"))],
                   AGBdata$agb,
                   method="rf",
                   importance=TRUE,
                   tuneGrid = expand.grid(mtry=2),
                   trControl = trainControl(method="cv",
                                            index=flds$index,
                                            indexOut=flds$indexOut,
                                            savePredictions = "final"))
    
    RMSEs[i_CV] <- global_validation(model)["RMSE"][[1]]
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
  save(RMSEs, file=f_out)
  
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    spatialCV(smpl, i, "AGB", startseed)
    # spatialCV(smpl, i, "OCS", startseed)
  }
}, mc.cores = cores)
