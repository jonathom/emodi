# *****************************************************************************
# R Script implementing conventional random f-fold cross-validation.  
# Related to the paper "Dealing with clustered samples for assessing map 
# accuracy by cross-validation".
# Contact: Sytze de Bruin, Wageningen University, Laboratory of Geo-information
# Science and Remote Sensing, email: sytze.debruin@wur.nl
# May 3, 2022
# *****************************************************************************

# sleep for random time so multiple jobs dont take the same job out of the csv
Sys.sleep(round(runif(1, min = 1, max = 240)))

# ****** load required library *******
library(ranger)
# library(NNDM)
source("~/R/NNDM/R/nndm.R")
library(sf)
library(raster)
library(caret)
library(parallel)

# ************ GLOBALS ***************
infolder <- "~/investigate_spatial_validation/debruin/samples"
outfolder <- "~/investigate_spatial_validation/debruin/CVresults"
# outfolder <- "~/iloek_job/wadoux/investigate_spatial_validation/debruin/CVresults"
datafolder <- "~/investigate_spatial_validation/debruin/data"
folder_name <- "nndm_single"

csv_file <- file.path(outfolder, folder_name, "nndm_processing.csv")
runs <- read.csv(csv_file)
i <- 0
found <- FALSE
while (found == FALSE && i < 301) {
  i <- i+1
  row <- runs[i,]
  if (row$job == 0) {
    # print(paste(row, i))
    row$job <- 1
    runs[i,] <- row
    found <- TRUE
  }
}

write.csv(runs, file = csv_file, row.names = FALSE)

# 20000 phi and 0.5 min train
n <- 5000 # usually 5000
n_CV <- 3

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, folder_name)))
  dir.create(paste0(outfolder, folder_name))

# ************ FUNCTIONS ***************

sumSquares <- function(ref, pred){
  muref <- mean(ref, na.rm=T)
  SSR <- sum((ref - pred)^2)
  SST <- sum((ref - muref)^2)
  return(c(SSR, SST))
}

err_fu <- function(obs, pred){
  rmse <- sqrt(mean((obs-pred)^2))
  muref <- mean(obs)
  SSR <- sum((obs - pred)^2)
  SST <- sum((obs - muref)^2)
  mec <- 1 - SSR/SST
  me <- mean(obs - pred)
  list(me = me, rmse = rmse, mec = mec)
}

nndmCV <- function(smpl, number, variate) {
  
  fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  MEC  <- numeric(n_CV)
  RMSE <- numeric(n_CV)
  
  # load sample file containing coordinates
  load(file.path(infolder,smpl,paste0(sprintf("%03d", number), "_coords", ".Rdata")))
  if(class(pts)[1] != "sf") {
    sample_df <- as.data.frame(pts)
    sample_sf <- st_as_sf(sample_df, coords = c("x", "y"))
  } else {
    sample_sf <- pts
  }
  
  agb_raster <- raster::raster(file.path(datafolder, "agb.tif")) # load agb raster

  for(i_CV in 1:n_CV) {
  
    raster_subset <- sampleRandom(agb_raster, n, sp=T) # subset raster
    raster_sf <- st_as_sf(raster_subset) # as sf
    sample_subset <- st_cast(st_sample(sample_sf, n), to = "POINT") # subset sample
    st_crs(sample_subset) <- st_crs(raster_sf) # set crs
    
    #st_as_sf(raster::rasterToPoints(agb_raster[[1]], spatial = TRUE))
    nndm <- nndm(tpoints = sample_subset, ppoints = raster_sf, phi = 20000, min_train = 0.5)
    # save(nndm, file="./nndm.Rdata")
    
    # Evaluate RF model using NDM CV
    trainControl_NNDM <- trainControl(method = "cv",
                                      index=nndm$indx_train,
                                      indexOut=nndm$indx_test,
                                      savePredictions = "final") # save predictions final to avoid writing CV myself
    
    paramGrid <-  data.frame(mtry = 2, min.node.size = 5, splitrule = "variance")
    
    if (variate == "AGB") {training_data <- AGBdata; rf_form <- agb~.}
    else {training_data <- OCSdata; rf_form <- ocs~.}
    
    # next: execute this and see where it goes from there
    mod_NNDM <- train(rf_form,
                      method = "ranger",
                      trControl = trainControl_NNDM,
                      tuneGrid = paramGrid, 
                      data = training_data)
    
    refs <- mod_NNDM$pred$obs
    preds <- mod_NNDM$pred$pred
    
    RMSE[i_CV] <- err_fu(refs, preds)["rmse"][[1]]
    MEC[i_CV] <- err_fu(refs, preds)["mec"][[1]]
  }

  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  fname2 <-  paste0("nndm", variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out  <- file.path(outfolder,folder_name, fname)
  f_out2 <- file.path(outfolder,folder_name, fname2)
  save(MEC, RMSE, file=f_out)
  save(nndm, file=f_out2)
}

# ************ CALL THE FUNCTIONS ************ 
nndmCV(smpl = row$sample, number = row$i_samp, variate = row$variate)
