# *****************************************************************************
# R Script implementing conventional random f-fold cross-validation.  
# Related to the paper "Dealing with clustered samples for assessing map 
# accuracy by cross-validation".
# Contact: Sytze de Bruin, Wageningen University, Laboratory of Geo-information
# Science and Remote Sensing, email: sytze.debruin@wur.nl
# May 3, 2022
# *****************************************************************************

# ****** load required library *******
library(ranger)
library(parallel)
source("/home/j/j_bahl03/R/CAST/R/global_validation.R")

# ************ GLOBALS ***************
samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")

infolder <- "~/deBruin_add_nndm/samples"
outfolder <- "~/emodi/CVresults"
startseed <- 1234567
n_CV   <- 3  # number of cross validation replications
n_samp <- 100  # number of sample replicates (for each design)
cores <- 10

# notes about the time:
# n_CV 10 x n_samp 20 with 20 cores in 1h: 3/5 samples
# n_CV  1 x n_samp 20 with 20 cores in 20m: all samples


# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/random")))
  dir.create(paste0(outfolder, "/random"))


# ************ FUNCTIONS ***************

sumSquares <- function(ref, pred){
  muref <- mean(ref, na.rm=T)
  SSR <- sum((ref - pred)^2)
  SST <- sum((ref - muref)^2)
  return(c(SSR, SST))
}

randomCV <- function(smpl, number, variate, seed){
  
  fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  MEC  <- numeric(n_CV)
  RMSE <- numeric(n_CV)
  
  for(i_CV in 1:n_CV){
    
    SSR <- 0
    SST <- 0
    folds <- list()
    
    # ************************
    # ****** 10 fold CV ******
    # ************************
    set.seed(seed)
    flds <- createFolds(AGBdata$agb, k = 10)
    
    model <- train(AGBdata[,!(names(AGBdata) %in% c("agb", "ID"))],
                   AGBdata$agb,
                   method="rf",
                   importance=TRUE,
                   tuneGrid = expand.grid(mtry=2),
                   trControl = trainControl(method="cv", index=flds,
                                            savePredictions = "final"))
    
    RMSE[i_CV] <- global_validation(model)["RMSE"][[1]]
    folds <- append(folds, flds)
    
  } # loop over i_CV
  
  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out  <- file.path(outfolder,"random", fname)
  save(RMSE, flds, file=f_out)
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    randomCV(smpl, i, "AGB", startseed)
    # randomCV(smpl, i, "OCS", startseed)
  }
}, mc.cores = cores)
