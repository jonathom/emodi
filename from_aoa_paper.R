# from https://github.com/HannaMeyer/MEE_AOA

# preparations

rm(list=ls())
#install_github("HannaMeyer/CAST")
library(virtualspecies)
library(caret)
library(CAST)
library(viridis)
library(gridExtra)
library(knitr)
library(grid)
library(latticeExtra)

rmse <- function(pred,obs){sqrt( mean((pred - obs)^2, na.rm = TRUE) )}

npoints <- 50 # number of training samples

design <- "random"  # either clustered,biased,random, biasedWithOutlier

countries <- c("Germany","Ireland","France", "Sweden") #if design==biased
countriesOutlier <- "Turkmenistan" #if design==biasedWithOutlier A single point is set here
nclusters <- 50 #number of clusters if design==clustered
maxdist <- 0.8 #maxdist for clustered samples if design==clustered
meansPCA <- c(3, -1) # means of the gaussian response functions to the 2 axes
sdPCA <- c(2, 2) # sd's of the gaussian response functions to the 2 axes
simulateResponse <- c("bio2","bio5","bio10", "bio13",
                      "bio14","bio19") # variables used to simulate the response
studyarea <- c(-15, 65, 30, 75) # extent of study area. Default: Europe
seed <- 10

# download data
predictors_global <- raster::getData('worldclim', var='bio', res=10, path='./data/')
wp <- extent(studyarea)
predictors <- crop(predictors_global,wp)

# create a mask for land area:
land_mask <- predictors[[1]]
values(land_mask)[!is.na(values(land_mask))] <- 1

# generate predictor / response
response_vs <- generateSpFromPCA(predictors[[simulateResponse]],
                                 means = meansPCA,sds = sdPCA, plot=F)
response <- response_vs$suitab.raster




mask <- rasterToPolygons(land_mask,dissolve=TRUE)

save(mask, predictors, response, file="./mask_preds_resp_europe.Rdata")


set.seed(seed)




p1 <- spplot(stretch(predictors[[simulateResponse]],0,1),col.regions=viridis(100),
             par.settings =list(strip.background=list(col="grey")))

p2 <-spplot(response,col.regions=viridis(100),
            sp.layout=list("sp.points", samplepoints, col = "red", first = FALSE))
grid.arrange(p1,p2,ncol=2,nrow=1,widths=c(1.25,1))

grid.text("a",x = unit(0.02, "npc"), 
          y = unit(0.96, "npc"),
          just = "left")
grid.text("b",x = unit(0.57, "npc"), 
          y = unit(0.96, "npc"),
          just = "left")

trainDat <- extract(predictors,samplepoints,df=TRUE)
trainDat$response <- extract (response,samplepoints)

if (design=="clustered"){
  trainDat <- merge(trainDat,samplepoints,by.x="ID",by.y="ID")
}

trainDat <- trainDat[complete.cases(trainDat),]


set.seed(seed)
if(design!="clustered"){
  model <- train(trainDat[,names(predictors)],
                 trainDat$response,
                 method="rf",
                 importance=TRUE,
                 tuneGrid = expand.grid(mtry = c(2:length(names(predictors)))),
                 trControl = trainControl(method="cv",savePredictions = TRUE))
  
  
  print("random cross-validation performance:")
  print(model)
}
#if data are clustered, clustered CV is used:
if(design=="clustered"){
  folds <- CreateSpacetimeFolds(trainDat, spacevar="clstrID",k=nclusters)
  model <- train(trainDat[,names(predictors)],
                 trainDat$response,
                 method="rf",
                 importance=TRUE,
                 tuneGrid = expand.grid(mtry = c(2:length(names(predictors)))),
                 trControl = trainControl(method="cv",index=folds$index,savePredictions = TRUE))
  print("leave-cluster-out cross-validation performance:")
  print(model)
  
  model_random <- train(trainDat[,names(predictors)],
                        trainDat$response,
                        method="rf",
                        importance=TRUE,
                        tuneGrid = expand.grid(mtry = c(2:length(names(predictors)))),
                        trControl = trainControl(method="cv",savePredictions = TRUE))
  
  print("random cross-validation performance:")
  print(model_random)
}

prediction <- predict(predictors,model)
truediff <- abs(prediction-response)