# make a clustered sample
csample <- function(x,n,nclusters,maxdist){
  cpoints <- sp::spsample(x, n = nclusters, type="random")
  result <- cpoints
  result$clstrID <- 1:length(cpoints)
  for (i in 1:length(cpoints)){
    ext <- rgeos::gBuffer(cpoints[i,], width = maxdist)
    newsamples <- sp::spsample(ext, n = (n-nclusters)/nclusters, 
                               type="random")
    newsamples$clstrID <- rep(i,length(newsamples))
    result <- rbind(result,newsamples)
    
  }
  result$ID <- 1:nrow(result)
  return(result)
}

make_sample <- function(design, mask, npoints, nclusters, maxdist) {
  
  if (design=="random"){
    samplepoints <- spsample(mask,npoints,"random")
  }
  
  if (design=="clustered"){
    samplepoints <- csample(mask,npoints,nclusters,maxdist=maxdist)
  } 
  
  if (design=="biased"){
    countryboundaries <- getData("countries", path='../data/')
    countryboundaries <- countryboundaries[countryboundaries$NAME_ENGLISH%in%c(countries),]
    samplepoints <- spsample(countryboundaries,npoints,"random")
  }
  
  if (design=="biasedWithOutlier"){
    countryboundaries <- getData("countries", path='../data/')
    countryboundariesOut <- countryboundaries[countryboundaries$NAME_ENGLISH%in%c(countriesOutlier),]
    countryboundaries <- countryboundaries[countryboundaries$NAME_ENGLISH%in%c(countries),]
    samplepoints <- spsample(countryboundaries,npoints,"random")
    samplepoints <- rbind(samplepoints,spsample(countryboundariesOut,1,"random"))
  }  
  
  smppts_sf <- st_as_sf(samplepoints)
  return(smppts_sf)
}
