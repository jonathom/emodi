samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")
n_samp <- 1:100  # number of sample replicates (for each design)
# 20000 phi and 0.5 min train

l <- data.frame(sample=NA, i_samp=NA, variate=NA)
x <- 1
for (smpl in samples) {
  for (i in n_samp) {
    for (var in c("AGB")) {
      l <- rbind(l, c(smpl, i, var))
      x <- x + 1
    }
  }
}

l <- l[2:nrow(l),]
l$job <- 0
write.csv(l, "~/iloek_job/wadoux/emodi/CVresults/nndm/nndm_processing.csv")
save(l, file="~/iloek_job/wadoux/investigate_spatial_validation/debruin/R/nndm_processing.Rdata")
