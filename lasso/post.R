postprocess <- function(G, loc2d, roiMask, testID)
  ## given graph and roiMask, find the connections only for this
  ## roi. roiMask can be a PxK matrix. P is number of voxels and K is
  ## number of roi's. testID is like "test2".
  {
    roiName = colnames(roiMask)
    for (roiIdx in 1:length(roiName)) {
      troiMask = roiMask[,roiIdx]
      tG = G
      tG[troiMask == 0,] = 0
      tG = as.numeric(tG | t(tG))
      tG = matrix(tG, nrow = nrow(roiMask), ncol = nrow(roiMask))
      tG1 <- network(tG)
      filename = paste("pics/", testID, "_", roiName[roiIdx], sep= "")
      pdf(file= filename, title = filename)
      plot(tG1, usearrows = FALSE, coord = loc2d,vertex.cex=0.5)
      dev.off()
    }
  }
