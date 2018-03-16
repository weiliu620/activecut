testID = "test2a"
extTSLength <- 50 ## num of time points extracted for each subject.
mydmax = 30
myK = 1.01
myEps = 3

library(GGMselect)
dyn.load("loopG0.so")
source("penaltyG0.R")
source("convertGraph.R")

# Chargement des donnees

load("fmriData.Rdata")
data <- fmriData$pcaObj$x
tsLength <- fmriData$tsLength # time series length
nSub <- fmriData$nSub # num of subjects
nNodes <- fmriData$nNodes # num of nodes in each

# extract some data
extractedData = matrix(0, nrow = nNodes, ncol = 0)

for ( s in 1:nSub) {
  extractedData <- cbind(extractedData, data[((s-1)*nNodes+1):(s*nNodes), 1:extTSLength])
}
extractedData = t(extractedData)

## loc3D <- matrix(scan("Loc3D.txt", 0), ncol=3, byrow=TRUE)
loc3d <- fmriData$linear2sub
loc2d <- loc3d[,1:2]
  
# Read in G0
load("G0fMRI.Rdata")

source("GGMselectG0.R")
source("ModLassoG0.R")

source("MyFamilyG0.R")
library(network)

source("GGMselectEnetG0.R")
source("ModLassoEnetG0.R")

## define a G0 matrix.
source("defineG0.R")
myG0 = defineG0(fmriData$linear2sub, myEps)

## Lasso EW. L1 penalty with light constrain.
mygraph<-selectFast(extractedData,family="EW", dmax = mydmax, G0=myG0, K=myK, gamma0=0)
mynetwork <- network(mygraph$EW$G)

pdf(file=paste("pics/", testID, ".pdf", sep = ""))
plot(mynetwork, usearrows = FALSE, coord = loc2d,vertex.cex=0.5)
dev.off()

## extract roi network.

roiName = colnames(fmriData$roiMask)
for (roiIdx in 1:length(roiName)) {
  troiMask = fmriData$roiMask[,roiIdx]
  tG = mygraph$EW$G
  tG[troiMask == 0,] = 0
  tG = as.numeric(tG | t(tG))
  tG = matrix(tG, nrow = nrow(fmriData$roiMask), ncol = nrow(fmriData$roiMask))
  tG1 <- network(tG)
  filename = paste("pics/", testID, "_", roiName[roiIdx], sep= "")
  pdf(file= filename, title = filename)
  plot(tG1, usearrows = FALSE, coord = loc2d,vertex.cex=0.5)
  dev.off()
}
