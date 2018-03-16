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
extTSLength <- 50 ## num of time points extracted for each subject.
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

mydmax = 30
