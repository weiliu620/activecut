## Lasso EW. L1 penalty with light constrain.
ptm <- proc.time()
graphHippoLasso<-selectFast(extractedData,family="LA", dmax = 30, G0=G0fMRI, K=1.1,gamma0=0)
proc.time() - ptm

graphHippoK1LassoG <- network(graphHippoLasso$LA$G)

pdf(file="pics/test1.pdf")
plot(graphHippoK1LassoG, usearrows = FALSE, coord = loc2d,vertex.cex=0.5)
dev.off()
