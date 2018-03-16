## Lasso EW with ridge.
ptm <- proc.time()
graphHippoLasso<-selectFast(extractedData,family="EW", dmax = mydmax, G0=G0fMRI, K=1.1, gamma0=0.1)
proc.time() - ptm

graphHippoK1LassoG <- network(graphHippoLasso$EW$G)

pdf(file="pics/test4.pdf")
plot(graphHippoK1LassoG, usearrows = FALSE, coord = loc2d,vertex.cex=0.5)
dev.off()
