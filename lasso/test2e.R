## Lasso EW. Enet
ptm <- proc.time()
graphHippoLasso<-selectFast(extractedData,family="EW", G0=G0fMRI, K=1.1, gamma0=1)
proc.time() - ptm

graphHippoK1LassoG <- network(graphHippoLasso$EW$G)

pdf(file="pics/test2e.pdf")
plot(graphHippoK1LassoG, usearrows = FALSE, coord = loc2d,vertex.cex=0.5)
dev.off()
