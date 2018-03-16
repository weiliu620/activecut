## ridge and enet.
ptm <- proc.time()
graphHippoLasso<-selectFast(extractedData,family="LA", G0=G0fMRI, K=1.1, gamma0=0.005, gamma=1)
proc.time() - ptm

graphHippoK1LassoG <- network(graphHippoLasso$EW$G)

pdf(file="pics/test3e.pdf")
plot(graphHippoK1LassoG, usearrows = FALSE, coord = loc2d,vertex.cex=0.5)
dev.off()
