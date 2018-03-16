## Lasso LA. Enet
ptm <- proc.time()
graphHippoLasso<-selectFast(extractedData,family="LA", G0=G0fMRI, K=1.1,gamma0=0, gamma=1)
proc.time() - ptm

graphHippoK1LassoG <- network(graphHippoLasso$LA$G)

pdf(file="pics/test1e.pdf")
plot(graphHippoK1LassoG, usearrows = FALSE, coord = loc2d,vertex.cex=0.5)
dev.off()
