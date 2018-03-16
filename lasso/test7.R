## Enet Ridge regularizaiton. LA

ptm <- proc.time()
mygraph<-selectFast(extractedData,family="LA", dmax = 30, G0=G0fMRI, K=1.1,gamma0=0, gamma=5)
proc.time() - ptm

mynetwork <- network(mygraph$LA$G)

pdf(file="pics/test7.pdf")
plot(mynetwork, usearrows = FALSE, coord = loc2d,vertex.cex=0.5)
dev.off()
