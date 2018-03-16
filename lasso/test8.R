## Enet ridge with EW

ptm <- proc.time()
mygraph<-selectFast(extractedData,family="EW", dmax = 30, G0=G0fMRI, K=1.1,gamma0=0, gamma=5)
proc.time() - ptm

mynetwork <- network(mygraph$EW$G)

pdf(file="pics/test8.pdf")
plot(mynetwork, usearrows = FALSE, coord = loc2d,vertex.cex=0.5)
dev.off()
