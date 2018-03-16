library(GGMselect)
dyn.load("loopG0.so")
source("penaltyG0.R")
source("convertGraph.R")

# Chargement des donnees
data <- matrix(scan("data.txt", 0), ncol=101, byrow=TRUE)
loc3D <- matrix(scan("Loc3D.txt", 0), ncol=3, byrow=TRUE)
loc2D <- loc3D[,1:2]

sub <- matrix(scan("NewMaille.txt", 0), ncol=1, byrow=TRUE)
subdata <- matrix(ncol = 101)
subLoc3D <- matrix(ncol = 3)
for (i in 1:length(sub)) {
	if (sub[i] != 0) {
		subdata <- rbind(subdata, data[i,])
		subLoc3D <- rbind(subLoc3D, loc3D[i,])
	}
}
subdata <- subdata[2:dim(subdata)[1],]
subLoc3D <- subLoc3D[2:dim(subLoc3D)[1],]
subLoc2D <- subLoc3D[,1:2]
data <- aperm(log(data), c(2, 1))
subdata <- aperm(log(subdata), c(2, 1))

# Chargement de G0
G0 <- matrix(scan("G0bis.txt", 0), ncol=713, byrow=TRUE)

#################################################################################
# appel de selectFast avec penalite l1 + G0
# G0 est une matrice p x p de 1 et 0, symmetric avec 0 sur la diagonale
#################################################################################

source("GGMselectG0.R")
source("ModLassoG0.R")

## Lasso LA

graphHippoLasso<-selectFast(subdata,family="LA",G0=G0,K=1.5,gamma0=0)

library(network)

# graphHippoK1LassoG <- network(graphHippoLasso$LA$G[,,1])
graphHippoK15LassoG <- network(graphHippoLasso$LA$G)
# graphHippoK2LassoG <- network(graphHippoLasso$LA$G[,,3])

# pdf(file="pics/LassowithG0K11gamma0_0.pdf")
# plot(graphHippoK1LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
# dev.off()
pdf(file="pics/LassowithG0bisK15gamma0_0.pdf")
plot(graphHippoK15LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
# pdf(file="pics/LassowithG0K2gamma0_0.pdf")
# plot(graphHippoK2LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
# dev.off()

## Lasso EW

graphHippoLasso<-selectFast(subdata,family="EW",G0=G0,K=1.5,gamma0=0)

#graphHippoK1EWG<- network(graphHippoLasso$EW$G)
#graphHippoK15EWG<- network(graphHippoLasso$EW$G)
#graphHippoK2EWG<- network(graphHippoLasso$EW$G[,,3])


#pdf(file="pics/EWwithG0K11gamma0_0.pdf")
#plot(graphHippoK1EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
##dev.off()
pdf(file="pics/EWwithG0bisK15gamma0_0.pdf")
plot(graphHippoK15EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
#pdf(file="pics/EWwithG0K2gamma0_0.pdf")
#plot(graphHippoK2EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()


## ridge
## Lasso LA


graphHippoLasso<-selectFast(subdata,family="LA",G0=G0,K=1.1,gamma0=0.01)

#graphHippoK1LassoG <- network(graphHippoLasso$LA$G[,,1])
graphHippoK15LassoG <- network(graphHippoLasso$LA$G)
#graphHippoK2LassoG <- network(graphHippoLasso$LA$G[,,3])

#pdf(file="pics/LassowithG0K11gamma0_0-005.pdf")
#plot(graphHippoK1LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
pdf(file="pics/LassowithG0bisK15gamma0_0-005.pdf")
plot(graphHippoK15LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
#pdf(file="pics/LassowithG0K2gamma0_0-005.pdf")
#plot(graphHippoK2LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

## Lasso EW

graphHippoLasso<-selectFast(subdata,family="EW",G0=G0,K=1.5,gamma0=0.005)

#graphHippoK1EWG<- network(graphHippoLasso$EW$G[,,1])
graphHippoK15EWG<- network(graphHippoLasso$EW$G)
#graphHippoK2EWG<- network(graphHippoLasso$EW$G[,,3])


#pdf(file="pics/EWwithG0K11gamma0_0-005.pdf")
#plot(graphHippoK1EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
pdf(file="pics/EWwithG0bisK15gamma0_0-005.pdf")
plot(graphHippoK15EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
#pdf(file="pics/EWwithG0K2gamma0_0-005.pdf")
#plot(graphHippoK2EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()


## Sans G0
## Lasso LA

graphHippoLasso<-selectFast(subdata,family="LA",K=1.5,gamma0=0)

#graphHippoK1LassoG <- network(graphHippoLasso$LA$G[,,1])
graphHippoK15LassoG <- network(graphHippoLasso$LA$G)
#graphHippoK2LassoG <- network(graphHippoLasso$LA$G[,,3])

#pdf(file="pics/LassowithoutG0K11.pdf")
#plot(graphHippoK1LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
pdf(file="pics/LassowithoutG0bisK15.pdf")
plot(graphHippoK15LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
#pdf(file="pics/LassowithoutG0K2.pdf")
#plot(graphHippoK2LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

## Lasso EW

graphHippoLasso<-selectFast(subdata,family="EW",K=1.5,gamma0=0)

#graphHippoK1EWG<- network(graphHippoLasso$EW$G[,,1])
graphHippoK15EWG<- network(graphHippoLasso$EW$G)
#graphHippoK2EWG<- network(graphHippoLasso$EW$G[,,3])


#pdf(file="pics/EWwithoutG0K11.pdf")
#plot(graphHippoK1EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
pdf(file="pics/EWwithoutG0bisK15.pdf")
plot(graphHippoK15EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
#pdf(file="pics/EWwithoutG0K2.pdf")
#plot(graphHippoK2EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()


###############################################################################
# appel de selectFast avec penalite l1 + l2 + G0
# G0 est une matrice p x p de 1 et 0, symmetric avec 0 sur la diagonale
# gamma >0 (scalaire)
##############################################################################

source("GGMselectEnetG0.R")
source("ModLassoEnetG0.R")

gamma=2.7

## LA

graphHippoEnet<-selectFast(subdata,family="LA",G0=G0,K=1.5,gamma=gamma,gamma0=0)

#graphHippoEnetGK11 <- network(graphHippoEnet$LA$G[,,1])
graphHippoEnetGK15 <- network(graphHippoEnet$LA$G)
#graphHippoEnetGK2 <- network(graphHippoEnet$LA$G[,,3])

#pdf(file="pics/EnetLAwithG0K11gamma1gamma0_0.pdf")
#plot(graphHippoEnetGK11, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
pdf(file="pics/EnetLAwithG0bisK15gamma27gamma0_0.pdf")
plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
#pdf(file="pics/EnetLAwithG0K2gamma1gamma0_0.pdf")
#plot(graphHippoEnetGK2, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

## EW

graphHippoEnet<-selectFast(subdata,family="EW",G0=G0,K=1.5,gamma=gamma,gamma0=0)

#graphHippoEnetGK11 <- network(graphHippoEnet$EW$G[,,1])
graphHippoEnetGK15 <- network(graphHippoEnet$EW$G)
#graphHippoEnetGK2 <- network(graphHippoEnet$EW$G[,,3])

#pdf(file="pics/EnetEWwithG0K11gamma1gamma0_0.pdf")
#plot(graphHippoEnetGK11, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
pdf(file="pics/EnetEWwithG0bisK15gamma27gamma0_0.pdf")
plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
#pdf(file="pics/EnetEWwithG0K2gamma1gamma0_0.pdf")
#plot(graphHippoEnetGK2, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

## gamma0 non nul
## LA

graphHippoEnet<-selectFast(subdata,family="LA",G0=G0,K=1.5,gamma=gamma,gamma0=0.05)

#graphHippoEnetGK11 <- network(graphHippoEnet$LA$G[,,1])
graphHippoEnetGK15 <- network(graphHippoEnet$LA$G)
#graphHippoEnetGK2 <- network(graphHippoEnet$LA$G[,,3])

#pdf(file="pics/EnetLAwithG0K11gamma1gamma0_0-005.pdf")
#plot(graphHippoEnetGK11, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
pdf(file="pics/EnetLAwithG0bisK15gamma1gamma0_0-005.pdf")
plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
#pdf(file="pics/EnetLAwithG0K2gamma1gamma0_0-005.pdf")
#plot(graphHippoEnetGK2, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

## EW

graphHippoEnet<-selectFast(subdata,family="EW",G0=G0,K=1.5,gamma=gamma,gamma0=0.005)

#graphHippoEnetGK11 <- network(graphHippoEnet$EW$G[,,1])
graphHippoEnetGK15 <- network(graphHippoEnet$EW$G)
#graphHippoEnetGK2 <- network(graphHippoEnet$EW$G[,,3])

#pdf(file="pics/EnetEWwithG0K11gamma1gamma0_0-005.pdf")
#plot(graphHippoEnetGK11, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
pdf(file="pics/EnetEWwithG0bisK15gamma27gamma0_0-005.pdf")
plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
#pdf(file="pics/EnetEWwithG0K2gamma1gamma0_0-005.pdf")
#plot(graphHippoEnetGK2, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

## Sans G0
## LA

graphHippoEnet<-selectFast(subdata,family="LA",K=1.5,gamma=gamma,gamma0=0)

#graphHippoEnetGK11 <- network(graphHippoEnet$LA$G[,,1])
graphHippoEnetGK15 <- network(graphHippoEnet$LA$G)
#graphHippoEnetGK2 <- network(graphHippoEnet$LA$G[,,3])

#pdf(file="pics/EnetLAwithoutG0K11gamma1.pdf")
#plot(graphHippoEnetGK11, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
pdf(file="pics/EnetLAwithoutG0bisK15gamma1.pdf")
plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
#pdf(file="pics/EnetLAwithG0outK2gamma1.pdf")
#plot(graphHippoEnetGK2, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

## EW

graphHippoEnet<-selectFast(subdata,family="EW",K=1.5,gamma=gamma,gamma0=0)

#graphHippoEnetGK11 <- network(graphHippoEnet$EW$G[,,1])
graphHippoEnetGK15 <- network(graphHippoEnet$EW$G)
#graphHippoEnetGK2 <- network(graphHippoEnet$EW$G[,,3])

#pdf(file="pics/EnetEWwithoutG0K11gamma1.pdf")
#plot(graphHippoEnetGK11, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
pdf(file="pics/EnetEWwithoutG0bisK15gamma1.pdf")
plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
#pdf(file="pics/EnetEWwithoutG0K2gamma1.pdf")
#plot(graphHippoEnetGK2, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
