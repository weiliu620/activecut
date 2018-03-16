library(GGMselect)
dyn.load("loopG0.so")
source("penaltyG0.R")
source("convertGraph.R")

# Chargement des donnees
dataX <- matrix(scan("X410.txt", 0), ncol=410, byrow=TRUE)
for (iK in 1:409){
	dataX[,iK]= dataX[,iK]-mean(dataX[,iK])
	dataX[,iK]<-dataX[,iK]/sqrt(var(dataX[,iK]))
	}

dataY <- matrix(scan("Y410.txt", 0), ncol=12, byrow=TRUE)

for (iK in 1:12){
	dataY[,iK]= dataY[,iK]-mean(dataY[,iK])
	dataY[,iK]<-dataY[,iK]/sqrt(var(dataY[,iK]))
	}

data <- array(0,c(410,421))
data[,1:409]<-dataX[,1:409]
data[,410:421]<-dataY

# Chargement de G0
#G0 <- matrix(scan("G0.txt", 0), ncol=713, byrow=TRUE)

#################################################################################
# appel de selectFast avec penalite l1 + G0
# G0 est une matrice p x p de 1 et 0, symmetric avec 0 sur la diagonale
#################################################################################

source("GGMselectG0.R")
source("ModLassoG0.R")

source("MyFamilyG0.R")

## Lasso LA

graphVitesseLasso<-selectFast(data,family="LA",K=c(1.1,1.5,2))

library(network)

graphVitesseK1LassoG <- network(graphVitesseLasso$LA$G[,,1])
graphVitesseK15LassoG <- network(graphVitesseLasso$LA$G)
graphVitesseK2LassoG <- network(graphVitesseLasso$LA$G[,,3])

pdf(file="pics/VitesseLassoK11.pdf")
plot(graphVitesseK1LassoG, usearrows = FALSE,vertex.cex=0.5)
dev.off()
#pdf(file="pics/LassowithG01K15gamma0_0.pdf")
#plot(graphHippoK15LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
# pdf(file="pics/LassowithG0K2gamma0_0.pdf")
# plot(graphHippoK2LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
# dev.off()

## Lasso EW
acrit<-NULL

graphHippoLasso<-selectFast(subdata,family="EW",G0=G0,K=1.1,gamma0=0)

acrit[1]<-graphHippoLasso$EW$crit.min

#graphHippoK1EWG<- network(graphHippoLasso$EW$G)
#graphHippoK15EWG<- network(graphHippoLasso$EW$G)
#graphHippoK2EWG<- network(graphHippoLasso$EW$G[,,3])


#pdf(file="pics/EWwithG0K11gamma0_0.pdf")
#plot(graphHippoK1EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
##dev.off()
#pdf(file="pics/EWwithG01K15gamma0_0.pdf")
#plot(graphHippoK15EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EWwithG0K2gamma0_0.pdf")
#plot(graphHippoK2EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()


## ridge
## Lasso LA

graphHippoLasso<-selectFast(subdata,family="LA",G0=G0,K=1.1,gamma0=0.005)
acrit[2]<-graphHippoLasso$LA$crit.min

#graphHippoK1LassoG <- network(graphHippoLasso$LA$G[,,1])
#graphHippoK15LassoG <- network(graphHippoLasso$LA$G)
#graphHippoK2LassoG <- network(graphHippoLasso$LA$G[,,3])

#pdf(file="pics/LassowithG0K11gamma0_0-005.pdf")
#plot(graphHippoK1LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/LassowithG01K15gamma0_0-01.pdf")
#plot(graphHippoK15LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/LassowithG0K2gamma0_0-005.pdf")
#plot(graphHippoK2LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

## Lasso EW

graphHippoLasso<-selectFast(subdata,family="EW",G0=G0,K=1.1,gamma0=0.005)
acrit[3]<-graphHippoLasso$EW$crit.min

#graphHippoK1EWG<- network(graphHippoLasso$EW$G[,,1])
#graphHippoK15EWG<- network(graphHippoLasso$EW$G)
#graphHippoK2EWG<- network(graphHippoLasso$EW$G[,,3])


#pdf(file="pics/EWwithG0K11gamma0_0-005.pdf")
#plot(graphHippoK1EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EWwithG01K15gamma0_0-005.pdf")
#plot(graphHippoK15EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EWwithG0K2gamma0_0-005.pdf")
#plot(graphHippoK2EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()


## Sans G0
## Lasso LA

#graphHippoLasso<-selectFast(subdata,family="LA",K=1.1,gamma0=0)

#graphHippoK1LassoG <- network(graphHippoLasso$LA$G[,,1])
#graphHippoK15LassoG <- network(graphHippoLasso$LA$G)
#graphHippoK2LassoG <- network(graphHippoLasso$LA$G[,,3])

#pdf(file="pics/LassowithoutG0K11.pdf")
#plot(graphHippoK1LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/LassowithoutG01K15.pdf")
#plot(graphHippoK15LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/LassowithoutG0K2.pdf")
#plot(graphHippoK2LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

## Lasso EW

#graphHippoLasso<-selectFast(subdata,family="EW",K=1.1,gamma0=0)


#graphHippoK1EWG<- network(graphHippoLasso$EW$G[,,1])
#graphHippoK15EWG<- network(graphHippoLasso$EW$G)
#graphHippoK2EWG<- network(graphHippoLasso$EW$G[,,3])


#pdf(file="pics/EWwithoutG0K11.pdf")
#plot(graphHippoK1EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EWwithoutG01K15.pdf")
#plot(graphHippoK15EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
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

gamma=1

## LA

graphHippoEnet<-selectFast(subdata,family="LA",G0=G0,K=1.1,gamma=gamma,gamma0=0)

acrit[4]<-graphHippoEnet$LA$crit.min
#graphHippoEnetGK11 <- network(graphHippoEnet$LA$G[,,1])
#graphHippoEnetGK15 <- network(graphHippoEnet$LA$G[,,2])
#graphHippoEnetGK2 <- network(graphHippoEnet$LA$G[,,3])

#pdf(file="pics/EnetLAwithG0K11gamma1gamma0_0.pdf")
#plot(graphHippoEnetGK11, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EnetLAwithG01K15gamma1gamma0_0.pdf")
#plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EnetLAwithG0K2gamma1gamma0_0.pdf")
#plot(graphHippoEnetGK2, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

## EW

graphHippoEnet<-selectFast(subdata,family="EW",G0=G0,K=1.1,gamma=gamma,gamma0=0)
acrit[5]<-graphHippoEnet$EW$crit.min

#graphHippoEnetGK11 <- network(graphHippoEnet$EW$G[,,1])
#graphHippoEnetGK15 <- network(graphHippoEnet$EW$G)
#graphHippoEnetGK2 <- network(graphHippoEnet$EW$G[,,3])

#pdf(file="pics/EnetEWwithG0K11gamma1gamma0_0.pdf")
#plot(graphHippoEnetGK11, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EnetEWwithG01K15gamma1gamma0_0.pdf")
#plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EnetEWwithG0K2gamma1gamma0_0.pdf")
#plot(graphHippoEnetGK2, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

## gamma0 non nul
## LA

graphHippoEnet<-selectFast(subdata,family="LA",G0=G0,K=1.1,gamma=gamma,gamma0=0.05)

acrit[6]<-graphHippoEnet$LA$crit.min

#graphHippoEnetGK11 <- network(graphHippoEnet$LA$G[,,1])
#graphHippoEnetGK15 <- network(graphHippoEnet$LA$G)
#graphHippoEnetGK2 <- network(graphHippoEnet$LA$G[,,3])

#pdf(file="pics/EnetLAwithG0K11gamma1gamma0_0-005.pdf")
#plot(graphHippoEnetGK11, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EnetLAwithG01K15gamma1gamma0_0-005.pdf")
#plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EnetLAwithG0K2gamma1gamma0_0-005.pdf")
#plot(graphHippoEnetGK2, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

## EW

graphHippoEnet<-selectFast(subdata,family="EW",G0=G0,K=1.1,gamma=gamma,gamma0=0.005)

acrit[7]<-graphHippoEnet$EW$crit.min

#graphHippoEnetGK11 <- network(graphHippoEnet$EW$G[,,1])
#graphHippoEnetGK15 <- network(graphHippoEnet$EW$G)
#graphHippoEnetGK2 <- network(graphHippoEnet$EW$G[,,3])

#pdf(file="pics/EnetEWwithG0K11gamma1gamma0_0-005.pdf")
#plot(graphHippoEnetGK11, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EnetEWwithG01K15gamma1gamma0_0-005.pdf")
#plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EnetEWwithG0K2gamma1gamma0_0-005.pdf")
#plot(graphHippoEnetGK2, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

## Sans G0
## LA

#graphHippoEnet<-selectFast(subdata,family="LA",K=1.5,gamma=gamma,gamma0=0)

#graphHippoEnetGK11 <- network(graphHippoEnet$LA$G[,,1])
#graphHippoEnetGK15 <- network(graphHippoEnet$LA$G)
#graphHippoEnetGK2 <- network(graphHippoEnet$LA$G[,,3])

#pdf(file="pics/EnetLAwithoutG0K11gamma1.pdf")
#plot(graphHippoEnetGK11, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EnetLAwithoutG01K15gamma1.pdf")
#plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EnetLAwithG0outK2gamma1.pdf")
#plot(graphHippoEnetGK2, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

## EW

#graphHippoEnet<-selectFast(subdata,family="EW",K=1.5,gamma=gamma,gamma0=0)

#graphHippoEnetGK11 <- network(graphHippoEnet$EW$G[,,1])
#graphHippoEnetGK15 <- network(graphHippoEnet$EW$G)
#graphHippoEnetGK2 <- network(graphHippoEnet$EW$G[,,3])

#pdf(file="pics/EnetEWwithoutG0K11gamma1.pdf")
#plot(graphHippoEnetGK11, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EnetEWwithoutG01K15gamma1.pdf")
#plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()
#pdf(file="pics/EnetEWwithoutG0K2gamma1.pdf")
#plot(graphHippoEnetGK2, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()

