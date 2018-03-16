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
G0 <- matrix(scan("G0.txt", 0), ncol=713, byrow=TRUE)

#################################################################################
# appel de selectFast avec penalite l1 + G0
# G0 est une matrice p x p de 1 et 0, symmetric avec 0 sur la diagonale
#################################################################################

library(network)

source("GGMselectG0.R")
source("ModLassoG0.R")

## Lasso EW

graphHippoLasso<-selectFast(subdata[1:44,],family="EW",G0=G0,K=1.1,gamma0=0)

graphHippoK15EWG<- network(graphHippoLasso$EW$G)

pdf(file="pics/Pop1_44_EWwithG0K11gamma0_0.pdf")
plot(graphHippoK15EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()


graphHippoLasso<-selectFast(subdata[58:101,],family="EW",G0=G0,K=1.5,gamma0=0)

graphHippoK15EWG<- network(graphHippoLasso$EW$G)

pdf(file="pics/Pop2_EWwithG0K15gamma0_0.pdf")
plot(graphHippoK15EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()

source("GGMselectG0.R")
source("ModLassoG0.R")

## Lasso EW

graphHippoLasso<-selectFast(subdata[1:57,],family="EW",G0=G0,K=1.1,gamma0=0)

graphHippoK15EWG<- network(graphHippoLasso$EW$G)

pdf(file="pics/Pop1_EWwithG0K11gamma0_0.pdf")
plot(graphHippoK15EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()


graphHippoLasso<-selectFast(subdata[58:101,],family="EW",G0=G0,K=1.1,gamma0=0)

graphHippoK15EWG<- network(graphHippoLasso$EW$G)

pdf(file="pics/Pop2_EWwithG0K11gamma0_0.pdf")
plot(graphHippoK15EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()


## ridge
## Lasso EW

graphHippoLasso<-selectFast(subdata[1:57,],family="EW",G0=G0,K=1.1,gamma0=0.005)

graphHippoK15EWG<- network(graphHippoLasso$EW$G)

pdf(file="pics/Pop1_EWwithG0K11gamma0_0-005.pdf")
plot(graphHippoK15EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()


graphHippoLasso<-selectFast(subdata[58:101,],family="EW",G0=G0,K=1.1,gamma0=0.005)

graphHippoK15EWG<- network(graphHippoLasso$EW$G)

pdf(file="pics/Pop2_EWwithG0K11gamma0_0-005.pdf")
plot(graphHippoK15EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()

## Sans G0
## Lasso LA

graphHippoLasso<-selectFast(subdata[1:57,],family="LA",K=1.1,gamma0=0)

graphHippoK15LassoG <- network(graphHippoLasso$LA$G)

pdf(file="pics/Pop1_LassowithoutG0K11.pdf")
plot(graphHippoK15LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()


graphHippoLasso<-selectFast(subdata[58:101,],family="LA",K=1.1,gamma0=0)

graphHippoK15LassoG <- network(graphHippoLasso$LA$G)

pdf(file="pics/Pop2_LassowithoutG0K11.pdf")
plot(graphHippoK15LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()


###############################################################################
# appel de selectFast avec penalite l1 + l2 + G0
# G0 est une matrice p x p de 1 et 0, symmetric avec 0 sur la diagonale
# gamma >0 (scalaire)
##############################################################################

source("GGMselectEnetG0.R")
source("ModLassoEnetG0.R")

gamma=1

## LA 

graphHippoEnet<-selectFast(subdata[12:55,],family="LA",G0=G0,K=1.1,gamma=gamma,gamma0=0)

graphHippoEnetGK15 <- network(graphHippoEnet$LA$G)

pdf(file="pics/Pop1_44_EnetLAwithG0K11gamma1gamma0_0.pdf")
plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()


# R crash complètement sans raison apparante

#graphHippoEnet<-selectFast(subdata[58:101,],family="LA",G0=G0,K=1.5,gamma=gamma,gamma0=0)

#graphHippoEnetGK15 <- network(graphHippoEnet$LA$G)

#pdf(file="pics/Pop2_EnetLAwithG0K15gamma1gamma0_0.pdf")
#plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
#dev.off()



## gamma0 non nul
## LA

graphHippoEnet<-selectFast(subdata[1:44,],family="LA",G0=G0,K=1.1,gamma=gamma,gamma0=0.005)

graphHippoEnetGK15 <- network(graphHippoEnet$LA$G)

pdf(file="pics/Pop1_44_EnetLAwithG0K11gamma1gamma0_0-005.pdf")
plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()


graphHippoEnet<-selectFast(subdata[58:101,],family="LA",G0=G0,K=1.1,gamma=gamma,gamma0=0.005)

graphHippoEnetGK15 <- network(graphHippoEnet$LA$G)

pdf(file="pics/Pop2_EnetLAwithG0K11gamma1gamma0_0-005.pdf")
plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()

######################################
#############STOP####################
######################################


## Sans G0
## LA

graphHippoEnet<-selectFast(subdata[1:57,],family="LA",K=1.1,gamma=gamma,gamma0=0)

graphHippoEnetGK15 <- network(graphHippoEnet$LA$G)

pdf(file="pics/Pop1_EnetLAwithoutG0K11gamma1.pdf")
plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()

graphHippoEnet<-selectFast(subdata[58:101,],family="LA",K=1.1,gamma=gamma,gamma0=0)

graphHippoEnetGK15 <- network(graphHippoEnet$LA$G)

pdf(file="pics/Pop2_EnetLAwithoutG0K11gamma1.pdf")
plot(graphHippoEnetGK15, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
