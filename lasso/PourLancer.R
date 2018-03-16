# creer un repertoire GGMselectGO et y mettre tous les fichiers
# dans ce repertoire faire:
# R CMD SHLIB -o loopG0.so loopG0.c scr.c
# ensuite lancer R
# puis lancer avec
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


source("GGMselectG0.R")
source("ModLassoG0.R")
graphHippoLasso<-selectFast(subdata,family="EW",G0=G0,K=c(1,1.5,2,2.5),gamma0=0)
# tester aussi family="EW" ; K>1
# 
library(network)
# graphHippoK1LassoG <- network(graphHippoLasso$LA$G[,,1])
# graphHippoK15LassoG <- network(graphHippoLasso$LA$G[,,2])
# graphHippoK2LassoG <- network(graphHippoLasso$LA$G[,,3])
graphHippoK1EWG<- network(graphHippoLasso$EW$G[,,1])
graphHippoK15EWG<- network(graphHippoLasso$EW$G[,,2])
graphHippoK2EWG<- network(graphHippoLasso$EW$G[,,3])
graphHippoK25EWG<- network(graphHippoLasso$EW$G[,,4])

# 
# pdf(file="LassowithG0K1gamma0_0.pdf")
# plot(graphHippoK1LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
# dev.off()
# pdf(file="LassowithG0K15gamma0_0.pdf")
# plot(graphHippoK15LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
# dev.off()
# pdf(file="LassowithG0K2gamma0_0.pdf")
# plot(graphHippoK2LassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
# dev.off()

pdf(file="EWwithG0K1gamma0_0.pdf")
plot(graphHippoK1EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
pdf(file="EWwithG0K15gamma0_0.pdf")
plot(graphHippoK15EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
pdf(file="EWwithG0K2gamma0_0.pdf")
plot(graphHippoK2EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()
pdf(file="EWwithG0K25gamma0_0.pdf")
plot(graphHippoK25EWG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()

###############################################################################
# appel de selectFast avec penalite l1 + l2 + G0
# G0 est une matrice p x p de 1 et 0, symmetric avec 0 sur la diagonale
# gamma >0 (scalaire)
##############################################################################

source("GGMselectEnetG0.R")
source("ModLassoEnetG0.R")

gamma=1
graphHippoEnet<-selectFast(subdata,family="LA",G0=G0,K=1.1,gamma=gamma,gamma0=0.001)
graphHippoEnetG <- network(graphHippoEnet$LA$G)

pdf(file="EnetK1-1withG0gamma1gamma0_0-001.pdf")
plot(graphHippoEnetG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()

stop()

########################################################
# stabilité par rapport au sous échantillonnage
#######################################################

sub <- matrix(scan("NewMaille_1.txt", 0), ncol=1, byrow=TRUE);

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
G0 <- matrix(scan("G0_1.txt", 0), ncol=ncol(subdata), byrow=TRUE)

# appel de selectFast avec penalite l1 + G0
# G0 est une matrice p x p de 1 et 0, symmetric avec 0 sur la diagonale

source("GGMselectG0.R")
source("ModLassoG0.R")
graphHippoLasso<-selectFast(subdata,family="LA",G0=G0,gamma0=1)
# 
library(network)
graphHippoLassoG <- network(graphHippoLasso$LA$G)
# 
pdf(file="LassowithG0_1gamma0_1.pdf")
plot(graphHippoLassoG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()


# appel de selectFast avec penalite l1 + l2 + G0
# G0 est une matrice p x p de 1 et 0, symmetric avec 0 sur la diagonale
# gamma >0 (scalaire)

source("GGMselectEnetG0.R")
source("ModLassoEnetG0.R")

gamma=1
graphHippoEnet<-selectFast(subdata,family="LA",G0=G0,gamma=gamma,gamma0=0)
graphHippoEnetG <- network(graphHippoEnet$LA$G)

pdf(file="EnetwithG0_1gamma1gamma0_0.pdf")
plot(graphHippoEnetG, usearrows = FALSE, coord = subLoc2D,vertex.cex=0.5)
dev.off()



