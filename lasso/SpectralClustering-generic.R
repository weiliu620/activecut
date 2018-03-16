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

library(network)

#################################################################################
# appel de selectFast avec penalite l1 + G0
# G0 est une matrice p x p de 1 et 0, symmetric avec 0 sur la diagonale
#################################################################################

source("GGMselectG0.R")
source("ModLassoG0.R")

graphHippo<-selectFast(subdata,family="EW",G0=G0,K=1.5,gamma0=0)


##########
# load info #
##########
library(scatterplot3d)
library(mclust)
#load("/Users/K/Desktop/hippo/graphHippoLARSK3.Rdata")
#DEFINIR Q!
Q <- graphHippo$EW$G
p<- dim(Q)[2]
#Z<-as.matrix(read.table("Loc3D.txt"))  # position 3D
Z <- subLoc3D
ZZ <- Z%*%t(Z)
D3<- diag(ZZ)%*%t(pmin(1,1:p))-ZZ
D3<-D3+t(D3)  # distance i-j square

############
# generate Q #
############
Q <- graphHippo$EW$G
Qseuil <- Q
Qseuil[D3<1000]<-0
#Q <- Qseuil+exp(-D3)
#Q <- D3*Q
#diag(Q) <- -Q%*% rep(1,p)

ecarttype <- 1.4128
G <- exp(-D3/(2* ecarttype^2))/(ecarttype * sqrt(2*pi))
	nbVois <- apply(Qseuil, 1, sum)
for(i in 1:p) {
	sumNormal <- 0
	for(j in i:p) {
		if(Qseuil[j,i] == 1) {
			G[j,i] == 1/(nbVois[[i]] + 2)
		}
		else {
			sumNormal <- sumNormal + G[j,i]
		}
	}
	for(j in 1:p) {
		if(G[j,i] != 1/(nbVois[[i]] + 2)) {
			#G[j, i] = G[j, i] / (4 * sumNormal)
			G[j, i] = G[j, i] / ((nbVois[[i]] + 2) * sumNormal)
		}
	}
	G[i,i] <- 1/(nbVois[[i]] + 2)
}
#diag(G) <- -G %*%rep(1,p)

# source("../neighbor.R")
# Qsans <- substractG0(Q, subLoc3D, 3)
# Qavec <- substractG0(Q, subLoc3D, 3) + buildG0(subLoc3D, 3)
# G0 <- buildG0(subLoc3D,3)

#Q_christophe <- G0 + 50 * Qsans
#diag(Q_christophe) <- -Q_christophe%*%rep(1,p)


#Q <- Q_christophe
Q <- G


#################
# spectral analysis #
#################
eigenQ <- eigen(Q,symmetric=TRUE)
pdf(file="pics/spectral.pdf")
par(mfrow=c(1,2))
plot(1:p,exp(5*eigenQ$values),main="eigenvalues")
data <- eigenQ$vectors[,2:4]
scatterplot3d(data[,1],data[,2],data[,3], main="eigenvectors")


#######################
# clustering avec kmeans #
#######################
datakmeans<-kmeans(data, centers=3, iter.max=2000,nstart=2000)
pdf(file="pics/spectralLassoEWwithG0.pdf")
classifkmeans <- datakmeans$cluster
par(mfrow=c(1,3))
scatterplot3d(data[,3],data[,1],data[,2],color=classifkmeans, main="eigenvectors")
scatterplot3d(Z[,3],Z[,1],Z[,2],color=classifkmeans,main="clustering")
scatterplot3d(Z[,3],-Z[,1],Z[,2],color=classifkmeans,main="clustering")
dev.off()


#
############################################################################
#                                           POUR MEMOIRE
############################################################################
#
#######################
## clustering avec Mclust #
#######################
#dataclust<-Mclust(data,G=3:3,modelNames="VVV") # fit un melange de 3 gaussiennes
#classif <- dataclust$classification
#par(mfrow=c(1,2))
#scatterplot3d(data[,1],data[,2],data[,3],color=classif, main="eigenvectors")
#scatterplot3d(Z[,3],Z[,1],Z[,2],color=classif,main="clustering")
#
#
#
#########
## divers #
#########
#diag(Q) <- -Q%*%pmin(1,1:p)
#eigenQ <- eigen(Q,symmetric=TRUE)
#par(mfrow=c(1,2)
#plot(1:p,exp(10*eigenQ$values))
#plot(1:p,1/(1-10*eigenQ$values))
#dev.off()
#plot(eigenQ$vectors[,2],eigenQ$vectors[,3])
#library(scatterplot3d)
#scatterplot3d(eigenQ$vectors[,2],eigenQ$vectors[,4],eigenQ$vectors[,3])
#
#Z<-as.matrix(read.table("/Users/K/Desktop/hippo/Loc3D.txt"))
#ZZ <- Z%*%t(Z)
#D3<- diag(ZZ)%*%t(pmin(1,1:p))-ZZ
#D3<-D3+t(D3)
#
#Q <- graphHippoLARSK2$LA$G
#Q <- D3*Q
#diag(Q) <- -Q%*%pmin(1,1:p)
#eigenQ <- eigen(Q,symmetric=TRUE)
#plot(1:p,exp(10*eigenQ$values))
#plot(eigenQ$vectors[,2],eigenQ$vectors[,3])
#library(scatterplot3d)
#scatterplot3d(eigenQ$vectors[,2],eigenQ$vectors[,3],eigenQ$vectors[,4])
#data <- eigenQ$vectors[,2:4]
#scatterplot3d(data[,1],data[,2],data[,3])
#dataclust<-Mclust(data,G=3:3,modelNames="VVV")
#classif <- dataclust$classification
#scatterplot3d(Z[,3],Z[,1],Z[,2],color=classif)
#
