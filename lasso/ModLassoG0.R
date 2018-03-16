

calcModLasso <- function(X, Dmax, method, max.st,
                          max.iter, eps,
                         beta,
                         tau,  h,  T0, G0, d0,gamma0) {
  # ---------------------------------------------------------------
  # FUNCTION
  #   Select a set of candidate edges, ranked by importance
  # INPUT
  #   X: n x p matrix
  #   Dmax: p dim. vector
  #   method: vector taking values in "EW", "LA", "C01"
  #   max.st: scalar ( min(n,p-1))
  #   max.iter : positive integer number (maximal number of
  #                iterations)
  #   eps  : positive real number (precision)
  #   G0 : p x p adjacency matrix
  #   d0 : p-dimensional vector of integers
  # OUTPUT
  #   GrGlob: matrix with 3 columns and a variable number of rows
  #   GrGlob[,1] = nodes
  #   GrGlob[,2] = neighbour of the node
  #   GrGlob[,3] = regularisation parameter
  #      (ranked in decreasing order)
  # CALLED BY
  #   calcLarsNEW when family=LA or EW
  # ---------------------------------------------------------------
  
  GrGlob  <- NULL
  Gr0 <- NULL
  p <- dim(X)[2]
  n <- dim(X)[1]
  Dmax <- pmax(Dmax-d0,0) # G0
  Nm0 <- 0
  if (method == "LA")
    betainit <- rep(1,p-1)
  
  for (a in 1:p) {
    modAgarder <- max.st
    Y <- X[,a]
    Z <- X[,-a]
    
    # indices m0 and projection
	m0 <- (1:p)[G0[a,]==1]  #!warning! m0 takes value in 1:p and not in 1:(p-1)    # G0 
	if (length(m0)>0) {
		Y <- Y-X[,m0]%*%solve(t(X[,m0])%*%X[,m0]+diag(gamma0,nrow=length(m0),ncol=length(m0)),t(X[,m0])%*%Y) # G0
		Nm0 <- Nm0 +1
		}
    if (method == "EW")
             betainit <- rep(1,p-1) # sapply(EW(x=Z,y=Y,beta,tau,h,T0, max.iter, eps),c)
    
    U <-  t(t(Z)*abs(betainit))
    ll <- lars(U, Y, normalize=FALSE, intercept=FALSE, max.steps=modAgarder)
    action <- unlist(ll$action)
    sign.act <- sign(action)
    Vois <- ((1:p)[-a])[abs(action)]*sign.act
    switch(method,
           EW =
           {
             modAgarder <-min(modAgarder,length((1:length(action))[cumsum(sign.act) <= Dmax[a]]))
           },
           LA =
           {
             modAgarder <- min(modAgarder,length(action))
           }) #end switch
    
    result <- cbind(rep(a,modAgarder),
                    Vois[1:modAgarder],
                    ll$lambda[1:modAgarder])
    GrGlob <- rbind(GrGlob,result)
    if (length(m0)>0) Gr0 <- rbind(Gr0,cbind(rep(a,d0[a]),m0)) # G0
  } # end a

  dimnames(GrGlob)[[2]] <- c("a","Vois(a)","lambda")
  # the rows of GrGlob are ordered according to lambda
  ind <- order(GrGlob[,3],decreasing=TRUE)
  GrGlob <- GrGlob[ind,]
  
  # Add Gr0
  if (Nm0>0) {
  	Gr0 <- cbind(Gr0,rep(2*GrGlob[1,3],dim(Gr0)[1])) # G0
  	GrGlob <- rbind(Gr0,GrGlob) # G0
  }
  return(GrGlob)
}

