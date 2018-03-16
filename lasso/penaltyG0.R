
penaltyG0 <- function(p,n, d0, dmax=min(3,n-3,p-1), K=2.5) {
  # ---------------------------------------------------------------
  # FUNCTION
  #   Compute the penalty  
  #   End-user function.
  # INPUT
  #   p,n: integers. Number of columns and number of rows of the data
  #   K : scalar or vector  (tuning parameter)
  #   dmax : scalar or p dimensional vector
  #         (maximum degree of the nodes of the graph)
  #         positive integers <=  n-3, p-1
  #         Default value: min(c(3,n-3,p-1))
  #	  d0 : p dimensional vector of integers = degrees of the nodes of G0 #G0
  # OUTPUT
  #  pen :  array (max(Dmax)+1, lK, p), with Dmax=max(dmax), lK=length(K) #G0
  # CALLED BY
  #  selectFast selectQE selectMyFam
  #   End-user function.
  # ---------------------------------------------------------------

  # Input checking
    if (any(K < 0))
    stop("K must be greater or equal to 0")
  # (K=0 works but is not advised)
  if (p<2)
    stop("p ", p, ", must be greater than 1")
  if (n<4)
    stop("n ",n ,", must be greater than 3")
  # d0 check #G0
  if (length(d0) != p)
	stop("d0 must be of size ",p) #G0
  if ( any(round(d0) != d0)) 
      stop("d0 must be a vector of integers") #G0
  if (any(d0>dmax))
	  stop("d0 must be not larger than dmax")  #G0
  # dmax verification
  ldmax <-length(dmax) 
  if ( any(dmax<1) ||  any(dmax > (n-3)) || any(dmax > (p-1)))
    stop("dmax must be greater than 0 and less than n - 3 and p - 1")
  if ( any(round(dmax) != dmax))
      stop("dmax must be integer")
  if ( (ldmax != 1) && (ldmax != p ))
      stop("dmax must be of length 1 or p")
  if (ldmax ==1) Dmax <- rep(dmax,p)
  else
    Dmax <- dmax
  # End of checking

  lK <- length(K)
  D <- max(Dmax)
  pen <- array(0,c(D+1,lK,p)) #G0
  Dm <- 0:D
  dimnames(pen) <- list(paste("D=",Dm),paste("K=",K),paste("p=",1:p)) #G0
  Nm <- n-Dm
  # dm <- pmin(Dm,rep(p/2,D+1))
  for (a in 1:p) { #G0
    Dm0 <- pmax(Dm-rep(d0[a],D+1),0) #G0  ie Dm0=(Dm-d0[a])+
	dm0 <- pmin(Dm0,rep((p-d0[a])/2,D+1)) #G0	ie dm0=min((Dm-d0[a])+,(p-d0[a])/2)
	Lm0 <- lgamma(p-d0[a])-lgamma(dm0+1)-lgamma(p-d0[a]-dm0)+2*log(Dm0+1) #G0	
	# en clair:  Lm0 = log(choose(p-1-d0[a],dm0))+2*log(Dm0+1)
	# old Lm <- lgamma(p)-lgamma(dm+1)-lgamma(p-dm)+2*log(Dm+1)
	EDkhi <- calcEDkhi(Dm+1,Nm-1,Lm0,p)
	if ( EDkhi$err != 0)
		{
		# error or warning cases
		mess <- paste(EDkhi$mess, "\n err=", EDkhi$err, "n=", n, "p=", p, "max(Dmax)=", max(Dmax))
		# On n'ecrit pas le warning quand EDkhi$err ==2:
		# c'est le programme C qui utilise les penalites
		# qui mettra un message dans le cas ou il
		# utilise au moins une des valeurs
		# de pen mises a Inf par calcEDkhi. Ceci afin que
		# le warning n'apparaisse pas si on n'utilise pas les valeurs
		# de pen concernees et n'apparaisse qu'une fois
		if (EDkhi$err !=2) {
			stop(mess)
		}
    } # fin err
	EDkhiR <- EDkhi$EDkhi2
	for (iK in 1:lK) {
		pen[,iK,a] <- K[iK]*(Nm/(Nm-1))*EDkhiR
	}
  } #G0	fin for a 
  return(pen)
  } # fin de penalty


calcEDkhi <- function(D,N,L,p) {
  # ---------------------------------------------------------------
  # FUNCTION
  #  For computation of the penalty
  # INPUT
  #   D, N vectors of integers
  #   K : vector
  #   p: integer
  # OUTPUT
  #   A list with components
  #     EDkhi2 value of EDkhi
  #     err : integer. Error or Warning code.
  #     mess : error message if err not equal 0
  # CALLED BY
  #   penalty
  # ---------------------------------------------------------------
     Bsup <- 1e+08
	 counter <- 0
    err <- 0
    mess <- NULL
    xq <- 0*D
    for (i in 2:length(D)) {
      if (xq[i-1] > Bsup) {
        xq[i:length(D)] <- Inf
        err <- 2
        mess <- paste("\nThe values", i,"to", length(D),
                      "of the penalty function greater than", Bsup, "are set to Inf")
        break
      }
	  if(L[i] != 0) { #G0
      counter <- counter+1 #G0
	  if (counter==1) { #G0
        xInf <- 0
        if (N[2]>5) {
          Delta <- (L[2]+log(5)+1/N[2])/(1-5/N[2])
          U <- sqrt((1+2*D[2]/(N[2]+2))*2*Delta/D[2])
          xSup <- D[2]*(1+exp(2*Delta/(N[2]+2))*U)**2*10
        }     else {
          xSup <- seq(10,100,by=10)*p
          fxSup <- EDkhi1(xSup,D[i],N[i],exp(-L[i]))
          if (max(fxSup) <0 ) {
            err <- 4
            mess <- paste("\nThe values of the penalty function cannot be calculated\n", xq)
          break
          } else {
            xSup <- min(xSup[fxSup>=0])
          }
        }
      }
      else {
#       (counter>1) #G0
        xInf <- xq[i-1]
        xSup <- xq[i-1]*seq(10,100,by=10)
        fxSup <- EDkhi1(xSup,D[i],N[i],exp(-L[i]))
        if (max(fxSup) <0 ) {
          err <- 3
          mess <- paste("\nThe values of the penalty function cannot be calculated\n", xq)
          break
        }
        xSup <- min(xSup[fxSup>=0])
      }
      if (L[i] < 50) {
        xx <-
          try(uniroot(EDkhi1,lower=xInf,upper=xSup,D=D[i],
                      N=N[i],q=exp(-L[i])))
        if (!is.list(xx)) {
          err <- 1
          mess <- xx
          break
        }
        xq[i] <- xx$root
      }
      else {
# L[i] >= 50
        fxSup <- EDkhi2(xSup,D[i],N[i],-L[i])
        if (fxSup <0) xSup<- 2*xSup
        xx <-
          try(uniroot(EDkhi2,lower=xInf,upper=xSup,D=D[i],
                      N=N[i],logq=-L[i]))
        if (!is.list(xx)) {
          err <- 1
          mess <- xx
          break
        }
        xq[i] <- xx$root
      }
	  } 
    }
    return(list(EDkhi2=xq,err=err, mess=mess))
  } # fin de EDkhi


EDkhi1 <- function(x,D,N,q) {
  # ---------------------------------------------------------------
  # FUNCTION
  #   Compute the penalty when L not too  small
  # INPUT
  #   x scalar or vector
  #   D, N : positive integer
  #   q : scalar
  # OUTPUT
  #   values in x : scalar or vector
  # CALLED BY
  #  EDkhi
  # ---------------------------------------------------------------
  fct1 <- pf(q=x/(D+2), df1=D+2, df2=N, lower.tail = FALSE)
  fct2 <- (x/D)*pf(q=(N+2)*x/(N*D), df1=D, df2=N+2, lower.tail = FALSE)
  return(q+fct2-fct1)
}

EDkhi2  <- function(x,D,N,logq) {
  # ---------------------------------------------------------------
  # FUNCTION
  #   Compute the penalty when L small
  # INPUT
  #   x scalar or vector
  #   D, N : positive integer
  #   logq : scalar
  # OUTPUT
  #   values in x : scalar or vector
  # CALLED BY
  #  EDkhi
 # ---------------------------------------------------------------
  t1 <- log(2*(2*x+N*D)/(N*(N+2)*x))-lbeta(1+D/2,N/2)
  t3 <- (N/2)*log(N/(N+x))+(D/2)*log(x/(N+x))
  return(logq-t1-t3)
}
