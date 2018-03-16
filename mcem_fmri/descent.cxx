#include "common.h"

double eval_dll(ImageType3D::Pointer samplePtr,
		int numClusters,
		double beta);

double descent(ImageType3D::Pointer samplePtr,
	       int numClusters,
	       double beta,
	       double t) // init step max step length.
{
     int S = 0, curS = 0; //  neighbors sum.
     double sumexps = 0, sumexp = 0;
     double dQ = 0; // derivative of log likelihood (prior).

     double delx = 1;
     int thisLabel = 0;
     int curLabel = 0;

     double alpha = 0.2;
     double rho = 0.2;
     double beta_old = -100;
     double beta_init = beta;

     ImageType3D::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;     
     ImageType3D::IndexType nidx;     

     // Descent iteration.

     while (abs((beta - beta_old)/beta) > 0.0001) {
	  // evaluate derivative of Q at beta.
	  dQ = -eval_dll(samplePtr, numClusters, beta);
	  // search direction.
	  delx = dQ>=0?(-1):1;

	  // backtracking line search, to find scale parameter t.

	  while (-eval_ll(samplePtr, numClusters, beta + t*delx)
		 > -eval_ll(samplePtr, numClusters, beta) 
		 + alpha * t * dQ * delx) {
	       t = rho * t;
	       printf("backtracking: t = %f\n", t);
	  }
	  
	  beta_old = beta;
	  beta = beta + t * delx;
	  printf("descent: beta = %f, dQ = %f\n", beta, dQ);
     }

     return beta;
}

// Evaluate derivative of prior log-likelihood at point beta.
double eval_dll(ImageType3D::Pointer samplePtr,
	       int numClusters,
	       double beta)
{
     int S = 0, curS = 0; //  neighbors sum.
     double sumexps = 0, sumexp = 0;
     double dQ = 0; // derivative of log likelihood (prior).
     int thisLabel = 0;
     int curLabel = 0;

     ImageType3D::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;     
     ImageType3D::IndexType nidx;     

     // compute derivative of Q function.
     dQ = 0;
     for (sampleIdx[2] = 0; sampleIdx[2] < sampleSize[2]; sampleIdx[2] ++) {
	  nidx[2] = sampleIdx[2];
	  for (sampleIdx[0] = 0; sampleIdx[0] < sampleSize[0]; sampleIdx[0] ++) {
	       for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1] ++) {
		    curLabel = samplePtr->GetPixel(sampleIdx);
		    sumexp = 0;
		    sumexps = 0;
		    for (thisLabel = 0; thisLabel < numClusters; thisLabel++) {
			 // Compute S_i
			 S = 0; 
			 if (sampleIdx[0] > 0) {
			      nidx[0] = sampleIdx[0] - 1;
			      nidx[1] = sampleIdx[1];
			      S = S + int(thisLabel != samplePtr->GetPixel(nidx));
			 }
			 
			 if  (sampleIdx[0] < sampleSize[0] - 1) {
			      nidx[0] = sampleIdx[0] +1;
			      nidx[1] = sampleIdx[1];
			      S = S + int(thisLabel != samplePtr->GetPixel(nidx));
			 }
			 
			 if  (sampleIdx[1] > 0) {
			      nidx[0] = sampleIdx[0];
			      nidx[1] = sampleIdx[1] - 1;
			      S = S + int(thisLabel != samplePtr->GetPixel(nidx));
			 }

			 if  (sampleIdx[1] < sampleSize[1] - 1) {
			      nidx[0] = sampleIdx[0];
			      nidx[1] = sampleIdx[1] + 1;
			      S = S + int(thisLabel != samplePtr->GetPixel(nidx));
			 }
			 sumexp = sumexp + exp(-beta * S);
			 sumexps = sumexps + S * exp(-beta * S);
			 if (thisLabel == curLabel) {
			      curS = S;
			 }
		    }
			 
		    dQ += (- curS + sumexps/sumexp)/double(sampleSize[2]);
	       }
	  }
     }
     return dQ;
}

