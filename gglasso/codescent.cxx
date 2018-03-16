#include <common.h>

float  softthresh(float z, float gamma);
float ObjValue(Eigen::VectorXf & R, // Nx1 residual vector.
	       Eigen::VectorXf & beta, 
	       Codescent_parameters & cdparam);

int NaiveUpdate(Eigen::VectorXf & Y, // Nx1 response vector.
		Eigen::MatrixXf & X, // NxP regressor matrix.
		Eigen::VectorXf & beta,
		Codescent_parameters & cdparam,
		unsigned verbose);

float NaiveUpdate_internal(Eigen::VectorXf & Y, // Nx1 response vector.
			   Eigen::MatrixXf & X, // NxP regressor matrix.
			   Eigen::VectorXf & beta,
			   Eigen::VectorXf & R, // Nx1 residual vector.
			   Codescent_parameters & cdparam);

int Updateresidual(Eigen::VectorXf & Y, // Nx1 response vector.
		   Eigen::MatrixXf & X, // NxP regressor matrix.
		   Eigen::VectorXf & beta, 
		   Eigen::VectorXf & R); // Nx1 residual vector.

int CovarianceUpdate(Eigen::VectorXf & Y, // Nx1 response vector.
		     Eigen::MatrixXf & X, // NxP regressor matrix.
		     Eigen::VectorXf & beta,
		     Codescent_parameters & cdparam,
		     unsigned verbose);

int CovarianceUpdate_internal(Eigen::VectorXf & XdotY,
			      Eigen::MatrixXf & XdotX,
			      Eigen::MatrixXf & X, // NxP regressor matrix.
			      Eigen::VectorXf & beta,
			      Codescent_parameters & cdparam,
			      unsigned verbose);

bool SweepCompleteSet(Eigen::VectorXf & XdotY,
		Eigen::MatrixXf & XdotX,
		Eigen::MatrixXf & X, // NxP regressor matrix.
		Eigen::VectorXf & grad,
		Eigen::VectorXf & beta,
		Eigen::VectorXf & actSet,
		Codescent_parameters & cdparam,
		      unsigned verbose);

bool SweepActiveSet(Eigen::VectorXf & XdotY,
		    Eigen::MatrixXf & XdotX,
		    Eigen::VectorXf & grad,
		    Eigen::VectorXf & beta,
		    Eigen::VectorXf & actSet,
		    Codescent_parameters & cdparam,
		    unsigned verbose);

int UpdateDotProd(Eigen::MatrixXf & XdotX, 
		  Eigen::MatrixXf & X, // NxP regressor matrix.
		  unsigned varIdx, 
		  unsigned numVar);

int UpdateGrad(Eigen::VectorXf & XdotY,
	       Eigen::MatrixXf & XdotX,
	       Eigen::VectorXf & grad,
	       Eigen::VectorXf & beta,
     	       Eigen::VectorXf & actSet);

/********************************************************/

int codescent(Eigen::VectorXf & Y, // Nx1 response vector.
	      Eigen::MatrixXf & X, // NxP regressor matrix.
	      Eigen::MatrixXf & beta,
	      Codescent_parameters & cdparam,
	      std::string method,
	      unsigned verbose)
{
     unsigned stepIdx = 0;
     Eigen::ArrayXf allLambda(cdparam.numSteps);
     Eigen::VectorXf curBeta(cdparam.numVar);
     curBeta.setZero();

     // gives a suggest initial value of lambda.
     float maxLambda = 0, thisMaxLambda = 0;
     for (unsigned varIdx = 0; varIdx < cdparam.numVar; varIdx ++) {
	  
	  thisMaxLambda = fabs(X.col(varIdx).dot(Y))/ (float)(cdparam.numObs * cdparam.alpha);
	  if (thisMaxLambda > maxLambda) {
	       maxLambda = thisMaxLambda;
	  }
     }

     printf("Suggest minimal lambda for initialization: %f.\n", maxLambda);

     // since we nee log(lambda), lambda need to be larger than zero.
     if (cdparam.finalLambda == 0) {
	  cdparam.finalLambda = 1e-37;
     }

     for (stepIdx = 0; stepIdx < cdparam.numSteps; stepIdx ++) {

	  allLambda[stepIdx] = exp( (log(cdparam.finalLambda) - log(cdparam.initLambda))*(stepIdx)/(cdparam.numSteps-1)  + log(cdparam.initLambda) );

     }

     if (verbose >= 1) {
	  std::cout << "All lambda array:\n" << allLambda << std::endl;
     }


     // run coordinate descent on this lambda value. warm start from previous
     // beta vector.
     for (stepIdx = 0; stepIdx < cdparam.numSteps; stepIdx ++) {
	  cdparam.lambda = allLambda[stepIdx];
	  if (method.compare("naive") == 0) {
	       NaiveUpdate(Y, X, curBeta, cdparam, verbose);
	  }
	  else if (method.compare("covupd") == 0) {
	       CovarianceUpdate(Y, X, curBeta, cdparam, verbose);
	  }
	  beta.col(stepIdx) = curBeta;
     }     
     
     return 0;
}

int NaiveUpdate(Eigen::VectorXf & Y, // Nx1 response vector.
		Eigen::MatrixXf & X, // NxP regressor matrix.
		Eigen::VectorXf & beta,
		Codescent_parameters & cdparam,
		unsigned verbose)
{
     Eigen::VectorXf R(cdparam.numObs);
     Updateresidual(Y, X, beta, R);

     unsigned sweepIdx = 0; // a sweep scan all X variables once. One step
			 // has many sweeps, and each step has same lambda.
     float deltaObjValue = 0;
     do {
	  // The naive update just update all variables once, and return the
	  // change of the obj function.
	  deltaObjValue = NaiveUpdate_internal(Y, X, beta, R, cdparam);
	  sweepIdx ++;
     } while (deltaObjValue > 0.0001 && sweepIdx < 5);

     return 0;
}

float NaiveUpdate_internal(Eigen::VectorXf & Y, // Nx1 response vector.
		  Eigen::MatrixXf & X, // NxP regressor matrix.
		  Eigen::VectorXf & beta,
		  Eigen::VectorXf & R, // Nx1 residual vector.
		  Codescent_parameters & cdparam)
{
     float lambda = cdparam.lambda;
     unsigned varIdx = 0;
     unsigned numObs = cdparam.numObs;
     unsigned numVar = cdparam.numVar;
     float newbeta = 0;
     float partialVar = 0;

     float oldObjValue = ObjValue(R, beta, cdparam);

     for (varIdx = 0; varIdx < numVar; varIdx ++) {
	  partialVar = (1/float(numObs)) * X.col(varIdx).dot(R) + beta(varIdx);
	  newbeta = softthresh(partialVar, lambda * cdparam.alpha)/ ( 1 + lambda*(1-cdparam.alpha) );

	  // if current coefficient changes, update coefficient vector, and also
	  // update residual to reflect this change.
	  if (fabs(newbeta - beta(varIdx)) > cdparam.eps ) {
	       beta(varIdx) = newbeta;
	       Updateresidual(Y, X, beta, R);
	  }
     } // varIdx

     // return the change of objective function.
     return fabs(oldObjValue - ObjValue(R, beta, cdparam));
     
}

int Updateresidual(Eigen::VectorXf & Y, // Nx1 response vector.
		   Eigen::MatrixXf & X, // NxP regressor matrix.
		   Eigen::VectorXf & beta, 
		   Eigen::VectorXf & R) // Nx1 residual vector.
{
     R = Y - X * beta;
     return 0;
}

float ObjValue(Eigen::VectorXf & R, // Nx1 residual vector.
	       Eigen::VectorXf & beta, 
	       Codescent_parameters & cdparam)
{
     float lambda = cdparam.lambda;

     // regularizer.
     float P = (1 - cdparam.alpha) * 0.5 * beta.squaredNorm() * cdparam.alpha * beta.lpNorm<1>();
     return (1/(float)(2*cdparam.numObs)) * R.squaredNorm() + lambda * P ;
}


float  softthresh(float z, float gamma)
{
     if (z > 0 && gamma < fabs(z) ) {
	  return (z - gamma);
     }
     else if (z < 0 && gamma < fabs(z) ) {
	  return (z + gamma);
     }
     else if (gamma >= fabs(z) ) {
	  return 0;
     }
     else {
	  // do nothing.
     }
     
     printf("softthresh(): reach end of the function. Something is wrong.\n");
     exit(1);
}
	      
int CovarianceUpdate(Eigen::VectorXf & Y, // Nx1 response vector.
		     Eigen::MatrixXf & X, // NxP regressor matrix.
		     Eigen::VectorXf & beta,
		     Codescent_parameters & cdparam,
		     unsigned verbose)
{
     unsigned varIdx = 0, varIdx2 = 0;
     unsigned numVar = cdparam.numVar;

     // Compute <x_j, y>
     Eigen::VectorXf XdotY(cdparam.numVar);
     for (varIdx = 0; varIdx < numVar; varIdx ++) {
	  XdotY(varIdx) = X.col(varIdx).dot(Y);
     }
     
     // compute <x_j, x_k>.
     Eigen::MatrixXf XdotX(cdparam.numVar, cdparam.numVar);
     for (varIdx = 0; varIdx < numVar; varIdx ++) {
	  for (varIdx2 = 0; varIdx2 <= varIdx; varIdx2 ++) {
	       XdotX(varIdx, varIdx2) = X.col(varIdx).dot(X.col(varIdx2));	       
	       XdotX(varIdx2, varIdx) = XdotX(varIdx, varIdx2);
	  }
     }
	       
     CovarianceUpdate_internal(XdotY, XdotX, X, beta, cdparam, verbose);

     return 0;
}

int CovarianceUpdate_internal(Eigen::VectorXf & XdotY,
			      Eigen::MatrixXf & XdotX,
			      Eigen::MatrixXf & X, // NxP regressor matrix.
			      Eigen::VectorXf & beta,
			      Codescent_parameters & cdparam,
			      unsigned verbose)
{

     unsigned actSweepIdx = 0, comSweepIdx = 0;; 

     // store the gradient: sum_i x_ij * r_i - \sum_{beta > 0} <x_j, x_k> *
     // beta_k. Also init beta to zero, since only when beta = 0, grad equals to
     // the value above.
     Eigen::VectorXf grad(cdparam.numVar);
     grad = XdotY;
     beta.setZero();

     // active set indicator and init to all zeros.
     Eigen::VectorXf actSet(cdparam.numVar);
     actSet.setZero();

     bool actSetConverge = 0;
     bool completeSetConverge = 0;

     // sweep complete set first.
     completeSetConverge = SweepCompleteSet(XdotY, XdotX, X, grad, beta, actSet, cdparam, verbose);
     while(!completeSetConverge) {
	  comSweepIdx ++;
	  actSweepIdx = 0;

	  // while(!actSetConverge) {
	  //      // sweep active set.
	  //      actSweepIdx ++;
	  //      actSetConverge = SweepActiveSet(XdotY, XdotX, grad, beta, actSet, cdparam, verbose);
	  //      if (verbose >= 2) {
	  // 	    printf("    Active set sweep %i.\n", actSweepIdx);
	  //      }
	  // };
	  
	  // sweep complete set.
	  completeSetConverge = SweepCompleteSet(XdotY, XdotX, X, grad, beta, actSet, cdparam, verbose);	  
	  if (verbose >= 1) {
	       printf("Complete set sweep %i.\n", comSweepIdx);
	  }
     }

     return 0;
}

// update variables  on whole set
bool SweepCompleteSet(Eigen::VectorXf & XdotY,
		Eigen::MatrixXf & XdotX,
		Eigen::MatrixXf & X, // NxP regressor matrix.
		Eigen::VectorXf & grad,
		Eigen::VectorXf & beta,
		Eigen::VectorXf & actSet,
		Codescent_parameters & cdparam,
		unsigned verbose)
{
     
     float lambda = cdparam.lambda;
     unsigned varIdx = 0;
     unsigned numVar = cdparam.numVar;
     float newbeta = 0;
     float partialVar = 0;
     bool converge = 1;


     for (varIdx = 0; varIdx < numVar; varIdx ++) {
	  partialVar = grad(varIdx) / cdparam.numObs + beta(varIdx);
	  newbeta = softthresh(partialVar, lambda * cdparam.alpha)/ ( 1 + lambda*(1-cdparam.alpha) );

	  if ( (fabs(newbeta) > cdparam.eps) && (fabs(beta(varIdx)) < cdparam.eps) ) {
	       // new variable enter the model.
	       actSet(varIdx) = 1;
	       if (verbose >= 1) {
		    printf("variable %i enter active set.\n", varIdx);
	       }
	  }

	  if (fabs(newbeta - beta(varIdx)) > cdparam.eps) {
	       // change of active variable. Update gradient.
	       beta(varIdx) = newbeta;

	       // updating gradient after updating beta.
	       UpdateGrad(XdotY, XdotX, grad, beta, actSet);
	       converge = 0;
		    
	       if (verbose >= 1) {
		    printf("    beta  %i is updated to %f.\n", varIdx, newbeta);
	       }
	  } // if (fabs
     } // varIdx
     return converge;
}


// update variables in active set.
bool SweepActiveSet(Eigen::VectorXf & XdotY,
		    Eigen::MatrixXf & XdotX,
		    Eigen::VectorXf & grad,
		    Eigen::VectorXf & beta,
		    Eigen::VectorXf & actSet,
		    Codescent_parameters & cdparam,
		    unsigned verbose)
{
     
     float lambda = cdparam.lambda;
     unsigned varIdx = 0;
     unsigned numVar = cdparam.numVar;
     float newbeta = 0;
     float partialVar = 0;

     bool converge = 1;

     for (varIdx = 0; varIdx < numVar; varIdx ++) {
	  // only update active set.
	  if (actSet(varIdx)) {
	       partialVar = grad(varIdx) / cdparam.numObs + beta(varIdx);
	       newbeta = softthresh(partialVar, lambda * cdparam.alpha)/ ( 1 + lambda*(1-cdparam.alpha) );

	       if (fabs(newbeta - beta(varIdx)) > cdparam.eps) {
		    // change of active variable. Update gradient.
		    beta(varIdx) = newbeta;
		    UpdateGrad(XdotY, XdotX, grad, beta, actSet);
		    converge = 0;
	       }	

	       if (fabs(newbeta) < cdparam.eps) {
		    // even variable in active set shink back to zero, do not
		    // deactivate it. Instead, keep it in active set. This may
		    // happen if coefficient go across zero to the othe other
		    // side of the y axis.
		    printf("warning: active coefficient %i shinks to %f.\n", varIdx, newbeta);
	       }
	  } // actSet
     } // varIdx

     return converge;
}


// new variable enter the model, need update XdotX. Populate <x_j, x_k>
int UpdateDotProd(Eigen::MatrixXf & XdotX, 
		  Eigen::MatrixXf & X, // NxP regressor matrix.
		  unsigned varIdx, 
		  unsigned numVar)
{
     unsigned k = 0;
     for (k = 0; k < numVar; k++) {
	  XdotX(varIdx, k) = X.col(varIdx).dot(X.col(k));
	  XdotX(k,varIdx) = XdotX(varIdx, k);

	  printf("XdotX(%i, %i) = %f, XdotX(%i, %i) = %f\n", varIdx, k, XdotX(varIdx, k), k, varIdx, XdotX(k, varIdx));
     }

     return 0;
}

int UpdateGrad(Eigen::VectorXf & XdotY,
	       Eigen::MatrixXf & XdotX,
	       Eigen::VectorXf & grad,
	       Eigen::VectorXf & beta,
     	       Eigen::VectorXf & actSet)

{
     unsigned numVar = XdotY.size();
     unsigned varIdx = 0, actIdx = 0;
     
     for (varIdx = 0; varIdx < numVar; varIdx ++) {
	  grad(varIdx) = XdotY(varIdx);

	  for (actIdx = 0; actIdx < numVar; actIdx ++) {
	       if (actSet(actIdx)) {
		    grad(varIdx) = grad(varIdx) - XdotX(varIdx, actIdx) * beta(actIdx);
	       }
	  } // end actIdx
     } // end varIdx.

     return 0;
}
