#include "commonalt.h"

#ifndef __COMMONALT_H__
#define __COMMONALT_H__
#warning "included more than once"


class MCModel
{
private:
     unsigned int m_numClusters;
     std::vector < CompType >  m_cls,  m_oldCls;
     double alpha_;
     double m_beta;
     unsigned int m_numScan;

//     double priorll_, m_oldPriorll;
//     double m_cll, m_oldCll;
     
     double m_oldBeta;

     unsigned short m_verbose;

     vnl_vector<float> m_priorLLNew;
     vnl_vector<float> m_condLLNew;

     vnl_vector<float> m_priorLLOld;
     vnl_vector<float> m_condLLOld;

     // neighborSum pointer, for log likelihood evaluation and also
     // for estimating beta.
     ImageType5DChar::Pointer m_neighborSumPtr;

public:
     MCModel(ImageType4DFloat::Pointer imagePtr,
	     ImageType4DChar::Pointer samplePtr,
	     unsigned int clusters,
	     unsigned int numScan,
	     double meanKappa,
	     double stdKappa,
	     std::string samplingMethod,
	     unsigned short verbose);


     ~MCModel();
     int printSelf(std::string format);

     // For pixels outside of the mask, we set it to -1. So later
     // routine will know this is outside of mask. No mask need later!
     int applyMask(ImageType4DChar::Pointer samplePtr,
		   ImageType3DFloat::Pointer maskPtr,
		   float maskThreshhold);


     int estep(ImageType4DChar::Pointer samplePtr,
		   ImageType4DFloat::Pointer imagePtr,
		   unsigned int burnin);

     int estimatePriorParameter(double a, double b,
     				ImageType4DChar::Pointer samplePtr,
				ImageType4DFloat::Pointer imagePtr);
     
     // evaluate prior log-likelihood for all samples, given beta.
     double evalPriorLL(ImageType4DChar::Pointer samplePtr, float beta);

     /* double eval_ll(ImageType3D::Pointer samplePtr, */
     /* 		 double curBeta); */

     // Given labels, estimate mu and numPoints.
     int estimateMu(ImageType4DChar::Pointer samplePtr,
		    ImageType4DFloat::Pointer imagePtr,
		    unsigned beginSampleIdx,
		    unsigned endSampleIdx);

     // Given labels and mean direction, estimate kappa. This is the
     // posterior mode of kappa given prior of kappa.
     int estimateKappaWithPrior(double meanKappa, 
				  double stdKappa, 
				  float timeSeriesLength,
				std::string initMethod);
     
     // wrapper for estimating mu and sigma. Can be called by other
     // non member functions.
     int estimateMuSigmaWrapper(ImageType4DChar::Pointer samplePtr,
				ImageType4DFloat::Pointer imagePtr,
				double meanKappa,
				double stdKappa);

     // wrapper of mu and sigma (or kappa) estimation.
     int MStep(ImageType4DChar::Pointer samplePtr,
	       ImageType4DFloat::Pointer imagePtr);

     /* double eval_cll(ImageType2D::Pointer imagePtr, */
     /* 		  ImageType3D::Pointer samplePtr); */
     
     /* // Test the model convergence. Will return 1 if iteration */
     /* // converge, otherwise return 0. */
     /* int testConvergence(float maxRelativeErr); */
     
     /* int updateLL(ImageType2D::Pointer imagePtr, */
     /* 		  ImageType3D::Pointer samplePtr); */
     
     /* int ll_over_beta(ImageType3D::Pointer samplePtr); */

     int testClusterCenterMove(float eps);

     double kmeans(ImageType4DFloat::Pointer imagePtr, 
		   ImageType4DChar::Pointer samplePtr);


     double estLabelsByMean(ImageType4DFloat::Pointer imagePtr,
     			    ImageType4DChar::Pointer samplePtr);


     int estPriorPramNewton(ImageType4DChar::Pointer samplePtr, 
			    ImageType4DFloat::Pointer imagePtr,
			    unsigned maxNewtonIter);


     float evalDerivativeQutient(ImageType4DChar::Pointer samplePtr,
				 float * firstDeriv,
				 float * secondDeriv);     


     int evalNeighborSum(ImageType4DChar::Pointer samplePtr);

     /* // Test first and second derivative over all beta. */
     /* int testDerivatives(ImageType3D::Pointer samplePtr, */
     /* 			 Image4DChar::Pointer neighborSumPtr); */

     /* // Find labeling by ICM in E step of EM. */

     int icmEstep(ImageType4DChar::Pointer samplePtr,
		  ImageType4DFloat::Pointer imagePtr,
		  std::string scope,
		  unsigned int maxScan);

     /* // Sort the labels such that small labels have small mean. Also */
     /* // change the labe image to be consistent with sorted Gaussian mu */
     /* // and sigma. */
     /* int sortLabels(ImageType3D::Pointer samplePtr); */

     // Given parametes mu and sigma, have maximum likelihood estimate
     // of the labels, ignoring the smoothness(Markov) paior, and
     // write back to last slcie of samplePtr. This can be used to
     // init labels before EM starts. Or, it can be used to init ICM.
     int estLabelML(ImageType4DChar::Pointer samplePtr,
		    ImageType4DFloat::Pointer imagePtr);

     /* // given current labels as initial labeling, use simulated */
     /* // annealing to get optimal labels. Save to last slice of */
     /* // samplePtr. */
     /* int annealing(ImageType2D::Pointer imagePtr,  */
     /* 		   ImageType3D::Pointer samplePtr,  */
     /* 		   unsigned maxRun); */

     /* // Graphcuts to estimate labels, given all parameters beta, mu and sigma. */
     /* int estimateLabelsGraphcuts(ImageType2D::Pointer imagePtr,  */
     /* 				 ImageType3D::Pointer samplePtr); */

     // Init parameters to true value.
     int initWithTruth(std::string cheatfile,
		       std::string trueimage,
		       ImageType4DChar::Pointer samplePtr,
		       ImageType4DFloat::Pointer imagePtr);

     // Draw samples from posterior p(f | d) by N independent chains.
     int nChainSampling(ImageType4DChar::Pointer samplePtr,
			ImageType4DFloat::Pointer imagePtr,
			unsigned int burnin);
     
     // Init (or re init) The N indepedent Markov chains by MLE.
     int nChainInit(ImageType4DChar::Pointer samplePtr, 
		    std::string method);

     /* // Ascent-based method to test EM convergence. */
     /* int ascentConvergeTest(ImageType3D::Pointer samplePtr,float signLevel, float ASE, float minllinc); */

     /* // Ascent-based method to test if more samples are needed. Will */
     /* // restore all parameters, including prior and conditional LL to */
     /* // old value. So, run this after ascentConvergenceTest(). */
     /* int ascentAcceptEstimation(ImageType3D::Pointer samplePtr, */
     /* 				float signLevel, */
     /* 				float ASE, */
     /* 				float addSamplePercent); */
     
     /* // Compute asymptotic standard errors. */
     /* double computeASE(ImageType2D::Pointer imagePtr, */
     /* 		       ImageType3D::Pointer samplePtr); */


     int updateSampleLL(ImageType4DFloat::Pointer imagePtr,
			    ImageType4DChar::Pointer samplePtr);

     // Init the cluster mean for K-means, by Furthest-First rule.
     int furthestFirstInit(ImageType4DFloat::Pointer imagePtr,
			   ImageType4DChar::Pointer samplePtr);

     // init kmeans by kmeans++ method. 
     int kmeanspp(ImageType4DFloat::Pointer imagePtr,
		  ImageType4DChar::Pointer samplePtr);

     int printMeanKappa();
};


#endif
