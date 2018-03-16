#include "commonalt.h"

#ifndef __COMMONALT_H__
#define __COMMONALT_H__
#warning "included more than once"


class MCModel
{
private:
     unsigned int m_numClusters;
     CompType * m_cls;
     double alpha_;
     double m_beta;
     unsigned int m_numScan;

     double priorll_, m_oldPriorll;
     double m_cll, m_oldCll;
     
     double m_oldBeta;
     CompType * m_oldCls;

     unsigned short m_verbose;

     vnl_vector<float> m_singleMCSampleLLNew;
     vnl_vector<float> m_singleMCSampleLLOld;
     
public:
     MCModel(ImageType2D::Pointer imagePtr,
	     ImageType3D::Pointer samplePtr,
	     unsigned int clusters,
	     unsigned int numScan,
	     double beta,
	     std::string samplingMethod,
	     unsigned short verbose);


     ~MCModel();
     int print(std::string format);
     int estep(ImageType2D::Pointer imagePtr,
	       ImageType3D::Pointer samplePtr,
	       unsigned int burnin);

     int estimatePriorParameter(double a, double b,
				ImageType3D::Pointer samplePtr);

     double eval_ll(ImageType3D::Pointer samplePtr,
		 double curBeta);

     int estimateMu(ImageType3D::Pointer samplePtr,
		    ImageType2D::Pointer imagePtr);

     int estimateSigma(ImageType3D::Pointer samplePtr,
		       ImageType2D::Pointer imagePtr);

     double eval_cll(ImageType2D::Pointer imagePtr,
		  ImageType3D::Pointer samplePtr);
     
     // Test the model convergence. Will return 1 if iteration
     // converge, otherwise return 0.
     int testConvergence(float maxRelativeErr);
     
     int updateLL(ImageType2D::Pointer imagePtr,
		  ImageType3D::Pointer samplePtr);
     
     int ll_over_beta(ImageType3D::Pointer samplePtr);

     int testClusterCenterMove();

     double kmeans(ImageType2D::Pointer imagePtr, ImageType3D::Pointer samplePtr);

     int initGaussianMixture(ImageType2D::Pointer imagePtr, 
			     ImageType3D::Pointer samplePtr, 
			     unsigned numRepeat);

     double estLabelsByMean(ImageType2D::Pointer imagePtr,
			    ImageType3D::Pointer samplePtr);

     int estPriorPramNewton(ImageType3D::Pointer samplePtr, unsigned maxNewtonIter);

     float evalFirstDerivative(Image4DChar::Pointer neighborSumPtr);

     float evalSecondDerivative(Image4DChar::Pointer neighborSumPtr);
     
     float evalDerivativeQutient(ImageType3D::Pointer samplePtr,
				 Image4DChar::Pointer neighborSumPtr,
				 float * firstDeriv,
				 float * secondDeriv);

     int evalNeighborSum(ImageType3D::Pointer samplePtr,
			 Image4DChar::Pointer neighborSumPtr);

     // Test first and second derivative over all beta.
     int testDerivatives(ImageType3D::Pointer samplePtr,
			 Image4DChar::Pointer neighborSumPtr);

     // Find labeling by ICM in E step of EM.
     int icmEstep(ImageType2D::Pointer imagePtr,
		  ImageType3D::Pointer samplePtr,
		  std::string scope,
		  unsigned int maxScan);


     // Sort the labels such that small labels have small mean. Also
     // change the labe image to be consistent with sorted Gaussian mu
     // and sigma.
     int sortLabels(ImageType3D::Pointer samplePtr);

     // Given parametes mu and sigma, have maximum likelihood estimate
     // of the labels, ignoring the smoothness(Markov) paior, and
     // write back to last slcie of samplePtr. This can be used to
     // init labels before EM starts. Or, it can be used to init ICM.
     int estLabelML(ImageType3D::Pointer samplePtr,
		    ImageType2D::Pointer imagePtr);

     // Given current labels as initial labeling, use simulated
     // annealing to get optimal labels. Save to last slice of
     // samplePtr.
     int annealing(ImageType2D::Pointer imagePtr, 
		   ImageType3D::Pointer samplePtr, 
		   unsigned maxRun);

     // Graphcuts to estimate labels, given all parameters beta, mu and sigma.
     int estimateLabelsGraphcuts(ImageType2D::Pointer imagePtr, 
				 ImageType3D::Pointer samplePtr);

     // Init parameters to true value.
     int initWithTruth(float sigma, float offset, float scale, float beta);

     // Draw samples from posterior p(f | d) by N independent chains. 
     int nChainSampling(ImageType2D::Pointer imagePtr,
			ImageType3D::Pointer samplePtr,
			unsigned int burnin);
     
     // Init (or re init) The N indepedent Markov chains by MLE.
     int nChainInit(ImageType3D::Pointer samplePtr, std::string method);

     // Ascent-based method to test EM convergence.
     int ascentConvergeTest(ImageType3D::Pointer samplePtr,float signLevel, float ASE, float minllinc);

     // Ascent-based method to test if more samples are needed. Will
     // restore all parameters, including prior and conditional LL to
     // old value. So, run this after ascentConvergenceTest().
     int ascentAcceptEstimation(ImageType3D::Pointer samplePtr,
				float signLevel,
				float ASE,
				float addSamplePercent);
     
     // Compute asymptotic standard errors.
     double computeASE(ImageType2D::Pointer imagePtr,
		       ImageType3D::Pointer samplePtr);

     int updateSingleSampleLL(ImageType2D::Pointer imagePtr,
			      ImageType3D::Pointer samplePtr);
};


#endif
