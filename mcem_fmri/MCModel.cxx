#include "commonalt.h"
#include "MCModel.h"
#include "GCoptimization.h"

#include <boost/math/distributions/students_t.hpp>
extern twister_base_gen_type mygenerator;

MCModel::MCModel(ImageType4DFloat::Pointer imagePtr,
		 ImageType4DChar::Pointer samplePtr,
		 unsigned int numClusters = 2,
		 unsigned int numScan = 10,
		 double meanKappa = 300,
		 double stdKappa = 1,
		 std::string samplingMethod = "mc",
		 unsigned short verbose = 0)


{
     if (numClusters <= 0) {
	  std::cout << "MCModel::MCModel: number of clusters must be greater than zero. Force set it to 2.\n";
	  m_numClusters = 2;
     }
     else {
	  m_numClusters = numClusters;
     }
     

     // init memberers.
     m_beta = 0;
     m_verbose = verbose;
     m_numScan = numScan;     

     ImageType4DFloat::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();
     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();

     ImageType4DFloat::IndexType imageIdx;     
     ImageType4DChar::IndexType sampleIdx;     

     // Apply mask, so labels outside the mask will be negative.
     // applyMask(samplePtr, maskPtr, maskThresh);

		    // // DBG
		    // saveimage4dchar(samplePtr, "samples.nii");
		    // saveExtractVolume(samplePtr, 0, "onesample.nii");
		    // // end of DBG

     // we define a bestCls to save the best mu so far.
     std::vector < CompType >  bestCls(m_numClusters);

     // init cls
     m_cls.resize(m_numClusters);
     m_oldCls.resize(m_numClusters);

     unsigned clsIdx = 0;
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_cls[clsIdx].mu.set_size(imageSize[3]);
	  m_oldCls[clsIdx].mu.set_size(imageSize[3]);
	  bestCls[clsIdx].mu.set_size(imageSize[3]);
     }


     // init m_singleMCSampleLLNew
     m_priorLLNew.set_size(m_numScan);
     m_priorLLNew.fill(0);
     m_condLLNew.set_size(m_numScan);
     m_condLLNew.fill(0);

     m_priorLLOld.set_size(m_numScan);
     m_priorLLOld.fill(0);
     m_condLLOld.set_size(m_numScan);
     m_condLLOld.fill(0);

     // Allocate memory for neighboring sum variable S.
     m_neighborSumPtr = ImageType5DChar::New();
     ImageType5DChar::SizeType neighborSumSize;
     neighborSumSize[0] = sampleSize[0];
     neighborSumSize[1] = sampleSize[1];
     neighborSumSize[2] = sampleSize[2];
     // Alloate memory for max number of MC samples.
     neighborSumSize[3] = sampleSize[3]; 
     neighborSumSize[4] = m_numClusters;
     ImageType5DChar::IndexType neighborSumIdx;
     neighborSumIdx.Fill(0);

     ImageType5DChar::RegionType neighborSumRegion;
     neighborSumRegion.SetSize(neighborSumSize);
     neighborSumRegion.SetIndex(neighborSumIdx);
     m_neighborSumPtr->SetRegions(neighborSumRegion);
     m_neighborSumPtr->Allocate();
     m_neighborSumPtr->FillBuffer(0);     

     unsigned repeatIdx = 0; // repeat index of K-means.
     double minMeanSSR = 1e10, meanSSR = 0;

     // // Init Gaussian Mixture parameters mu, numPoints and sigma
     // // Repeat kmeans clustering a few times and choose the best one
     // // by sum of square error.
     // for (repeatIdx = 0; repeatIdx < kmeansIter; repeatIdx ++) {
     // 	  printf("kmeans iteration %i begin:\n", repeatIdx);

     // 	  meanSSR = kmeans(imagePtr, samplePtr);

     // 	  if (meanSSR < minMeanSSR) {
     // 	       minMeanSSR = meanSSR;
     // 	       // Save best mu in to bestCls
     // 	       bestCls = m_cls;

     // 	       if (m_verbose >=1) {
     // 		    printf("initClusterCenter: meanSSR %f < minMeanSSR . best cluster centers are:\n", meanSSR);
     // 		    printSelf("normal");
     // 	       }
     // 	  }

     // }

     // // Found the best mean vectors. Restore the best mean vectors and numPoints
     // // into cls.
     // m_cls = bestCls;

     // // Given mu, estimate labels. Save to last slice of samplePtr.
     // estLabelsByMean(imagePtr, samplePtr);     

     // // back up current parameters.
     // m_oldCls = m_cls;

     // estimate mu.

     printSelf("normal");
     estimateMu(samplePtr, imagePtr, sampleSize[3]-1, sampleSize[3]-1);
     printSelf("normal");

     // Init Kappa. For big stdKappa, we hae diffuse prior and we use
     // Banerjee's method for estimation. For small stdKappa, we have
     // strong prior, the init value of kappa should be close to meanKappa.
     if (stdKappa > 2) {
	  estimateKappaWithPrior(meanKappa, 1e10, imageSize[3], "banerjee");
     }
     else {
	  for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	       m_cls[clsIdx].kappa = meanKappa;
	  }
     }

     // MLE: estimate labels and save in last MC samples.
     estLabelML(samplePtr, imagePtr);

     // Manually set beta to a small value. This is not different with
     // estimating beta using labels from MLE above. Anyway when EM
     // starts, beta will be estimated from MC (or ICM) samples. Only
     // difference is one the first iteration of E step sampling, we
     // ignore (because very small beta) spatial smoothness prior, and
     // sample from posterior is equivalent to sampling from ML.


     // // With known parameters, Use ICM to estimate labels and
     // // save it in last slice of samplePtr.
     // icmEstep(imagePtr, samplePtr, "lastone", 20);

     // Init N indepedent Markov chains with same results from MLE or ICM.
     if (samplingMethod.compare("nchain") == 0) {
     	  nChainInit(samplePtr, "ml");
     }



}
     
MCModel::~MCModel()
{

}

int MCModel::printSelf(std::string format)
{
     unsigned short clsIdx;
//     double 
     double priorLLNew = m_priorLLNew.mean();
     double condLLNew = m_condLLNew.mean();
     double priorLLOld = m_priorLLOld.mean();
     double condLLOld = m_condLLOld.mean();

     if (format.compare("normal") == 0) {
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
		    printf("cluster[%2i]: # of points: %ld, kappa = %4.2f, mu = ", clsIdx, m_cls[clsIdx].numPoints, m_cls[clsIdx].kappa); 
		    printVnlVector(m_cls[clsIdx].mu, 3);
	       }
	       printf("beta = %8f,   prior ll = %E,   cond ll = %E,   joint ll = %E, old joint ll = %E\n", m_beta, priorLLNew, condLLNew, priorLLNew + condLLNew, priorLLOld + condLLOld);
	  }
     else if (format.compare("table") == 0) {
	  printf("%2.4f %.2f %.2f %.2f   %.2f %.2f %.2f ", m_beta, priorLLNew, condLLNew, priorLLNew + condLLNew, priorLLOld, condLLOld, priorLLOld +  condLLOld);
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
		    printf("%6.4f %10i  ", m_cls[clsIdx].kappa, m_cls[clsIdx].numPoints);
	       }
	       printf("\n");
	  }
     else {
	  printf("MCModel::print(): print format not recognized.\n");
     }
     
     return 0;
}

int MCModel::applyMask(ImageType4DChar::Pointer samplePtr,
		       ImageType3DFloat::Pointer maskPtr,
		       float maskThreshhold)
{
     
     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();

     // mask.
     ImageType3DFloat::RegionType maskRegion;
     maskRegion = maskPtr->GetLargestPossibleRegion();
     ImageType3DFloat::SizeType maskSize = maskRegion.GetSize();
     IteratorType3DFloat maskIt(maskPtr, maskRegion);
     ImageType3DFloat::IndexType maskIdx;    

     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  maskIdx = maskIt.GetIndex();
	  sampleIdx[0] = maskIdx[0];
	  sampleIdx[1] = maskIdx[1];
	  sampleIdx[2] = maskIdx[2];

	  // For pixels outside of the mask, we set it to -1. So later
	  // routine will know this is outside of mask. No mask need later!
	  if (maskIt.Get() > maskThreshhold) {
	       for (sampleIdx[3] = 0; sampleIdx[3] < sampleSize[3]; sampleIdx[3] ++) {
		    samplePtr->SetPixel(sampleIdx, 0);
	       }
	  }
	  else {
	       for (sampleIdx[3] = 0; sampleIdx[3] < sampleSize[3]; sampleIdx[3] ++) {
		    samplePtr->SetPixel(sampleIdx, -1);
	       }
	  }
     } // end for
}

// Given labels, estimate mu and numPoints.
int MCModel::estimateMu(ImageType4DChar::Pointer samplePtr,
			     ImageType4DFloat::Pointer imagePtr,
			     unsigned beginSampleIdx,
			     unsigned endSampleIdx)
{
     unsigned int clsIdx = 0;


     // images
     ImageType4DFloat::IndexType imageIdx;    
     ImageType4DFloat::RegionType imageRegion;
     imageRegion = imagePtr->GetLargestPossibleRegion();
     ImageType4DFloat::SizeType imageSize = imageRegion.GetSize();
     IteratorType4DFloat imageIt(imagePtr, imageRegion);

     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     // define sample region.
     sampleIdx[0] = 0;
     sampleIdx[1] = 0;
     sampleIdx[2] = 0;
     sampleIdx[3] = beginSampleIdx;
     sampleSize[3] = endSampleIdx - beginSampleIdx + 1;

     ImageType4DChar::RegionType sampleRegion;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);
     
     IteratorType4DChar sampleIt(samplePtr, sampleRegion);

     vnl_vector <float> timeSeries(imageSize[3], 0);

     // reset all mu to zero.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_cls[clsIdx].numPoints = 0;
	  m_cls[clsIdx].mu = 0;
     }

     for (sampleIt.GoToBegin(); !sampleIt.IsAtEnd(); ++ sampleIt) {
	  if (sampleIt.Get() >= 0) {
	       sampleIdx = sampleIt.GetIndex();

	       imageIdx[0] = sampleIdx[0];
	       imageIdx[1] = sampleIdx[1];
	       imageIdx[2] = sampleIdx[2];

	       // First save time series at this position into a vector.
	       for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
		    timeSeries[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
	       }

	       clsIdx = sampleIt.Get();
	       m_cls[clsIdx].mu += timeSeries;
	       m_cls[clsIdx].numPoints ++;
	  } 
     }

     // compute mu_k and meanNorm.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {

	  // meanNorm see page 1350 of "Clustering on the Unit
	  // Hypersphere using von Mises-Fisher distributions" by
	  // Banerjee.
	  m_cls[clsIdx].meanNorm = m_cls[clsIdx].mu.two_norm() / m_cls[clsIdx].numPoints;
	  // m_cls[clsIdx].mu = m_cls[clsIdx].mu/m_cls[clsIdx].numPoints;
	  m_cls[clsIdx].mu.normalize();
     }


     return 0;
}


int MCModel::estimateKappaWithPrior(double meanKappa, 
				  double stdKappa, 
				  float timeSeriesLength,
				  std::string initMethod)
{

     // estimate kappa with prior information. Assume kappa also a
     // random variable, and find the posterior mode of kappa,
     // i.e. P(kappa | data) \prop P(kappa) x P(data | kappa). prior
     // is assume to be Gaussian with mu as meanKappa, and std as
     // stdKappa.

     int clsIdx = 0;
     using boost::math::cyl_bessel_i;

     float kappa = 0;
     double Ap = 0;

     float kappaOldMethod = 0;
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  if (initMethod.compare("current") == 0) {
	       // Just using current kappa as init value.
	       if (m_verbose >= 1) {
		    printf("estimateKappaWithPrior(): cls[%i]. Use current value to init kappa.  ", clsIdx);
	       }
	  }
	  else if (initMethod.compare("banerjee") == 0) {
	       if (m_verbose >= 1) {
		    printf("estimateKappaWithPrior(): cls[%i]. Use Banerjee to init kappa.  ", clsIdx);
	       }

	       m_cls[clsIdx].kappa = m_cls[clsIdx].meanNorm * (timeSeriesLength - m_cls[clsIdx].meanNorm * m_cls[clsIdx].meanNorm) / (1 - m_cls[clsIdx].meanNorm * m_cls[clsIdx].meanNorm);
	  }

	  if (m_verbose >= 1) {
	       printf("kappa_0 = %f, ", m_cls[clsIdx].kappa);
	  }

	  // first Newton, and get K1.
	  Ap = exp(logBesselI(timeSeriesLength/2, double(m_cls[clsIdx].kappa)) - logBesselI(timeSeriesLength/2 - 1, double(m_cls[clsIdx].kappa)));

	  ///// Debug ////
	  kappaOldMethod = m_cls[clsIdx].kappa - (Ap - m_cls[clsIdx].meanNorm)/(1 - Ap*Ap - ( (timeSeriesLength - 1) / (m_cls[clsIdx].kappa) )*Ap);
	  if (m_verbose >= 1) {
	       printf("kappa_1_noPrior= %f, ", kappaOldMethod);
	  }	       
	  ///// End of debug ////

	  kappa = m_cls[clsIdx].kappa;

	  m_cls[clsIdx].kappa = kappa - (Ap - m_cls[clsIdx].meanNorm + (kappa - meanKappa)/(m_cls[clsIdx].numPoints * stdKappa * stdKappa) ) 
	       / (1 - Ap*Ap - ( (timeSeriesLength - 1) / kappa ) * Ap + 1/(m_cls[clsIdx].numPoints * stdKappa*stdKappa) );

	  if (m_verbose >= 1) {
	       printf("kappa_1 = %f, ", m_cls[clsIdx].kappa);
	  }

	  // 2nd Newton, and get K2.
	  Ap = exp(logBesselI(timeSeriesLength/2, double(m_cls[clsIdx].kappa)) - logBesselI(timeSeriesLength/2 - 1, double(m_cls[clsIdx].kappa)));

	  kappa = m_cls[clsIdx].kappa;
	  m_cls[clsIdx].kappa = kappa - (Ap - m_cls[clsIdx].meanNorm + (kappa - meanKappa)/(m_cls[clsIdx].numPoints * stdKappa * stdKappa) ) 
	       / (1 - Ap*Ap - ( (timeSeriesLength - 1) / kappa ) * Ap + 1/(m_cls[clsIdx].numPoints * stdKappa*stdKappa) );


	  if (m_verbose >= 1) {
	       printf("kappa_2 = %f.\n ", m_cls[clsIdx].kappa);

	  }


     } // clsIdx


}

// wrapper of mu and sigma (or kappa) estimation.
int MCModel::estimateMuSigmaWrapper(ImageType4DChar::Pointer samplePtr,
				    ImageType4DFloat::Pointer imagePtr,
				    double meanKappa,
				    double stdKappa)
{
     ImageType4DFloat::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();
     
     // back up current parameters.
     m_oldCls = m_cls;

     // estimate mu and kappa.
     estimateMu(samplePtr, imagePtr, 0, m_numScan - 1);

     // kappa should not be too big or too small. 
     if (stdKappa > 2) {
	  // diffuse prior. Use Banerjee estimation as init value for
	  // Newton.
	  estimateKappaWithPrior(meanKappa, stdKappa, imageSize[3], "banerjee");
     }
     else {
	  // strong prior. Use current value as prior.
	  estimateKappaWithPrior(meanKappa, stdKappa, imageSize[3], "current");	 } 

     // update conditional log-likelihood P(d | f).
     updateSampleLL(imagePtr, samplePtr);
}


double MCModel::evalPriorLL(ImageType4DChar::Pointer samplePtr, float beta)
{
     vnl_vector<float> priorLL;
     priorLL.set_size(m_numScan);
     unsigned mcSampleIdx = 0;
     for(mcSampleIdx = 0; mcSampleIdx < m_numScan; mcSampleIdx ++) {
	  priorLL[mcSampleIdx] = singleSamplePriorLL(samplePtr,
						     m_neighborSumPtr,
						     mcSampleIdx,
						     beta,
						     m_numClusters,
						     m_numScan);
     }
     return priorLL.mean();
}

// Maximization prior log-likelihood over beta. return beta.
int  MCModel::estimatePriorParameter(double a, double b,
				     ImageType4DChar::Pointer samplePtr,
				     ImageType4DFloat::Pointer imagePtr)

{
     if (a >= b) {
	  printf("estimatePriorParameter: a must be smaller than b.\n");
	  exit(1);
     }


     double tau = (sqrt(5) - 1)/2;
     double x1 = 0, x2 = 0;
     double f1 = 0, f2 = 0;
     x1 = a + (1 - tau) * (b - a);
     f1 =  - evalPriorLL(samplePtr, x1);
     x2 = a + tau * (b - a);

     // compute the minimal value of negagive log-likelihood.
     f2 =  - evalPriorLL(samplePtr, x2);
     while (fabs(b - a) > 0.001) {
	  if (f1 > f2) {
	       a = x1;
	       x1 = x2;
	       f1 = f2;
	       x2 = a + tau * (b - a);
	       f2 =  - evalPriorLL(samplePtr, x2);
	  }
	  else {
	       b = x2;
	       x2 = x1;
	       f2 = f1;
	       x1 = a + (1 - tau) * (b - a);
	       f1 =  - evalPriorLL(samplePtr, x1);
	  }
	  printf("a = %.3f, x1 = %.3f, x2 = %.3f, b = %.3f, f1 = %f, f2 = %f\n", a, x1, x2, b, f1, f2);
     }
     // Save old beta in m_oldBeta.
     m_oldBeta = m_beta;
     m_beta = (a + b)/2;


     // Update prior log-likelihood.
     updateSampleLL(imagePtr, samplePtr);     

     // // DBG
     // for (float fake_beta = 0.01; fake_beta < 30; fake_beta += 0.1)
     // {
     // 	  printf("fake_beta = %f, prior LL = %f.\n", fake_beta, evalPriorLL(samplePtr, fake_beta));
     // }
     
     return 0;
}


// Init cluster center by K-means.
double MCModel::kmeans(ImageType4DFloat::Pointer imagePtr, 
		       ImageType4DChar::Pointer samplePtr)

{
     unsigned int clsIdx = 0;

     // images
     ImageType4DFloat::IndexType imageIdx;    
     ImageType4DFloat::RegionType imageRegion;
     imageRegion = imagePtr->GetLargestPossibleRegion();
     ImageType4DFloat::SizeType imageSize = imageRegion.GetSize();
     IteratorType4DFloat imageIt(imagePtr, imageRegion);

     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();

     vnl_vector <float> timeSeries(imageSize[3], 0);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);
	  
     // // Randomly init the cluster center for K-means clustering.
     // sampleIdx[3] = 0;
     // for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
     // 	  // Randomly find a point in the mask.
     // 	  do {
     // 	       sampleIdx[0] = floor(uni() * sampleSize[0]);
     // 	       sampleIdx[1] = floor(uni() * sampleSize[1]);
     // 	       sampleIdx[2] = floor(uni() * sampleSize[2]);
     // 	  }
     // 	  while(samplePtr->GetPixel(sampleIdx) < 0);
     // 	  // OK I find the point, and assign it to cluster center.
     // 	  imageIdx[0] = sampleIdx[0];
     // 	  imageIdx[1] = sampleIdx[1];
     // 	  imageIdx[2] = sampleIdx[2];
     // 	  for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
     // 	       m_cls[clsIdx].mu[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
     // 	  }
     // }

     // Init kmeans by kmeans++ method.
     kmeanspp(imagePtr, samplePtr);

     double meanSSR = 0;
     do {
	  // update label assignment.
	  meanSSR = estLabelsByMean(imagePtr, samplePtr);

	  // back up current parameters.
	  m_oldCls = m_cls;

	  // estimate mu.
	  estimateMu(samplePtr, imagePtr, sampleSize[3]-1, sampleSize[3]-1);

     }
     while(testClusterCenterMove(1e-6));
     
     return meanSSR;
}

int MCModel::furthestFirstInit(ImageType4DFloat::Pointer imagePtr,
				ImageType4DChar::Pointer samplePtr)
{

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);


     unsigned int clsIdx = 0;

     // images
     ImageType4DFloat::IndexType imageIdx;    
     ImageType4DFloat::RegionType imageRegion;
     imageRegion = imagePtr->GetLargestPossibleRegion();
     ImageType4DFloat::SizeType imageSize = imageRegion.GetSize();
     IteratorType4DFloat imageIt(imagePtr, imageRegion);

     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::RegionType sampleRegion;

     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     // we only iterate over the last volume of the samplePtr, so we
     // need define the region such that it points to the last volume.
     sampleIdx[0] = 0;
     sampleIdx[1] = 0;
     sampleIdx[2] = 0;
     sampleIdx[3] = sampleSize[3] - 1;

     sampleSize[3] = 1;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);

     IteratorType4DChar sampleIt(samplePtr, sampleRegion);

     vnl_vector <float> timeSeries(imageSize[3], 0);

     double minDistance = 0;
     double maxMinDistance = 0;
     double thisDistance = 0;
     int prevIdx = 0;
     ImageType4DChar::IndexType bestSampleIdx;

     // Find a random bogus center.
     sampleIdx[3] = 0;
     do {
	  sampleIdx[0] = floor(uni() * sampleSize[0]);
	  sampleIdx[1] = floor(uni() * sampleSize[1]);
	  sampleIdx[2] = floor(uni() * sampleSize[2]);
     }
     while(samplePtr->GetPixel(sampleIdx) < 0);
     // OK I find the point, and assign it to mu_0
     imageIdx[0] = sampleIdx[0];
     imageIdx[1] = sampleIdx[1];
     imageIdx[2] = sampleIdx[2];
     for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
	  m_cls[0].mu[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
     }

     // Find all cluster center by Furthest-First rule. This also
     // includes first cluster. At last we throw away the bogus
     // center.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  // assume we have 1:clsIdx-1 mu available.
	  maxMinDistance = -1e10;
	  for (sampleIt.GoToBegin(); !sampleIt.IsAtEnd(); ++ sampleIt) {
	       if (sampleIt.Get() >= 0) {

		    // First save time series at this position into a vector.
		    sampleIdx = sampleIt.GetIndex();
		    imageIdx[0] = sampleIdx[0];
		    imageIdx[1] = sampleIdx[1];
		    imageIdx[2] = sampleIdx[2];	  
		    for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
			 timeSeries[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
		    }
	       
		    minDistance = 1e10;
		    if (clsIdx == 0) {
			 minDistance = (timeSeries - m_cls[0].mu).two_norm();
		    }
		    else {
			 for (prevIdx = 0; prevIdx < clsIdx; prevIdx ++) {
			      thisDistance = (timeSeries - m_cls[prevIdx].mu).two_norm();
			      if (thisDistance < minDistance) {
				   minDistance = thisDistance;
			      }
			 }
		    }
		    // got the minDistance for current points. Now check if
		    // this is max among all points.
		    if (minDistance > maxMinDistance) {
			 maxMinDistance = minDistance;
			 bestSampleIdx = sampleIdx;
		    }

	       } // in mask.
	  } // sampleIt.

	  // OK I find the best point, and assign mu to that point.
	  imageIdx[0] = bestSampleIdx[0];
	  imageIdx[1] = bestSampleIdx[1];
	  imageIdx[2] = bestSampleIdx[2];
	  for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
	       m_cls[clsIdx].mu[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
	  }

     } // clsIdx

}
	       

int MCModel::kmeanspp(ImageType4DFloat::Pointer imagePtr,
		      ImageType4DChar::Pointer samplePtr)
{

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);


     unsigned int clsIdx = 0;
     float randNum = 0;
     float cdfValue = 0;
     float sumOfPdf = 0;

     // images
     ImageType4DFloat::IndexType imageIdx;    
     ImageType4DFloat::RegionType imageRegion;
     imageRegion = imagePtr->GetLargestPossibleRegion();
     ImageType4DFloat::SizeType imageSize = imageRegion.GetSize();

     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::RegionType sampleRegion;

     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();

     // Allocate memory for pdf of each data point.

     ImageType3DFloat::Pointer pdfPtr = ImageType3DFloat::New();
     ImageType3DFloat::IndexType pdfStart;
     ImageType3DFloat::IndexType pdfIdx;
     pdfStart.Fill(0);
     ImageType3DFloat::SizeType pdfSize;
     pdfSize[0] = sampleSize[0];
     pdfSize[1] = sampleSize[1];
     pdfSize[2] = sampleSize[2];

     ImageType3DFloat::RegionType pdfRegion;
     pdfRegion.SetSize(pdfSize);
     pdfRegion.SetIndex(pdfStart);
     pdfPtr->SetRegions(pdfRegion);
     pdfPtr->Allocate();
     pdfPtr->FillBuffer(0);

     // we only iterate over the last volume of the samplePtr, so we
     // need define the region such that it points to the last volume.
     sampleIdx[0] = 0;
     sampleIdx[1] = 0;
     sampleIdx[2] = 0;
     sampleIdx[3] = sampleSize[3] - 1;

     sampleSize[3] = 1;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);

     IteratorType4DChar sampleIt(samplePtr, sampleRegion);

     vnl_vector <float> timeSeries(imageSize[3], 0);

     double minDistance = 0;
     double thisDistance = 0;
     int prevIdx = 0;

     // Find a random bogus center.
     sampleIdx[3] = 0;
     do {
	  sampleIdx[0] = floor(uni() * sampleSize[0]);
	  sampleIdx[1] = floor(uni() * sampleSize[1]);
	  sampleIdx[2] = floor(uni() * sampleSize[2]);
     }
     while(samplePtr->GetPixel(sampleIdx) < 0);
     // OK I find the point, and assign it to mu_0
     imageIdx[0] = sampleIdx[0];
     imageIdx[1] = sampleIdx[1];
     imageIdx[2] = sampleIdx[2];
     for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
	  m_cls[0].mu[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
     }

     // Find all cluster center by Furthest-First rule. This also
     // includes first cluster. At last we throw away the bogus
     // center.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  // assume we have 1:clsIdx-1 mu available.
	  sumOfPdf = 0;

	  for (sampleIt.GoToBegin(); !sampleIt.IsAtEnd(); ++ sampleIt) {
	       if (sampleIt.Get() >= 0) {

		    // First save time series at this position into a vector.
		    sampleIdx = sampleIt.GetIndex();
		    imageIdx[0] = sampleIdx[0];
		    imageIdx[1] = sampleIdx[1];
		    imageIdx[2] = sampleIdx[2];	  
		    for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
			 timeSeries[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
		    }
	       
		    minDistance = 1e10;
		    if (clsIdx == 0) {
			 minDistance = (timeSeries - m_cls[0].mu).two_norm();
		    }
		    else {
			 for (prevIdx = 0; prevIdx < clsIdx; prevIdx ++) {
			      thisDistance = (timeSeries - m_cls[prevIdx].mu).two_norm();
			      if (thisDistance < minDistance) {
				   minDistance = thisDistance;
			      }
			 }
		    }
		    
		    // Got the minDistance for current point. Save it in pdfPtr.
		    pdfIdx[0] = sampleIdx[0];
		    pdfIdx[1] = sampleIdx[1];
		    pdfIdx[2] = sampleIdx[2];
		    pdfPtr->SetPixel(pdfIdx, minDistance);
		    sumOfPdf += minDistance;
	       } // in mask.
	  } // sampleIt.

	  // normalize pdf.
	  for (sampleIt.GoToBegin(); !sampleIt.IsAtEnd(); ++ sampleIt) {
	       if (sampleIt.Get() >= 0) {
		    sampleIdx = sampleIt.GetIndex();
		    pdfIdx[0] = sampleIdx[0];
		    pdfIdx[1] = sampleIdx[1];
		    pdfIdx[2] = sampleIdx[2];
		    pdfPtr->SetPixel(pdfIdx, pdfPtr->GetPixel(pdfIdx) / sumOfPdf);
	       }
	  }

	  randNum = uni();
	  cdfValue = 0;
	  sampleIt.GoToBegin();
	       
	  while (randNum > cdfValue) {
	       while(sampleIt.Get() < 0) {
		    ++ sampleIt;
	       }
	       sampleIdx = sampleIt.GetIndex();

	       pdfIdx[0] = sampleIdx[0];
	       pdfIdx[1] = sampleIdx[1];
	       pdfIdx[2] = sampleIdx[2];
	       cdfValue += pdfPtr->GetPixel(pdfIdx);

	       ++ sampleIt;
	  }
	  
	  // found the data point.
	  imageIdx[0] = pdfIdx[0];
	  imageIdx[1] = pdfIdx[1];
	  imageIdx[2] = pdfIdx[2];

	  for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
	       m_cls[clsIdx].mu[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
	  }

     } // clsIdx

}
	       


double MCModel::estLabelsByMean(ImageType4DFloat::Pointer imagePtr,
				ImageType4DChar::Pointer samplePtr)

{
     unsigned int clsIdx = 0;

     // images
     ImageType4DFloat::IndexType imageIdx;    
     ImageType4DFloat::RegionType imageRegion;
     imageRegion = imagePtr->GetLargestPossibleRegion();
     ImageType4DFloat::SizeType imageSize = imageRegion.GetSize();
     IteratorType4DFloat imageIt(imagePtr, imageRegion);

     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::RegionType sampleRegion;

     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     // we only iterate over the last volume of the samplePtr, so we
     // need define the region such that it points to the last volume.
     sampleIdx[0] = 0;
     sampleIdx[1] = 0;
     sampleIdx[2] = 0;
     sampleIdx[3] = sampleSize[3] - 1;

     sampleSize[3] = 1;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);

     IteratorType4DChar sampleIt(samplePtr, sampleRegion);


     unsigned nearestClsLabel = 0;
     unsigned timeIdx = 0;
     double nearestDistance = 1000000;
     float  thisDistance = 0;
     double meanSSR = 0;
     unsigned long numAllPoints = 0;

     vnl_vector <float> timeSeries(imageSize[3], 0);
//     float pixelIntensity = 0;

     
     for (sampleIt.GoToBegin(); !sampleIt.IsAtEnd(); ++ sampleIt) {
	  if (sampleIt.Get() >= 0) {
	       numAllPoints ++;
	       // First save time series at this position into a vector.
	       sampleIdx = sampleIt.GetIndex();
	       imageIdx[0] = sampleIdx[0];
	       imageIdx[1] = sampleIdx[1];
	       imageIdx[2] = sampleIdx[2];	  
	       for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
		    timeSeries[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
	       }
	  
	       nearestDistance = 1000000;
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
		    thisDistance = (timeSeries - m_cls[clsIdx].mu).two_norm();
		    if (thisDistance < nearestDistance) {
			 nearestClsLabel = clsIdx;
			 nearestDistance = thisDistance;
		    }
	       }

	       // Found the nearest cluster label.
	       samplePtr->SetPixel(sampleIdx, nearestClsLabel);
	       meanSSR += nearestDistance;
	  }
     }
     meanSSR = meanSSR / numAllPoints;
     return meanSSR;
}

int MCModel::testClusterCenterMove(float eps)
{
     unsigned int clsIdx = 0;
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
	  if ((m_cls[clsIdx].mu - m_oldCls[clsIdx].mu).two_norm() > eps) {
	       return 0;
	  }
     }

     return 1;
}

// Compute the qutient of the first and second derivative.
float MCModel::evalDerivativeQutient(ImageType4DChar::Pointer samplePtr,
				     float * firstDeriv,
				     float * secondDeriv)
{
     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     // define sample region.
     sampleIdx.Fill(0);
     sampleSize[3] = m_numScan;
     ImageType4DChar::RegionType sampleRegion;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);
     
     IteratorType4DChar sampleIt(samplePtr, sampleRegion);

     double MZero = 0, MOne = 0, MTwo = 0;

     double expBetaSum = 0;

     ImageType5DChar::IndexType neighborSumIdx;
     unsigned int clsIdx = 0;
     unsigned char thisNeighborSum = 0;

     * firstDeriv = 0;
     * secondDeriv = 0;

     for (sampleIt.GoToBegin(); !sampleIt.IsAtEnd(); ++ sampleIt) {
	  if (sampleIt.Get() >= 0) {
	       sampleIdx = sampleIt.GetIndex();
	       
	       neighborSumIdx[0] = sampleIdx[0];
	       neighborSumIdx[1] = sampleIdx[1];
	       neighborSumIdx[2] = sampleIdx[2];
	       neighborSumIdx[3] = sampleIdx[3];
	       
	       // Compute M0, M1, and M2
	       MZero = 0, MOne = 0, MTwo = 0;
	       for (neighborSumIdx[4] = 0; neighborSumIdx[4] < m_numClusters; neighborSumIdx[4]++) {

		    thisNeighborSum = m_neighborSumPtr->GetPixel(neighborSumIdx);
		    expBetaSum = exp(- m_beta * thisNeighborSum);
		    MZero += expBetaSum;
		    expBetaSum = expBetaSum * thisNeighborSum; // Now it's delta M1.
		    MOne += expBetaSum;
		    expBetaSum = expBetaSum * thisNeighborSum; // Now it's delta M2.
		    MTwo += expBetaSum;
	       }

	       neighborSumIdx[4] = sampleIt.Get();
	       thisNeighborSum = m_neighborSumPtr->GetPixel(neighborSumIdx);

	       // now the number of MC samples is m_clusters, not the
	       // 3rd dim of neighboruSum data block. This is because
	       // neighborSum is allocated with 3rd dim as the max
	       // number of MC samples.
	       *firstDeriv += (double)(-thisNeighborSum + MOne/MZero)/ m_numScan;
	       *secondDeriv += (double)(MOne*MOne - MZero * MTwo)/(double)(MZero * MZero * m_numScan);

	  } // in mask.
     } // iterator
     return (*firstDeriv) / (*secondDeriv);
}

int MCModel::evalNeighborSum(ImageType4DChar::Pointer samplePtr)

{
     signed short currentLabel = 0;
     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::RegionType sampleRegion;
     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     sampleIdx.Fill(0);
     // only sample on first m_numScan volumes.
     sampleSize[3] = m_numScan;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);

     // Define neighborhood iterator
     typedef itk::ConstantBoundaryCondition< ImageType4DChar >  BoundaryConditionType;
     BoundaryConditionType constCondition;
     constCondition.SetConstant(-1);
     typedef itk::NeighborhoodIterator< ImageType4DChar, BoundaryConditionType > NeighborhoodIteratorType;
     
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);

     NeighborhoodIteratorType neighborSampleIt( radius, samplePtr, sampleRegion );
     neighborSampleIt.OverrideBoundaryCondition(&constCondition);

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1, 0}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1, 0}};

     ImageType5DChar::IndexType neighborSumIdx;
     unsigned int neighborSum = 0;
     
     signed thisLabel = 0;
     for (neighborSampleIt.GoToBegin(); !neighborSampleIt.IsAtEnd(); ++ neighborSampleIt) {
	  currentLabel = neighborSampleIt.GetCenterPixel();
	  if (currentLabel >= 0) {
	       sampleIdx = neighborSampleIt.GetIndex();

	       neighborSumIdx[0] = sampleIdx[0];
	       neighborSumIdx[1] = sampleIdx[1];	  
	       neighborSumIdx[2] = sampleIdx[2];
	       neighborSumIdx[3] = sampleIdx[3];

	       for (thisLabel = 0; thisLabel < m_numClusters; thisLabel++) {
		    neighborSum = 0;
		    if (neighborSampleIt.GetPixel(xminus) >= 0) {
			 neighborSum += int(thisLabel != neighborSampleIt.GetPixel(xminus));
		    }
		    
		    if (neighborSampleIt.GetPixel(xplus) >= 0) {
			 neighborSum += int(thisLabel != neighborSampleIt.GetPixel(xplus));
		    }
		    
		    if (neighborSampleIt.GetPixel(yminus) >= 0) {
			 neighborSum += int(thisLabel != neighborSampleIt.GetPixel(yminus));
		    }

		    if (neighborSampleIt.GetPixel(yplus) >= 0) {
			 neighborSum += int(thisLabel != neighborSampleIt.GetPixel(yplus));
		    }

		    if (neighborSampleIt.GetPixel(zminus) >= 0) {
			 neighborSum += int(thisLabel != neighborSampleIt.GetPixel(zminus));
		    }

		    if (neighborSampleIt.GetPixel(zplus) >= 0) {
			 neighborSum += int(thisLabel != neighborSampleIt.GetPixel(zplus));
		    }

		    neighborSumIdx[4] = thisLabel;
		    m_neighborSumPtr->SetPixel(neighborSumIdx, neighborSum);
	       } // thisLabel
	  } // in mask
     } // iterator.
     return 0;
}

int MCModel::estPriorPramNewton(ImageType4DChar::Pointer samplePtr, 
				ImageType4DFloat::Pointer imagePtr,
				unsigned maxNewtonIter)

{

     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     
     // Newton's method.
     float stepLength = 1000;
     float firstDeriv = 0, secondDeriv = 0;
     unsigned newtonIter = 0;

     // save beta into oldBeta.
     m_oldBeta = m_beta;
     do 
     {
	  evalDerivativeQutient(samplePtr, &firstDeriv, &secondDeriv);
	  stepLength = firstDeriv/secondDeriv;
	  
	  if (m_verbose >= 1) {
	       printf("estPriorPramNewton(): current beta = %f, stepLength = %f, 1st deriv = %f, 2nd deriv = %f.\n", m_beta, stepLength, firstDeriv, secondDeriv);
	  }
//	  m_beta = m_beta - vnl_math_sgn(stepLength) * vnl_math_min(vnl_math_abs(stepLength), float(0.1));

	  m_beta = m_beta -stepLength;

	  newtonIter++;
     }
     while(fabs(stepLength) > 0.0001 && newtonIter < maxNewtonIter);
     
     // 
     updateSampleLL(imagePtr, samplePtr);

     return (0);
}


int MCModel::icmEstep(ImageType4DChar::Pointer samplePtr,
		      ImageType4DFloat::Pointer imagePtr,
		      std::string scope,
		      unsigned int maxScan)

{
     using boost::math::cyl_bessel_i;

     unsigned int scan = 0;
     double denergy = 0;
     double smallestDenergy = 10000;
     unsigned clsIdx = 0;
     signed char bestLabel = 0;

     unsigned anyChangePixel = 0;
     int currentLabel = 0;


     // images
     ImageType4DFloat::IndexType imageIdx;    
     ImageType4DFloat::RegionType imageRegion;
     imageRegion = imagePtr->GetLargestPossibleRegion();
     ImageType4DFloat::SizeType imageSize = imageRegion.GetSize();
     IteratorType4DFloat imageIt(imagePtr, imageRegion);

     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::RegionType sampleRegion;
     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     sampleIdx.Fill(0);

     // only sample on first m_numScan volumes.
     if (scope.compare("lastone") == 0) {
	  sampleIdx[3] = sampleSize[3] - 1;
	  sampleSize[3] = 1;
     }
     else if (scope.compare("all") == 0) {
	  sampleIdx[3] = 0;
	  sampleSize[3] = m_numScan;
     }
     else {
	  printf("icmEstep(): scope must be either all or lastone. exit.\n");
	  exit(1);
     }

     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);
     IteratorType4DChar sampleIt(samplePtr, sampleRegion);

     vnl_vector <float> timeSeries(imageSize[3], 0);

     // Define neighborhood iterator
     typedef itk::ConstantBoundaryCondition< ImageType4DChar >  BoundaryConditionType;
     BoundaryConditionType constCondition;
     constCondition.SetConstant(-1);
     typedef itk::NeighborhoodIterator< ImageType4DChar, BoundaryConditionType > NeighborhoodIteratorType;
     
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);

     NeighborhoodIteratorType neighborSampleIt( radius, samplePtr, sampleRegion);
     neighborSampleIt.OverrideBoundaryCondition(&constCondition);

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1, 0}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1, 0}};


     ///////        Define random number generator.        ////////
     // Uniform integer generator.
     boost::uniform_int<> uni_int(1, m_numClusters); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // compute normalization constant of von Mises
     // Fisher. vmfLogConst[i] is the log (c_d (kappa) ) for the i'th
     // clusters. See page 1350 of "Clustering on the Unit-Sphere
     // using von Mises Fisher..."
     vnl_vector<double> vmfLogConst(m_numClusters);
     float myD = imageSize[3];
     double const Pi = 4 * atan(1);
     double allNumPoints = 0;
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
//	  vmfLogConst[clsIdx] = (myD/2 - 1) * log (m_cls[clsIdx].kappa) - myD/2 * log(2*Pi) -  log(cyl_bessel_i(double(myD/2 - 1), double(m_cls[clsIdx].kappa)));
	  vmfLogConst[clsIdx] = (myD/2 - 1) * log (m_cls[clsIdx].kappa) - myD/2 * log(2*Pi) -  logBesselI(myD/2 - 1, m_cls[clsIdx].kappa);
	  allNumPoints += m_cls[clsIdx].numPoints;

//	  printf("icmEstep(): clsIdx[%i]: cyl = %f, approx = %f\n", clsIdx, log(cyl_bessel_i(double(myD/2 - 1), double(m_cls[clsIdx].kappa))), logBesselI(myD/2 - 1, m_cls[clsIdx].kappa));
     }


     // for ICM on last volume, we only need count points in last
     // volume.
     if (scope.compare("lastone") == 0) {
	  allNumPoints = allNumPoints/m_numScan;
     }

     // sampling Markov Random Field. 

     scan = 0;

     do {
	  anyChangePixel = 0;
	  for (neighborSampleIt.GoToBegin(); !neighborSampleIt.IsAtEnd(); ++ neighborSampleIt) {
	       currentLabel = neighborSampleIt.GetCenterPixel();
	       if (currentLabel >= 0) {
		    smallestDenergy = 1e10;
		    for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {

			 denergy 
			      = int(clsIdx != neighborSampleIt.GetPixel(xminus))
			      - int(currentLabel != neighborSampleIt.GetPixel(xminus))

			      + int(clsIdx != neighborSampleIt.GetPixel(xplus)) 
			      - int(currentLabel != neighborSampleIt.GetPixel(xplus))

			      + int(clsIdx != neighborSampleIt.GetPixel(yminus)) 
			      - int(currentLabel != neighborSampleIt.GetPixel(yminus))

			      + int(clsIdx != neighborSampleIt.GetPixel(yplus)) 
			      - int(currentLabel != neighborSampleIt.GetPixel(yplus))

			      + int(clsIdx != neighborSampleIt.GetPixel(zminus)) 
			      - int(currentLabel != neighborSampleIt.GetPixel(zminus))

			      + int(clsIdx != neighborSampleIt.GetPixel(zplus)) 
			      - int(currentLabel != neighborSampleIt.GetPixel(zplus));

			 denergy = m_beta * denergy;

			 // likelihood energy.
			 sampleIdx = neighborSampleIt.GetIndex();
			 imageIdx[0] = sampleIdx[0];
			 imageIdx[1] = sampleIdx[1];
			 imageIdx[2] = sampleIdx[2];

			 // First save time series at this position into a vector.
			 for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
			      timeSeries[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
			 }

			 // The energy = - Log(likelihood).
			 denergy = denergy
			      + (- vmfLogConst[clsIdx]
				 - m_cls[clsIdx].kappa * inner_product(timeSeries, m_cls[clsIdx].mu))

			      - (- vmfLogConst[currentLabel]
				 - m_cls[currentLabel].kappa * inner_product(timeSeries, m_cls[currentLabel].mu));

			      if (denergy < smallestDenergy) {
				   bestLabel = clsIdx;
				   smallestDenergy = denergy;
			      }
		    } // end clsIdx

		    // By ICM, we just assign the bestLabel to current pixel.
		    if (bestLabel != currentLabel) {
			 neighborSampleIt.SetCenterPixel(bestLabel);
			 anyChangePixel ++;
		    }
	       } // in mask.
	  } // iterator.
	  scan ++;

	  if (m_verbose >= 1) {
	       printf("icmEstep(): scan = %i, anyChaingePixel = %i.\n", scan, anyChangePixel);
	  }

     } // scans
     while (anyChangePixel > allNumPoints * 0.005  && scan <= maxScan);

     // Since labels have changed. we need update whatever that
     // depends on labels, including the neighborSum data structure.

     // Compute neighborSum.
     evalNeighborSum(samplePtr);
     // compute prior LL and cond LL. In fact only prior LL
     // changed. Anyway...
     updateSampleLL(imagePtr, samplePtr);

     return (0);
}


int MCModel::estLabelML(ImageType4DChar::Pointer samplePtr,
			ImageType4DFloat::Pointer imagePtr)
{
     using boost::math::cyl_bessel_i;
     unsigned int clsIdx = 0;

     // images
     ImageType4DFloat::IndexType imageIdx;    
     ImageType4DFloat::RegionType imageRegion;
     imageRegion = imagePtr->GetLargestPossibleRegion();
     ImageType4DFloat::SizeType imageSize = imageRegion.GetSize();
     IteratorType4DFloat imageIt(imagePtr, imageRegion);

     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType4DChar::RegionType sampleRegion;
     // we only iterate over the last volume of the samplePtr, so we
     // need define the region such that it points to the last volume.
     sampleIdx[0] = 0;
     sampleIdx[1] = 0;
     sampleIdx[2] = 0;
     sampleIdx[3] = sampleSize[3] - 1;

     sampleSize[3] = 1;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);

     IteratorType4DChar sampleIt(samplePtr, sampleRegion);


     // compute normalization constant of von Mises
     // Fisher. vmfLogConst[i] is the log (c_d (kappa) ) for the i'th
     // clusters. See page 1350 of "Clustering on the Unit-Sphere
     // using von Mises Fisher..."
     vnl_vector<double> vmfLogConst(m_numClusters);
     float myD = imageSize[3];
     double const Pi = 4 * atan(1);

     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
//	  vmfLogConst[clsIdx] = (myD/2 - 1) * log (m_cls[clsIdx].kappa) - myD/2 * log(2*Pi) -  log(cyl_bessel_i(double(myD/2 - 1), double(m_cls[clsIdx].kappa)));

	  vmfLogConst[clsIdx] = (myD/2 - 1) * log (m_cls[clsIdx].kappa) - myD/2 * log(2*Pi) -  logBesselI(myD/2 - 1, m_cls[clsIdx].kappa);

//	  printf("estLbelML(): clsIdx[%i]: cyl = %f, approx = %f\n", clsIdx, log(cyl_bessel_i(double(myD/2 - 1), double(m_cls[clsIdx].kappa))), logBesselI(myD/2 - 1, m_cls[clsIdx].kappa));

     }

     vnl_vector <float> timeSeries(imageSize[3], 0);

     signed char bestLabel = 0;
     float thisLL = 0,  largestLL= -100000;
     sampleIdx[3] = sampleSize[3] - 1; // write to last volume.

     for (sampleIt.GoToBegin(); !sampleIt.IsAtEnd(); ++ sampleIt) {
	  if (sampleIt.Get() >= 0) {
	       sampleIdx = sampleIt.GetIndex();

	       // First save time series at this position into a vector.
	       imageIdx[0] = sampleIdx[0];
	       imageIdx[1] = sampleIdx[1];
	       imageIdx[2] = sampleIdx[2];	  
	       for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
		    timeSeries[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
	       }

	       // Which cluster is best for current time series.
	       largestLL = -1000000;
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
		    thisLL = vmfLogConst[clsIdx]
			 + m_cls[clsIdx].kappa * inner_product(timeSeries, m_cls[clsIdx].mu);
		    if (thisLL > largestLL) {
			 largestLL = thisLL;
			 bestLabel = clsIdx;
		    }
	       }

	       // write back.
	       samplePtr->SetPixel(sampleIdx, bestLabel);
	  }
     }

     // Since labels have changed. we need update whatever that
     // depends on labels, including the neighborSum data structure.

     // compute prior LL and cond LL. In fact only prior LL
     // changed. Anyway...
     updateSampleLL(imagePtr, samplePtr);

     return (0);
}

int MCModel::initWithTruth(std::string cheatfile,
			   std::string trueimage,
			   ImageType4DChar::Pointer samplePtr,
			   ImageType4DFloat::Pointer imagePtr)
{

     unsigned clsIdx = 0;
     std::ifstream fileStream(cheatfile.c_str(), std::ifstream::in);
     vnl_vector<float>::iterator it;

     float temp = 0;
     if (fileStream.is_open())
     {
	  // First line is beta
	  fileStream >> m_beta;
	  for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	       // first element , kappa
	       fileStream >> temp;
	       m_cls[clsIdx].kappa = temp;
	       for (it = m_cls[clsIdx].mu.begin(); it < m_cls[clsIdx].mu.end(); it++) {
		    
		    fileStream >> * it;
	       }
	  }
	  
     }
     else {
	  std::cout << "Unable to open file: " << cheatfile << std::endl;
	  exit(1);
     }


     // read in true label image
     ReaderType3DChar::Pointer trueImageReader = ReaderType3DChar::New();
     trueImageReader->SetFileName(trueimage);
     trueImageReader->Update();
     ImageType3DChar::Pointer trueImagePtr = trueImageReader->GetOutput();

     ImageType3DChar::RegionType trueImageRegion;
     trueImageRegion = trueImageReader->GetOutput()->GetLargestPossibleRegion();
     ImageType3DChar::SizeType trueImageSize = trueImageRegion.GetSize();
     IteratorType3DChar trueImageIt(trueImagePtr, trueImageRegion);
     ImageType3DChar::IndexType trueImageIdx;


     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType4DChar::SizeType regionSize;
     // define sample region.
     sampleIdx[0] = 0;
     sampleIdx[1] = 0;
     sampleIdx[2] = 0;
     sampleIdx[3] = sampleSize[3] - 1;

     regionSize = sampleSize;
     regionSize[3] = 1;

     ImageType4DChar::RegionType sampleRegion;
     sampleRegion.SetSize(regionSize);
     sampleRegion.SetIndex(sampleIdx);
     
     IteratorType4DChar sampleIt(samplePtr, sampleRegion);

     save3dchar(trueImagePtr, "mytrueth.nii");

     for (sampleIt.GoToBegin(), trueImageIt.GoToBegin(); !sampleIt.IsAtEnd(); ++ sampleIt, ++ trueImageIt) {
	  if (sampleIt.Get() >= 0) {
	       sampleIdx = sampleIt.GetIndex();
	       trueImageIdx = trueImageIt.GetIndex();

	       sampleIt.Set(trueImageIt.Get());
	  }
     }


     // After init all parameters (kappa, beta, mu) with truth, we
     // need to estimate labels.

     // estimate mu.
     estimateMu(samplePtr, imagePtr, sampleSize[3]-1, sampleSize[3]-1);

     // Given mu, estimate labels. Save to last slice of samplePtr.
     estLabelsByMean(imagePtr, samplePtr);     

     // MLE: estimate labels and save in last MC samples.
     estLabelML(samplePtr, imagePtr);

     // back up current parameters.
     m_oldCls = m_cls;

     nChainInit(samplePtr, "ml");	  

}


int MCModel::nChainInit(ImageType4DChar::Pointer samplePtr, 
			std::string method)
{

     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType4DChar::IndexType sampleIdx;     
     ImageType4DChar::IndexType sampleDestIdx;          
     ImageType4DChar::RegionType sampleRegion;


     sampleIdx[0] = 0;
     sampleIdx[1] = 0;
     sampleIdx[2] = 0;

     vnl_vector<float> singleLL;
     singleLL.set_size(m_numScan);

     // choose which volume is the source sample.
     if (method.compare("ml") == 0) {
	  sampleIdx[3] = sampleSize[3] - 1;
     }
     else if (method.compare("largestll") == 0) {
	  // choose the one with largest likelihood.
	  singleLL = m_priorLLNew + m_condLLNew;
	  unsigned short maxLLIdx = 0;
	  double maxLL = -10000000;
	  for (unsigned mcSampleIdx = 0; mcSampleIdx < m_numScan; mcSampleIdx++) {
	       if (singleLL[mcSampleIdx] > maxLL) {
		    maxLL = singleLL[mcSampleIdx];
		    maxLLIdx = mcSampleIdx;
	       }
	  }
	  sampleIdx[3] = maxLLIdx;

	  if (m_verbose >=1) {
	       printf("nChainInit(): MC sample %i is used for init value of all samples. \n", sampleIdx[3]);
	  }
     }
     else {
	  printf("nChainInit(): option must be either ml or largestll. Exit().");
	  exit(1);
     }

     sampleSize[3] = 1;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);
     IteratorType4DChar sampleIt(samplePtr, sampleRegion);
     
     for (sampleIt.GoToBegin(); !sampleIt.IsAtEnd(); ++ sampleIt) {
	  if (sampleIt.Get() >= 0) {
	       sampleIdx = sampleIt.GetIndex();

	       sampleDestIdx[0] = sampleIdx[0];
	       sampleDestIdx[1] = sampleIdx[1];
	       sampleDestIdx[2] = sampleIdx[2];

	       for (sampleDestIdx[3] = 0; sampleDestIdx[3] < m_numScan; sampleDestIdx[3] ++) {
		    samplePtr->SetPixel(sampleDestIdx, sampleIt.Get());
	       } // end time series.
	  } // end mask > 0
     } // end all maskIt

     // DBG
     saveimage4dchar(samplePtr, "samples.nii");

     return (0);
}

int MCModel::nChainSampling(ImageType4DChar::Pointer samplePtr,
			    ImageType4DFloat::Pointer imagePtr,
			    unsigned int burnin)
{

     using boost::math::cyl_bessel_i;

     int scan = 0;
     double p_acpt = 0;
     double denergy = 0;
     int cand = 1;
     int currentLabel = 0;
     unsigned clsIdx = 0;

     // images
     ImageType4DFloat::IndexType imageIdx;    
     ImageType4DFloat::RegionType imageRegion;
     imageRegion = imagePtr->GetLargestPossibleRegion();
     ImageType4DFloat::SizeType imageSize = imageRegion.GetSize();
     IteratorType4DFloat imageIt(imagePtr, imageRegion);

     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::RegionType sampleRegion;
     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     sampleIdx.Fill(0);
     // only sample on first m_numScan volumes.
     sampleSize[3] = m_numScan;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);
     IteratorType4DChar sampleIt(samplePtr, sampleRegion);

     vnl_vector <float> timeSeries(imageSize[3], 0);

     // Define neighborhood iterator
     typedef itk::ConstantBoundaryCondition< ImageType4DChar >  BoundaryConditionType;
     BoundaryConditionType constCondition;
     constCondition.SetConstant(-1);
     typedef itk::NeighborhoodIterator< ImageType4DChar, BoundaryConditionType > NeighborhoodIteratorType;
     
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);

     NeighborhoodIteratorType neighborSampleIt( radius, samplePtr, sampleRegion );
     neighborSampleIt.OverrideBoundaryCondition(&constCondition);

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1, 0}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1, 0}};


     ///////        Define random number generator.        ////////
     // Uniform integer generator.
     boost::uniform_int<> uni_int(0, m_numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     // compute normalization constant of von Mises
     // Fisher. vmfLogConst[i] is the log (c_d (kappa) ) for the i'th
     // clusters. See page 1350 of "Clustering on the Unit-Sphere
     // using von Mises Fisher..."
     vnl_vector<double> vmfLogConst(m_numClusters);
     float myD = imageSize[3];
     double const Pi = 4 * atan(1);

     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
//	  vmfLogConst[clsIdx] = (myD/2 - 1) * log (m_cls[clsIdx].kappa) - myD/2 * log(2*Pi) -  log(cyl_bessel_i(double(myD/2 - 1), double(m_cls[clsIdx].kappa)));

	  vmfLogConst[clsIdx] = (myD/2 - 1) * log (m_cls[clsIdx].kappa) - myD/2 * log(2*Pi) -  logBesselI(myD/2 - 1, m_cls[clsIdx].kappa);

//	  printf("nChainSampling(): clsIdx[%i]: cyl = %f, approx = %f\n", clsIdx, log(cyl_bessel_i(double(myD/2 - 1), double(m_cls[clsIdx].kappa))), logBesselI(myD/2 - 1, m_cls[clsIdx].kappa));

     }


     // sampling Markov Random Field. 
     scan = 0;
     for (scan = 0; scan < burnin + 1; scan ++) {
	  for (neighborSampleIt.GoToBegin(); !neighborSampleIt.IsAtEnd(); ++ neighborSampleIt) {
	       currentLabel = neighborSampleIt.GetCenterPixel();
	       if (currentLabel >= 0) {

		    cand = roll_die();
		    denergy 
			 = int(cand != neighborSampleIt.GetPixel(xminus))
			 - int(currentLabel != neighborSampleIt.GetPixel(xminus))

			 + int(cand != neighborSampleIt.GetPixel(xplus)) 
			 - int(currentLabel != neighborSampleIt.GetPixel(xplus))

			 + int(cand != neighborSampleIt.GetPixel(yminus)) 
			 - int(currentLabel != neighborSampleIt.GetPixel(yminus))

			 + int(cand != neighborSampleIt.GetPixel(yplus)) 
			 - int(currentLabel != neighborSampleIt.GetPixel(yplus))

			 + int(cand != neighborSampleIt.GetPixel(zminus)) 
			 - int(currentLabel != neighborSampleIt.GetPixel(zminus))

			 + int(cand != neighborSampleIt.GetPixel(zplus)) 
			 - int(currentLabel != neighborSampleIt.GetPixel(zplus));

		    denergy = m_beta * denergy;

		    // likelihood energy.
		    sampleIdx = neighborSampleIt.GetIndex();
		    imageIdx[0] = sampleIdx[0];
		    imageIdx[1] = sampleIdx[1];
		    imageIdx[2] = sampleIdx[2];

		    // First save time series at this position into a vector.
		    for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
			 timeSeries[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
		    }


		    // The energy = - Log(likelihood).
		    denergy = denergy
			 + (- vmfLogConst[cand] - m_cls[cand].kappa * inner_product(timeSeries, m_cls[cand].mu))

			 - (- vmfLogConst[currentLabel] - m_cls[currentLabel].kappa * inner_product(timeSeries, m_cls[currentLabel].mu));
		    
		    // if energy change less than zero, just accept
		    // candidate. otherwise accept with exp(- energy
		    // change).
			 
		    if (denergy <= 0) {
			 neighborSampleIt.SetCenterPixel(cand);
		    }
		    else {
			 p_acpt = exp(-denergy);
			 if (uni() < p_acpt) {
			      neighborSampleIt.SetCenterPixel(cand);
			 }
		    }
	       } // in mask
	  } // iterators.
	  if (m_verbose >= 1) {
	       printf("nChainSamping(): scan %u done.\n", scan);
	  }

     } // for scan
	  

	  
     // Since labels have changed. we need update whatever that
     // depends on labels, including the neighborSum data structure.

     // Compute neighborSum.
     evalNeighborSum(samplePtr);
     // compute prior LL and cond LL. In fact only prior LL
     // changed. Anyway...
     updateSampleLL(imagePtr, samplePtr);

     // DBG
     saveimage4dchar(samplePtr, "samples.nii");
     saveExtractVolume(samplePtr, sampleSize[3] - 1, "onesample.nii");

     // end of DBG
     return 0;    
}

int MCModel::updateSampleLL(ImageType4DFloat::Pointer imagePtr,
			    ImageType4DChar::Pointer samplePtr)
{
     unsigned clsIdx = 0;
     unsigned mcSampleIdx = 0;
     double  lambda = 0;

     for(mcSampleIdx = 0; mcSampleIdx < m_numScan; mcSampleIdx ++) {
	  m_priorLLNew[mcSampleIdx] = singleSamplePriorLL(samplePtr,
							  m_neighborSumPtr,
							  mcSampleIdx,
							  m_beta,
							  m_numClusters,
							  m_numScan);
							  
	  m_priorLLOld[mcSampleIdx] = singleSamplePriorLL(samplePtr,
							  m_neighborSumPtr,
							  mcSampleIdx,
							  m_oldBeta,
							  m_numClusters,
							  m_numScan);

	  m_condLLNew[mcSampleIdx] = singleSampleCondLL(imagePtr,
							samplePtr,
							mcSampleIdx,
							m_cls);

	  m_condLLOld[mcSampleIdx] = singleSampleCondLL(imagePtr,
							samplePtr,
							mcSampleIdx,
							m_oldCls);

	  lambda = (m_priorLLNew[mcSampleIdx] + m_condLLNew[mcSampleIdx])
	       - (m_priorLLOld[mcSampleIdx] - m_condLLOld[mcSampleIdx]);
	       
	  if (m_verbose >= 2) {
	       printf("prioreLLNew = %f, priorLLOld = %f, condLLNew = %f, condLLOld = %f, lambda = %f\n", m_priorLLNew[mcSampleIdx], m_priorLLOld[mcSampleIdx], m_condLLNew[mcSampleIdx], m_condLLOld[mcSampleIdx], lambda);
	  }
     }

     return (0);
}

double singleSamplePriorLL(ImageType4DChar::Pointer samplePtr,
			   ImageType5DChar::Pointer neighborSumPtr,
			   unsigned short mcSampleIdx,
			   float beta,
			   unsigned short numClusters,
			   unsigned short numScan)
{

     double MZero = 0;
     double priorll = 0;
     int clsIdx = 0;
     unsigned char thisNeighborSum = 0;
     double expBetaSum = 0;

     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();

     // neighborSum.
     ImageType5DChar::IndexType neighborSumIdx;

     // define sample region. Only compute mcSampleIdx volume.
     sampleIdx.Fill(0);
     sampleIdx[3] = mcSampleIdx;

     sampleSize[3] = 1;
     ImageType4DChar::RegionType sampleRegion;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);
     
     IteratorType4DChar sampleIt(samplePtr, sampleRegion);



     for (sampleIt.GoToBegin(); !sampleIt.IsAtEnd(); ++ sampleIt) {
	  if (sampleIt.Get() >= 0) {
	       sampleIdx = sampleIt.GetIndex();
	       
	       neighborSumIdx[0] = sampleIdx[0];
	       neighborSumIdx[1] = sampleIdx[1];
	       neighborSumIdx[2] = sampleIdx[2];
	       neighborSumIdx[3] = sampleIdx[3];
	       
	       // Compute M0
	       MZero = 0;
	       for (neighborSumIdx[4] = 0; neighborSumIdx[4] < numClusters; neighborSumIdx[4]++) {

		    thisNeighborSum = neighborSumPtr->GetPixel(neighborSumIdx);
		    expBetaSum = exp(- beta * thisNeighborSum);
		    MZero += expBetaSum;
	       }
	       
	       neighborSumIdx[4] = sampleIt.Get();
	       thisNeighborSum = neighborSumPtr->GetPixel(neighborSumIdx);

	       priorll += (- beta * thisNeighborSum - log(MZero));

	  } // in mask.
     } // iterator

     return (priorll);
}

double singleSampleCondLL(ImageType4DFloat::Pointer imagePtr,
			  ImageType4DChar::Pointer samplePtr,
			  unsigned short mcSampleIdx,
			  std::vector<CompType> & cls)
			  

{
     using boost::math::cyl_bessel_i;
     double cll = 0;
     unsigned int clsIdx = 0;
     unsigned numClusters = cls.size();

     // images
     ImageType4DFloat::IndexType imageIdx;    
       ImageType4DFloat::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();

     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();

     sampleIdx.Fill(0);
     sampleIdx[3] = mcSampleIdx;
     sampleSize[3] = 1;

     ImageType4DChar::RegionType sampleRegion;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);

     IteratorType4DChar sampleIt(samplePtr, sampleRegion);

     vnl_vector <float> timeSeries(imageSize[3], 0);

     // compute normalization constant of von Mises
     // Fisher. vmfLogConst[i] is the log (c_d (kappa) ) for the i'th
     // clusters. See page 1350 of "Clustering on the Unit-Sphere
     // using von Mises Fisher..."
     vnl_vector<double> vmfLogConst(numClusters);
     float myD = imageSize[3];
     double const Pi = 4 * atan(1);

     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
//	  vmfLogConst[clsIdx] = (myD/2 - 1) * log (cls[clsIdx].kappa) - myD/2 * log(2*Pi) -  log(cyl_bessel_i(double(myD/2 - 1), double(cls[clsIdx].kappa)));

	  vmfLogConst[clsIdx] = (myD/2 - 1) * log (cls[clsIdx].kappa) - myD/2 * log(2*Pi) -  logBesselI(myD/2 - 1, cls[clsIdx].kappa);

//	  printf("SingleSampleCondLL(): clsIdx[%i]: cyl = %f, approx = %f\n", clsIdx, log(cyl_bessel_i(double(myD/2 - 1), double(cls[clsIdx].kappa))), logBesselI(myD/2 - 1, cls[clsIdx].kappa));

     }


     for (sampleIt.GoToBegin(); !sampleIt.IsAtEnd(); ++ sampleIt) {
	  if (sampleIt.Get() >= 0) {
	       sampleIdx = sampleIt.GetIndex();

	       // First save time series at this position into a vector.
	       imageIdx[0] = sampleIdx[0];
	       imageIdx[1] = sampleIdx[1];
	       imageIdx[2] = sampleIdx[2];	  
	       for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
		    timeSeries[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
	       }
	       
	       clsIdx = sampleIt.Get();
	       cll += vmfLogConst[clsIdx]
			 + cls[clsIdx].kappa * inner_product(timeSeries, cls[clsIdx].mu);
	  } // in mask.
     } // iterator
     return (cll);
}

int MCModel::estep(ImageType4DChar::Pointer samplePtr,
		   ImageType4DFloat::Pointer imagePtr,
		   unsigned int burnin)
{

     using boost::math::cyl_bessel_i;

     int scan = 0;
     double p_acpt = 0;
     double denergy = 0;
     int cand = 1;
     int currentLabel = 0;
     unsigned clsIdx = 0;
     double dlogll = 0;

     // images
     ImageType4DFloat::IndexType imageIdx;    
     ImageType4DFloat::RegionType imageRegion;
     imageRegion = imagePtr->GetLargestPossibleRegion();
     ImageType4DFloat::SizeType imageSize = imageRegion.GetSize();
     IteratorType4DFloat imageIt(imagePtr, imageRegion);

     // samples
     ImageType4DChar::IndexType sampleIdx;
     ImageType4DChar::RegionType sampleRegion;
     ImageType4DChar::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     sampleIdx.Fill(0);
     // only sample on last sample
     sampleIdx[3] = sampleSize[3] - 1;
     sampleSize[3] = 1;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);
     IteratorType4DChar sampleIt(samplePtr, sampleRegion);

     // Destination sample
     ImageType4DChar::RegionType sampleDestRegion;
     ImageType4DChar::IndexType sampleDestIdx;
     sampleDestIdx.Fill(0);
     sampleDestRegion.SetSize(sampleSize); // size is same.
     sampleDestRegion.SetIndex(sampleDestIdx);
//     IteratorType4DChar sampleDestIt(samplePtr, sampleDestRegion);



     vnl_vector <float> timeSeries(imageSize[3], 0);

     // Define neighborhood iterator
     typedef itk::ConstantBoundaryCondition< ImageType4DChar >  BoundaryConditionType;
     BoundaryConditionType constCondition;
     constCondition.SetConstant(-1);
     typedef itk::NeighborhoodIterator< ImageType4DChar, BoundaryConditionType > NeighborhoodIteratorType;
     
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);

     NeighborhoodIteratorType neighborSampleIt( radius, samplePtr, sampleRegion );
     neighborSampleIt.OverrideBoundaryCondition(&constCondition);

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1, 0}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1, 0}};


     ///////        Define random number generator.        ////////
     // Uniform integer generator.
     boost::uniform_int<> uni_int(0, m_numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     // compute normalization constant of von Mises
     // Fisher. vmfLogConst[i] is the log (c_d (kappa) ) for the i'th
     // clusters. See page 1350 of "Clustering on the Unit-Sphere
     // using von Mises Fisher..."
     vnl_vector<double> vmfLogConst(m_numClusters);
     float myD = imageSize[3];
     double const Pi = 4 * atan(1);

     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  vmfLogConst[clsIdx] = (myD/2 - 1) * log (m_cls[clsIdx].kappa) - myD/2 * log(2*Pi) -  logBesselI(myD/2 - 1, m_cls[clsIdx].kappa);
     }


     // sampling Markov Random Field. 
     scan = 0;
     for (scan = 0; scan < burnin + m_numScan; scan ++) {
	  for (neighborSampleIt.GoToBegin(); !neighborSampleIt.IsAtEnd(); ++ neighborSampleIt) {
	       currentLabel = neighborSampleIt.GetCenterPixel();
	       if (currentLabel >= 0) {

		    cand = roll_die();
		    denergy 
			 = int(cand != neighborSampleIt.GetPixel(xminus))
			 - int(currentLabel != neighborSampleIt.GetPixel(xminus))

			 + int(cand != neighborSampleIt.GetPixel(xplus)) 
			 - int(currentLabel != neighborSampleIt.GetPixel(xplus))

			 + int(cand != neighborSampleIt.GetPixel(yminus)) 
			 - int(currentLabel != neighborSampleIt.GetPixel(yminus))

			 + int(cand != neighborSampleIt.GetPixel(yplus)) 
			 - int(currentLabel != neighborSampleIt.GetPixel(yplus))

			 + int(cand != neighborSampleIt.GetPixel(zminus)) 
			 - int(currentLabel != neighborSampleIt.GetPixel(zminus))

			 + int(cand != neighborSampleIt.GetPixel(zplus)) 
			 - int(currentLabel != neighborSampleIt.GetPixel(zplus));

		    denergy = m_beta * denergy;

		    // likelihood energy.
		    sampleIdx = neighborSampleIt.GetIndex();
		    imageIdx[0] = sampleIdx[0];
		    imageIdx[1] = sampleIdx[1];
		    imageIdx[2] = sampleIdx[2];

		    // First save time series at this position into a vector.
		    for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
			 timeSeries[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
		    }


		    // The energy = - Log(likelihood).
		    dlogll = (- vmfLogConst[cand] - m_cls[cand].kappa * inner_product(timeSeries, m_cls[cand].mu)) - (- vmfLogConst[currentLabel] - m_cls[currentLabel].kappa * inner_product(timeSeries, m_cls[currentLabel].mu));
		    dlogll = dlogll / imageSize[3];
		    dlogll = dlogll * 5;
		    denergy = denergy + dlogll;

		    // if energy change less than zero, just accept
		    // candidate. otherwise accept with exp(- energy
		    // change).
			 
		    if (denergy <= 0) {
			 neighborSampleIt.SetCenterPixel(cand);
		    }
		    else {
			 p_acpt = exp(-denergy);
			 if (uni() < p_acpt) {
			      neighborSampleIt.SetCenterPixel(cand);
			 }
		    }
	       } // in mask
	  } // iterators.
	  if (m_verbose >= 1) {
	       printf("nChainSamping(): scan %u done.\n", scan);
	  }

	  // Save it to correct place.
	  if (scan >= burnin) {
	       sampleDestIdx.Fill(0);
	       sampleDestIdx[3] = scan - burnin;

	       sampleDestRegion.SetSize(sampleSize); // size is same.
	       sampleDestRegion.SetIndex(sampleDestIdx);
	       IteratorType4DChar sampleDestIt(samplePtr, sampleDestRegion);

	       for (sampleIt.GoToBegin(), sampleDestIt.GoToBegin(); !sampleIt.IsAtEnd(); ++ sampleIt, ++ sampleDestIt) {
		    sampleDestIt.Set(sampleIt.Get());
	       }

	  }
     } // for scan
	  

	  
     // Since labels have changed. we need update whatever that
     // depends on labels, including the neighborSum data structure.

     // Compute neighborSum.
     evalNeighborSum(samplePtr);
     // compute prior LL and cond LL. In fact only prior LL
     // changed. Anyway...
     updateSampleLL(imagePtr, samplePtr);

     // DBG
     saveimage4dchar(samplePtr, "samples.nii");
     saveExtractVolume(samplePtr, sampleSize[3] - 1, "onesample.nii");

     // end of DBG
     return 0;    
}

int MCModel::printMeanKappa()
{
     unsigned int timeSeriesLength = m_cls[0].mu.size();
     unsigned clsIdx = 0;
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  std::cout << "POST: " << m_cls[clsIdx].kappa << " ";
	  printVnlVector(m_cls[clsIdx].mu, timeSeriesLength);     
     }
}
