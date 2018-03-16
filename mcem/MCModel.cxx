#include "commonalt.h"
#include "MCModel.h"
#include "GCoptimization.h"

#include <boost/math/distributions/students_t.hpp>
twister_base_gen_type mygenerator(42u);

MCModel::MCModel(ImageType2D::Pointer imagePtr,
		 ImageType3D::Pointer samplePtr,
		 unsigned int numClusters = 2,
		 unsigned int numScan = 10,
		 double beta = 0.7,
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
     
     m_beta = beta;
     priorll_ = 0;
     m_cll = 0;
     m_oldPriorll = 0;
     m_cll = 0;
     m_verbose = verbose;
     

     ImageType2D::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     m_numScan = numScan;
     ImageType2D::IndexType imageIdx;     
     ImageType3D::IndexType sampleIdx;     
     const unsigned numRepeat = 10;


//     mygenerator.seed(static_cast<unsigned int>(std::time(0)));
     mygenerator.seed(static_cast<unsigned int>(2011));

     // Uniform integer generator.
     boost::uniform_int<> uni_int(1, m_numClusters); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     m_cls  = new CompType[m_numClusters];
     m_oldCls  = new CompType[m_numClusters];

     unsigned repeatIdx = 0;
     unsigned clsIdx = 0;
     double minMeanSSR = 1e10, meanSSR = 0;

     // we define a bestCls to save the best mu so far.
     CompType * bestCls  = new CompType[m_numClusters];

     // Init Gaussian Mixture parameters mu, numPoints and sigma
     // Repeat kmeans clustering a few times and choose the best one
     // by sum of square error.
     for (repeatIdx = 0; repeatIdx < numRepeat; repeatIdx ++) {
	  printf("kmeans iteration %i begin:\n", repeatIdx);
	  meanSSR = kmeans(imagePtr, samplePtr);
	  if (meanSSR < minMeanSSR) {
	       minMeanSSR = meanSSR;
	       // Save best mu in to bestCls
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
		    bestCls[clsIdx].mu = m_cls[clsIdx].mu;
		    bestCls[clsIdx].numPoints = m_cls[clsIdx].numPoints;
	       }

	       if (m_verbose >=1) {
		    printf("initClusterCenter: meanSSR %f < minMeanSSR . best cluster centers are:\n", meanSSR);
		    print("normal");
	       }
	  }
     }

     // Found the best mu. Restore the best mu and numPoints into cls.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
	  m_cls[clsIdx].mu = bestCls[clsIdx].mu;
	  m_cls[clsIdx].numPoints = bestCls[clsIdx].numPoints;
     }

     // Sort labels such that small labels have small mean.
     sortLabels(samplePtr);    

     // Given mu, estimate labels. Save to last slice of samplePtr.
     estLabelsByMean(imagePtr, samplePtr);     
     
     // Based on the mu, numPoints and labels, estimate sigma. First we clear
     // all sigma to zero.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_cls[clsIdx].sigma = 0;
     }
     
     // Compute sigma_k.
     sampleIdx[2] = sampleSize[2] - 1;
     for (sampleIdx[0] = 0;sampleIdx[0] < sampleSize[0]; sampleIdx[0] ++) {
	  for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1] ++) {
		    imageIdx[0] = sampleIdx[0];
		    imageIdx[1] = sampleIdx[1];
		    clsIdx = samplePtr->GetPixel(sampleIdx);
		    m_cls[clsIdx].sigma += pow((imagePtr->GetPixel(imageIdx) - m_cls[clsIdx].mu), 2);
	       } 
	  } 
     
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_cls[clsIdx].sigma = sqrt(m_cls[clsIdx].sigma/m_cls[clsIdx].numPoints);
     }


     // MLE: estimate labels and save in last slice of samplePtr.
     estLabelML(samplePtr, imagePtr);

     // Manually set beta to a small value. This is not different with
     // estimating beta using labels from MLE above. Anyway when EM
     // starts, beta will be estimated from MC (or ICM) samples. Only
     // difference is one the first iteration of E step sampling, we
     // ignore (because very small beta) spatial smoothness prior, and
     // sample from posterior is equivalent to sampling from ML.
     m_beta = beta;

     // With known parameters, Use ICM to estimate labels and
     // save it in last slice of samplePtr.
     icmEstep(imagePtr, samplePtr, "lastone", 20);

     // Init N indepedent Markov chains with same results from MLE or ICM.
     if (samplingMethod.compare("nchain") == 0) {
	  nChainInit(samplePtr, "ml");
     }

     saveExtractSlice(samplePtr, sampleSize[2]-1, "init.nii");

     // init m_singleMCSampleLLNew
     m_singleMCSampleLLNew.set_size(m_numScan);
     m_singleMCSampleLLOld.set_size(m_numScan);

     delete [] bestCls;
}
     
MCModel::~MCModel()
{
     delete [] m_cls;
     delete [] m_oldCls;
}

int MCModel::print(std::string format)
{
     unsigned short clsIdx;
     if (format.compare("normal") == 0) {
	       for (int clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
		    printf("cluster[%2i]: # of points: %ld,   mu = %6.2f, sigma = %4.2f.\n", clsIdx, m_cls[clsIdx].numPoints, m_cls[clsIdx].mu, m_cls[clsIdx].sigma); 
	       }
	       printf("beta = %8f,   prior ll = %10.2f,   cond ll = %10.2f,   joint ll = %10.2f, old joint ll = %10.2f\n", m_beta, priorll_, m_cll, priorll_ + m_cll, m_oldPriorll + m_oldCll);
	  }
     else if (format.compare("table") == 0) {
	  printf("%2.4f %.2f %.2f %.2f   ", m_beta, priorll_, m_cll, priorll_ + m_cll);
	  for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	       printf("%6.4f %6.4f %10i  ", m_cls[clsIdx].mu, m_cls[clsIdx].sigma, m_cls[clsIdx].numPoints);
	  }
	  printf("\n");
	  printf("NEWTAG %i %f %f \n", m_numScan, m_beta, priorll_ + m_cll);
	       
	  }
     else {
	  printf("MCModel::print(): print format not recognized.\n");
     }
     
     return 0;
}

int MCModel::estep(ImageType2D::Pointer imagePtr,
		   ImageType3D::Pointer samplePtr,
		   unsigned int burnin)
{
     int scan = 0;
     double p_acpt = 0;
     double denergy = 0;
     int cand = 1;
     int label = 0;
     int nLabel = 0; // neighobr label
     float intensity = 0; // current pixel intensity.

     ImageType3D::IndexType nidx;    
     ImageType2D::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::SizeType sampleSize =  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;
     ImageType3D::IndexType sampleDestIdx;     
     ImageType2D::IndexType imageIdx;


     //////////////////////////////////////////////////////////////
     ///////        Define random number generator.        ////////
     //////////////////////////////////////////////////////////////

     // Uniform integer generator.
     boost::uniform_int<> uni_int(1, m_numClusters); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);


     // sampling Markov Random Field. 
     scan = 0;
     sampleIdx[2] = sampleSize[2] - 1;
     nidx[2] = sampleSize[2] - 1;

     while (scan < burnin + m_numScan) {
     	  for (sampleIdx[0] = 0; sampleIdx[0] < sampleSize[0]; sampleIdx[0]++) {
     	       for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1]++) {
		    cand = roll_die() - 1;
		    label = samplePtr->GetPixel(sampleIdx);
     		    denergy = 0;
     		    if (sampleIdx[0] > 0) {
			 nidx[0] = sampleIdx[0] - 1;
			 nidx[1] = sampleIdx[1];
			 nLabel = samplePtr->GetPixel(nidx); // neighbor label
			 denergy = denergy + int(cand != nLabel) - int(label != nLabel);
		    }

     		    if (sampleIdx[0] < sampleSize[0]-1) {
			 nidx[0] = sampleIdx[0] + 1;
			 nidx[1] = sampleIdx[1];
			 nLabel = samplePtr->GetPixel(nidx); // neighbor label
			 denergy = denergy + int(cand != nLabel) - int(label != nLabel);
		    }

     		    if (sampleIdx[1] > 0) {
			 nidx[0] = sampleIdx[0];
			 nidx[1] = sampleIdx[1] - 1;
			 nLabel = samplePtr->GetPixel(nidx); // neighbor label
			 denergy = denergy + int(cand != nLabel) - int(label != nLabel);
		    }

     		    if (sampleIdx[1] < sampleSize[1]-1) {
			 nidx[0] = sampleIdx[0];
			 nidx[1] = sampleIdx[1] + 1;
			 nLabel = samplePtr->GetPixel(nidx); // neighbor label
			 denergy = denergy + int(cand != nLabel) - int(label != nLabel);
		    }

     		    denergy = m_beta * denergy;

		    // likelihood energy.
		    imageIdx[0] = sampleIdx[0];
		    imageIdx[1] = sampleIdx[1];
		    intensity = imagePtr->GetPixel(imageIdx);
		    denergy = denergy
			 + (intensity - m_cls[cand].mu)*(intensity - m_cls[cand].mu)/(2*m_cls[cand].sigma * m_cls[cand].sigma) + log(m_cls[cand].sigma) // Candidate's likelihood energy.
			 - (intensity - m_cls[label].mu)*(intensity - m_cls[label].mu)/(2*m_cls[label].sigma * m_cls[label].sigma) - log(m_cls[label].sigma); // Current label's likelihood energy.
		    
     		    // if energy change less than zero, just accept
     		    // candidate. otherwise accept with exp(- energy
     		    // change).

     		    if (denergy <= 0) {
			 
			 samplePtr->SetPixel(sampleIdx, cand);
     		    }
     		    else {
     			 p_acpt = exp(-denergy);
     			 if (uni() < p_acpt) {
			      samplePtr->SetPixel(sampleIdx, cand);
     			 }
     		    }
     	       }
     	  }
	  
	  if (scan >= burnin) {
	       // Save last slice. 
	       sampleDestIdx[2] = scan - burnin;
	       
	       for (sampleIdx[0] = 0; sampleIdx[0] < sampleSize[0]; sampleIdx[0] ++) {
		    for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1]++) {
			 sampleDestIdx[0] = sampleIdx[0];
			 sampleDestIdx[1] = sampleIdx[1];
			 samplePtr->SetPixel(sampleDestIdx, samplePtr->GetPixel(sampleIdx));
		    }
	       }
	  }

	  if (scan <= 1) {
	       printf("Scan %i...", ++scan);
	  }
	  else {
	       printf("%i...", ++scan);
	  }
     }
     printf("\n");


     m_oldCll = m_cll;
     m_cll = eval_cll(imagePtr, samplePtr);

     m_oldPriorll = priorll_;
     priorll_ = eval_ll(samplePtr, m_beta);

     return 0;    
}

// Given labels, estimate mu and numPoints.
int MCModel::estimateMu(ImageType3D::Pointer samplePtr,
			ImageType2D::Pointer imagePtr)
{
     int k = 0; // cluster number.

     ImageType3D::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType2D::IndexType idx;     
     ImageType3D::IndexType sampleIdx;     

     // Compute number of points in each cluster. Save old parameters in m_oldCls;
     for (k = 0; k < m_numClusters; k ++) {
	  m_oldCls[k].numPoints = m_cls[k].numPoints;
	  m_oldCls[k].mu = m_cls[k].mu;
	  m_cls[k].numPoints = 0;
	  m_cls[k].mu = 0;
     }
     
     for (sampleIdx[0] = 0; sampleIdx[0] < sampleSize[0]; sampleIdx[0] ++) {
	  for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1] ++) {
	       for (sampleIdx[2] = 0; sampleIdx[2] < m_numScan; sampleIdx[2] ++) {
		    idx[0] = sampleIdx[0];
		    idx[1] = sampleIdx[1];
		    k = samplePtr->GetPixel(sampleIdx);
		    m_cls[k].mu += imagePtr->GetPixel(idx);
		    m_cls[k].numPoints ++;
	       } 
	  } 
     }
     // compute mu_k
     for (k = 0; k < m_numClusters; k ++) {
	  m_cls[k].mu = m_cls[k].mu/m_cls[k].numPoints;
     }
     
     return 0;
}

// Given labels, mu and numPoints, estimate sigma.
int MCModel::estimateSigma(ImageType3D::Pointer samplePtr,
			   ImageType2D::Pointer imagePtr)
{
     unsigned int k = 0; // cluster number.

     ImageType3D::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType2D::IndexType idx;     
     ImageType3D::IndexType sampleIdx;     

     // Compute number of points in each cluster. Save old parameters in m_oldCls;
     for (k = 0; k < m_numClusters; k ++) {
	  m_oldCls[k].sigma = m_cls[k].sigma;
	  m_cls[k].sigma = 0;
     }
     
     // Compute sigma_k.
     for (sampleIdx[0] = 0;sampleIdx[0] < sampleSize[0]; sampleIdx[0] ++) {
	  for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1] ++) {
	       for (sampleIdx[2] = 0; sampleIdx[2] < m_numScan; sampleIdx[2] ++) {
		    idx[0] = sampleIdx[0];
		    idx[1] = sampleIdx[1];
		    k = samplePtr->GetPixel(sampleIdx);
		    m_cls[k].sigma += (imagePtr->GetPixel(idx) - m_cls[k].mu) * (imagePtr->GetPixel(idx) - m_cls[k].mu);
	       } 
	  } 
     }
     
     for (k = 0; k < m_numClusters; k ++) {
	  m_cls[k].sigma = sqrt(m_cls[k].sigma/m_cls[k].numPoints);
     }

     // update conditional log-likelihood P(d | f).
     m_oldCll = m_cll;
     m_cll = eval_cll(imagePtr, samplePtr);

     return 0;
}


// Given a beta, return the prior log-likelihood.
double MCModel::eval_ll(ImageType3D::Pointer samplePtr,
		     double curBeta)
{
     int k = 0; // cluster number.
     int S = 0, curS = 0; //  neighbors sum.
     double sumexp = 0;
     double Q = 0; // log likelihood (prior).
     int thisLabel = 0;
     int curLabel = 0;

     ImageType3D::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;     
     ImageType3D::IndexType nidx;     

     Q = 0;

     for (sampleIdx[2] = 0; sampleIdx[2] < m_numScan; sampleIdx[2] ++) {
	  nidx[2] = sampleIdx[2];
	  for (sampleIdx[0] = 0; sampleIdx[0] < sampleSize[0]; sampleIdx[0] ++) {
	       for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1] ++) {
		    curLabel = samplePtr->GetPixel(sampleIdx);
		    sumexp = 0;
		    for (thisLabel = 0; thisLabel < m_numClusters; thisLabel++) {
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
			 sumexp = sumexp + exp(-curBeta * S);
			 if (thisLabel == curLabel) {
			      curS = S;
			 }
		    }
			 
		    Q += (-curBeta * curS - log(sumexp))/double(m_numScan);
	       }
	  }
     }

     return Q;
}



// Maximization prior log-likelihood over beta. return beta.
int  MCModel::estimatePriorParameter(double a, double b,
				     ImageType3D::Pointer samplePtr)

{
     if (a >= b) {
	  printf("estimatePriorParameter: a must be smaller than b.\n");
	  exit(1);
     }
     double tau = (sqrt(5) - 1)/2;
     double x1 = 0, x2 = 0;
     double f1 = 0, f2 = 0;
     x1 = a + (1 - tau) * (b - a);
     f1 =  - eval_ll(samplePtr, x1);
     x2 = a + tau * (b - a);

     // compute the minimal value of negagive log-likelihood.
     f2 =  - eval_ll(samplePtr, x2);
     while (fabs(b - a) > 0.001) {
	  if (f1 > f2) {
	       a = x1;
	       x1 = x2;
	       f1 = f2;
	       x2 = a + tau * (b - a);
	       f2 =  - eval_ll(samplePtr, x2);
	  }
	  else {
	       b = x2;
	       x2 = x1;
	       f2 = f1;
	       x1 = a + (1 - tau) * (b - a);
	       f1 =  - eval_ll(samplePtr, x1);
	  }
	  printf("a = %f, x1 = %f, x2 = %f, b = %f, f1 = %f, f2 = %f\n", a, x1, x2, b, f1, f2);
     }
     // Save old beta in m_oldBeta.
     m_oldBeta = m_beta;
     m_beta = (a + b)/2;

     // Update prior log-likelihood.
     m_oldPriorll = priorll_;
     priorll_ = eval_ll(samplePtr, m_beta);
     
     return 0;
}


// Given mu and sigma, evaluate the conditional log-likelihood P(d | f).
double MCModel::eval_cll(ImageType2D::Pointer imagePtr,
			 ImageType3D::Pointer samplePtr)
{
     ImageType3D::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;     
     ImageType2D::IndexType idx;     
     double cll = 0;
     int k = 0;
     
     for (sampleIdx[2] = 0; sampleIdx[2] < m_numScan; sampleIdx[2] ++) {
	  for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1] ++) {
	       for (sampleIdx[0] = 0; sampleIdx[0] < sampleSize[0]; sampleIdx[0] ++) {
		    idx[0] = sampleIdx[0];
		    idx[1] = sampleIdx[1];
		    k = samplePtr->GetPixel(sampleIdx);
		    cll = cll + (1/double(m_numScan)) * (
			 pow((imagePtr->GetPixel(idx) - m_cls[k].mu), 2)/(-2*pow(m_cls[k].sigma, 2))
			 - log(m_cls[k].sigma)
			 );
	       } 
	  } 
     }
     return cll;
}


// Test the model convergence. Will return 1 if iteration converge, otherwise return 0.
int MCModel::testConvergence(float maxRelativeErr)
{
     double relerr = 0;
     double thisMaxErr = 0;
     double thisRelErrMu = 0;
     double thisRelErrSigma = 0;
     for (int k = 0; k < m_numClusters; k++) {
//	  printf("m_cls[%i].mu = %f, OldM_Cls[%i].mu = %f, relerr.mu = %f, m_cls[%i].sigma = %f, OldM_Cls[%i].sigma = %f, relerr.sigma = %f\n", k, m_cls[k].mu, k, m_oldCls[k].mu, fabs((m_cls[k].mu - m_oldCls[k].mu)/m_cls[k].mu), k, m_cls[k].sigma, k, m_oldCls[k].sigma, fabs((m_cls[k].sigma - m_oldCls[k].sigma)/m_cls[k].sigma));

	  if (fabs((m_cls[k].mu - m_oldCls[k].mu)/m_cls[k].mu) > maxRelativeErr ||
	      fabs((m_cls[k].sigma - m_oldCls[k].sigma)/m_cls[k].sigma) > maxRelativeErr) {
	       return 0;
	  }
	       
     }
     if (fabs((m_beta - m_oldBeta)/m_beta) > maxRelativeErr) {
	  return 0;
     }

     return 1;
}


// update prior LL and conditional LL.
int MCModel::updateLL(ImageType2D::Pointer imagePtr,
		      ImageType3D::Pointer samplePtr)
{
     priorll_ = eval_ll(samplePtr, m_beta);
     m_cll = eval_cll(imagePtr, samplePtr);
}

int MCModel::ll_over_beta(ImageType3D::Pointer samplePtr)

{
     for (double beta = 0.001; beta < 5; beta = beta + 0.2) {
	  
	  printf("%f %f\n", beta, eval_ll(samplePtr, beta));
     }
     return 0;
}

// Init cluster center by K-means.
double MCModel::kmeans(ImageType2D::Pointer imagePtr, ImageType3D::Pointer samplePtr)
{
     unsigned int clsIdx = 0;
     ImageType2D::IndexType imageIdx;    
     ImageType2D::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;
     ImageType3D::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);
	  
     // Randomly init the cluster center for K-means clustering.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
	  imageIdx[0] = floor(uni() * imageSize[0]);
	  imageIdx[1] = floor(uni() * imageSize[1]);
	  m_cls[clsIdx].mu = imagePtr->GetPixel(imageIdx);
     }

     sampleIdx[2] = sampleSize[2] - 1;
     float pixelIntensity = 0;
     double meanSSR = 0;
     do {
	  // update label assignment.
	  meanSSR = estLabelsByMean(imagePtr, samplePtr);

	  // update cluster center.
	  // reset mu and numPoints. But save it before erasing it.
	  for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
	       m_oldCls[clsIdx].mu = m_cls[clsIdx].mu;
	       m_cls[clsIdx].mu = 0;
	       m_cls[clsIdx].numPoints = 0;
	  }

	  for (imageIdx[0] = 0; imageIdx[0] < imageSize[0]; imageIdx[0] ++) {
	       for (imageIdx[1] = 0; imageIdx[1] < imageSize[1]; imageIdx[1] ++) {

		    sampleIdx[0] = imageIdx[0];
		    sampleIdx[1] = imageIdx[1];
		    clsIdx = samplePtr->GetPixel(sampleIdx);
		    m_cls[clsIdx].mu = m_cls[clsIdx].mu * (double(m_cls[clsIdx].numPoints)/double(m_cls[clsIdx].numPoints + 1)) + imagePtr->GetPixel(imageIdx)/double(m_cls[clsIdx].numPoints + 1);
		    m_cls[clsIdx].numPoints ++;
	       }
	  }	       
     }
     while(testClusterCenterMove());
     
     return meanSSR;
}

double MCModel::estLabelsByMean(ImageType2D::Pointer imagePtr,
			     ImageType3D::Pointer samplePtr)
{
     unsigned int clsIdx = 0;
     ImageType2D::IndexType imageIdx;    
     ImageType2D::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;
     ImageType3D::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();

     unsigned nearestClsLabel = 0;
     double nearestDistance = 1000000;
     double meanSSR = 0;

     sampleIdx[2] = sampleSize[2] - 1;
     float pixelIntensity = 0;

     for (imageIdx[0] = 0; imageIdx[0] < imageSize[0]; imageIdx[0] ++) {
	  for (imageIdx[1] = 0; imageIdx[1] < imageSize[1]; imageIdx[1] ++) {
	       nearestDistance = 1000000;
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
		    pixelIntensity = imagePtr->GetPixel(imageIdx);
		    if (pow(pixelIntensity - m_cls[clsIdx].mu, 2) < nearestDistance) {
			 nearestClsLabel = clsIdx;
			 nearestDistance = pow(pixelIntensity - m_cls[clsIdx].mu, 2);
		    }
	       }
	       // Found the nearest cluster label.
	       sampleIdx[0] = imageIdx[0];
	       sampleIdx[1] = imageIdx[1];
	       samplePtr->SetPixel(sampleIdx, nearestClsLabel);
	       meanSSR += nearestDistance;
	  }
     }
     meanSSR = meanSSR / double(imageSize[0] * imageSize[1]);
     return meanSSR;
}

int MCModel::testClusterCenterMove()
{
     unsigned int clsIdx = 0;
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
	  if (m_cls[clsIdx].mu != m_oldCls[clsIdx].mu) {
	       return 1;
	  }
     }
     return 0;
}

// Compute the qutient of the first and second derivative.
float MCModel::evalDerivativeQutient(ImageType3D::Pointer samplePtr,
				     Image4DChar::Pointer neighborSumPtr,
				     float * firstDeriv,
				     float * secondDeriv)
{
     Image4DChar::RegionType neiborSumRegion = neighborSumPtr->GetLargestPossibleRegion();
     ImageType3D::IndexType sampleIdx;     
//     double secondDeriv = 0;
//     double firstDeriv = 0;
     double MZero = 0, MOne = 0, MTwo = 0;

     double expBetaSum = 0;
     Image4DChar::SizeType neighborSumSize = 
	  neighborSumPtr->GetLargestPossibleRegion().GetSize();
     Image4DChar::IndexType neighborSumIdx;
     unsigned int clsIdx = 0;
     unsigned char thisNeighborSum = 0;

     * firstDeriv = 0;
     * secondDeriv = 0;
     for (neighborSumIdx[0] = 0; neighborSumIdx[0] < neighborSumSize[0]; neighborSumIdx[0]++) {
	  for (neighborSumIdx[1] = 0; neighborSumIdx[1] < neighborSumSize[1]; neighborSumIdx[1]++) {
	       for (neighborSumIdx[2] = 0; neighborSumIdx[2] < neighborSumSize[2]; neighborSumIdx[2]++) {
		    // Compute M0, M1, and M2
		    MZero = 0, MOne = 0, MTwo = 0;
		    for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
			 neighborSumIdx[3] = clsIdx;
			 thisNeighborSum = neighborSumPtr->GetPixel(neighborSumIdx);
			 expBetaSum = exp(- m_beta * thisNeighborSum);
			 MZero += expBetaSum;
			 expBetaSum = expBetaSum * thisNeighborSum; // Now it's delta M1.
			 MOne += expBetaSum;
			 expBetaSum = expBetaSum * thisNeighborSum; // Now it's delta M2.
			 MTwo += expBetaSum;

		    }

		    sampleIdx[0] = neighborSumIdx[0];
		    sampleIdx[1] = neighborSumIdx[1];
		    sampleIdx[2] = neighborSumIdx[2];
		    neighborSumIdx[3] = samplePtr->GetPixel(sampleIdx);
		    thisNeighborSum = neighborSumPtr->GetPixel(neighborSumIdx);
		    *firstDeriv += (double)(-thisNeighborSum + MOne/MZero)/ (double)(neighborSumSize[2]);
		    *secondDeriv += (double)(MOne*MOne - MZero * MTwo)/(double)(MZero * MZero * neighborSumSize[2]);

	       }
	  }
     }
     return (*firstDeriv) / (*secondDeriv);
}

int MCModel::evalNeighborSum(ImageType3D::Pointer samplePtr,
			     Image4DChar::Pointer neighborSumPtr)
{
     ImageType3D::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     sampleSize[2] = m_numScan;
     ImageType3D::IndexType sampleIdx;     
     sampleIdx[0] = 0;
     sampleIdx[1] = 0;
     sampleIdx[2] = 0;
     ImageType3D::RegionType sampleRegion;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);

     ConstIteratorType3D sampleIt(samplePtr, sampleRegion);

     Image4DChar::IndexType neighborSumIdx;
     unsigned int neighborSum = 0;
     
     sampleIt.GoToBegin();
     unsigned thisLabel = 0;
     while(!sampleIt.IsAtEnd()) {
	  sampleIdx = sampleIt.GetIndex();
	  neighborSumIdx[0] = sampleIdx[0];
	  neighborSumIdx[1] = sampleIdx[1];	  
	  neighborSumIdx[2] = sampleIdx[2];

	  for (thisLabel = 0; thisLabel < m_numClusters; thisLabel++) {

	       neighborSumIdx[3] = thisLabel;
	       neighborSum = 0;

	       sampleIdx = sampleIt.GetIndex();			 
	       if (sampleIdx[0] > 0) {
		    sampleIdx[0] --;
		    neighborSum += int(thisLabel != samplePtr->GetPixel(sampleIdx));
	       }

	       sampleIdx = sampleIt.GetIndex();			 
	       if  (sampleIdx[0] < sampleSize[0] - 1) {
		    sampleIdx[0] ++;
		    neighborSum += int(thisLabel != samplePtr->GetPixel(sampleIdx));
	       }

	       sampleIdx = sampleIt.GetIndex();
	       if  (sampleIdx[1] > 0) {
		    sampleIdx[1] --;
		    neighborSum += int(thisLabel != samplePtr->GetPixel(sampleIdx));
	       }
	       
	       sampleIdx = sampleIt.GetIndex();
	       if  (sampleIdx[1] < sampleSize[1] - 1) {
		    sampleIdx[1] ++;
		    neighborSum += int(thisLabel != samplePtr->GetPixel(sampleIdx));
	       }
	  
	       neighborSumPtr->SetPixel(neighborSumIdx, neighborSum);
	  }
	  ++sampleIt;
     }
     return 0;
}

int MCModel::estPriorPramNewton(ImageType3D::Pointer samplePtr, unsigned maxNewtonIter)
{
     ImageType3D::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();

     // Allocate memory for neighboring sum variable S.
     Image4DChar::Pointer neighborSumPtr = Image4DChar::New();
     Image4DChar::SizeType neighborSumSize;
     neighborSumSize[0] = sampleSize[0];
     neighborSumSize[1] = sampleSize[1];
     // only use current num of scans.
     neighborSumSize[2] = m_numScan; 
     neighborSumSize[3] = m_numClusters;
     Image4DChar::IndexType neighborSumIdx;
     neighborSumIdx[0] = 0;
     neighborSumIdx[1] = 0;
     neighborSumIdx[2] = 0;
     neighborSumIdx[3] = 0;
     Image4DChar::RegionType region4d;
     region4d.SetSize(neighborSumSize);
     region4d.SetIndex(neighborSumIdx);
     neighborSumPtr->SetRegions(region4d);
     neighborSumPtr->Allocate();

     // Compute S.
     evalNeighborSum(samplePtr, neighborSumPtr);
     
     // Newton's method.
     float stepLength = 1000;
     float firstDeriv = 0, secondDeriv = 0;
     unsigned newtonIter = 0;
     m_oldBeta = m_beta;
     do 
     {
	  evalDerivativeQutient(samplePtr, neighborSumPtr, &firstDeriv, &secondDeriv);
	  stepLength = firstDeriv/secondDeriv;
	  
	  if (m_verbose >= 2) {
	       printf("estPriorPramNewton(): current beta = %f, stepLength = %f, 1st deriv = %f, 2nd deriv = %f.\n", m_beta, stepLength, firstDeriv, secondDeriv);
	  }
	  m_beta = m_beta - vnl_math_sgn(stepLength) * vnl_math_min(vnl_math_abs(stepLength), 0.1);

	  newtonIter++;
     }
     while(fabs(stepLength) > 0.0001 && newtonIter < maxNewtonIter);
     m_oldPriorll = priorll_;
     priorll_ = eval_ll(samplePtr, m_beta);

//     testDerivatives(samplePtr, neighborSumPtr);     
     
}


int MCModel::testDerivatives(ImageType3D::Pointer samplePtr,
			     Image4DChar::Pointer neighborSumPtr)

{
     float firstDer = 0;
     float secondDer = 0;
     double betaBackup = m_beta;
     
     printf("beta firstDer secondDer\n");
     for (m_beta = 0.01; m_beta < 6; m_beta = m_beta + 0.01) {
	  evalDerivativeQutient(samplePtr, neighborSumPtr, &firstDer, &secondDer);
	  printf("beta = %f, f(beta) = %f, 1st = %f, 2nd = %f\n", m_beta, eval_ll(samplePtr, m_beta), firstDer, secondDer);
     }
     m_beta = betaBackup;
}


int MCModel::icmEstep(ImageType2D::Pointer imagePtr,
		      ImageType3D::Pointer samplePtr,
		      std::string scope,
		      unsigned int maxScan)

{
     unsigned int scan = 0;
     double postEnergy = 0;
     double smallestPostEnergy = 10000;
     unsigned char label = 0;
     unsigned char nLabel = 0; // neighobr label
     unsigned clsIdx = 0;
     unsigned char bestLabel = 0;
     float intensity = 0; // current pixel intensity.

     unsigned short beginSliceIdx;
     unsigned short endSliceIdx;


     ImageType3D::IndexType sampleIdx;
     ImageType2D::IndexType imageIdx;
     ImageType3D::IndexType nidx;    
     ImageType2D::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();

     if (scope.compare("lastone") == 0) {
	  beginSliceIdx = sampleSize[2] - 1;
	  endSliceIdx = sampleSize[2] - 1;
     }
     else if (scope.compare("all") == 0) {
	  beginSliceIdx = 0;
	  endSliceIdx = m_numScan - 1;
     }
     else {
	  printf("icmEstep(): scope must be either all or lastone. exit.\n");
	  exit(1);
     }

     ///////        Define random number generator.        ////////
     // Uniform integer generator.
     boost::uniform_int<> uni_int(1, m_numClusters); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);
     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     // sampling Markov Random Field. 

     unsigned anyChangePixel = 0;
     for (sampleIdx[2] = beginSliceIdx; sampleIdx[2] <= endSliceIdx; sampleIdx[2] ++) {
	  scan = 0;
	  nidx[2] = sampleIdx[2];
	  do {
	       anyChangePixel = 0;
	       for (sampleIdx[0] = 0; sampleIdx[0] < imageSize[0]; sampleIdx[0]++) {
		    for (sampleIdx[1] = 0; sampleIdx[1] < imageSize[1]; sampleIdx[1]++) {
			 smallestPostEnergy = 10000;
			 for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
			      postEnergy = 0;
			      if (sampleIdx[0] > 0) {
				   nidx[0] = sampleIdx[0] - 1;
				   nidx[1] = sampleIdx[1];
				   nLabel = samplePtr->GetPixel(nidx); // neighbor label
				   postEnergy += int(clsIdx != nLabel);
			      }

			      if (sampleIdx[0] < imageSize[0]-1) {
				   nidx[0] = sampleIdx[0] + 1;
				   nidx[1] = sampleIdx[1];
				   nLabel = samplePtr->GetPixel(nidx); // neighbor label
				   postEnergy += int(clsIdx != nLabel);
			      }

			      if (sampleIdx[1] > 0) {
				   nidx[0] = sampleIdx[0];
				   nidx[1] = sampleIdx[1] - 1;
				   nLabel = samplePtr->GetPixel(nidx); // neighbor label
				   postEnergy += int(clsIdx != nLabel);
			      }

			      if (sampleIdx[1] < imageSize[1]-1) {
				   nidx[0] = sampleIdx[0];
				   nidx[1] = sampleIdx[1] + 1;
				   nLabel = samplePtr->GetPixel(nidx); // neighbor label
				   postEnergy += int(clsIdx != nLabel);
			      }

			      postEnergy = m_beta * postEnergy;

			      // likelihood energy.
			      imageIdx[0] = sampleIdx[0];
			      imageIdx[1] = sampleIdx[1];
			      intensity = imagePtr->GetPixel(imageIdx);
			      postEnergy += (intensity - m_cls[clsIdx].mu)*(intensity - m_cls[clsIdx].mu)/(2*m_cls[clsIdx].sigma * m_cls[clsIdx].sigma) + log(m_cls[clsIdx].sigma); // Candidate's likelihood energy.

			      // Now we get posterior energy. We want
			      // to find the smallest one among all
			      // labels from 1 to numClusters.
			      if (postEnergy < smallestPostEnergy) {
				   bestLabel = clsIdx;
				   smallestPostEnergy = postEnergy;
			      }
			 }
		    
			 // By ICM, we just assign the bestLabel to current pixel.
			 if (bestLabel != samplePtr->GetPixel(sampleIdx)) {
			      if (m_verbose >= 3) {
			      printf("[%i,%i, %i->%i] ", sampleIdx[0], sampleIdx[1], samplePtr->GetPixel(sampleIdx), bestLabel);
			      }
			      samplePtr->SetPixel(sampleIdx, bestLabel);
			      anyChangePixel ++;
			 }
		    } // sampleIdx[1]
	       } // sampleIdx[0]

	       scan++;

	       if (m_verbose >= 2) {
		    printf("icmEstep(): MC Sample = %i, scan = %i, anyChaingePixel = %i.\n", sampleIdx[2], scan, anyChangePixel);
		    saveExtractSlice(samplePtr, beginSliceIdx, "onesample.nii");
	       }
	  } // do many scans.
	  while (anyChangePixel > 0 && scan <= maxScan);
     } // sampleIdx[2]
     return 0;    

}


int MCModel::sortLabels(ImageType3D::Pointer samplePtr)
{
     unsigned i = 0, j = 0;
     unsigned clsIdx = 0;
     float smallestMean = 0;
     unsigned smallestLabel = 0;
     float tempMean = 0;
     float tempSigma = 0;
     unsigned int tempNumPoints = 0;

     for (i = 0; i < m_numClusters; i ++) {
	  smallestMean = 1000000;
	  for (j = i; j < m_numClusters; j ++) {
	       if (m_cls[j].mu < smallestMean) {
		    smallestMean = m_cls[j].mu;
		    smallestLabel = j;
	       }
	       
	  }

	  tempMean = m_cls[smallestLabel].mu;
	  tempSigma = m_cls[smallestLabel].sigma;
	  tempNumPoints = m_cls[smallestLabel].numPoints;

	  m_cls[smallestLabel].mu = m_cls[i].mu;
	  m_cls[smallestLabel].sigma = m_cls[i].sigma;
	  m_cls[smallestLabel].numPoints = m_cls[i].numPoints;

	  m_cls[i].mu = tempMean;
	  m_cls[i].sigma = tempSigma;
	  m_cls[i].numPoints = tempNumPoints;


     }

}

int MCModel::estLabelML(ImageType3D::Pointer samplePtr,
	       ImageType2D::Pointer imagePtr)
{
     float pixelInt = 0;
     ImageType2D::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::SizeType sampleSize =  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;
     ImageType2D::IndexType imageIdx;
     
     unsigned  clsIdx = 0;
     
     unsigned char bestLabel = 0;
     float thisDistance = 0, smallestDistance = 100000;
     sampleIdx[2] = sampleSize[2] - 1; // write to last slice.
     for (imageIdx[0] = 0; imageIdx[0] < imageSize[0]; imageIdx[0] ++) {
	  for (imageIdx[1] = 0; imageIdx[1] < imageSize[1]; imageIdx[1] ++) {	  
	       smallestDistance = 100000;
	       pixelInt = imagePtr->GetPixel(imageIdx);
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
		    thisDistance = pow((pixelInt - m_cls[clsIdx].mu), 2)/(2*pow(m_cls[clsIdx].sigma, 2)) + log(m_cls[clsIdx].sigma);
		    if (thisDistance < smallestDistance) {
			 smallestDistance = thisDistance;
			 bestLabel = clsIdx;
		    }
	       }
	       // write back.
	       sampleIdx[0] = imageIdx[0];
	       sampleIdx[1] = imageIdx[1];
	       samplePtr->SetPixel(sampleIdx, bestLabel);
	  }
     }

}

int MCModel::annealing(ImageType2D::Pointer imagePtr, ImageType3D::Pointer samplePtr, unsigned maxRun)
{
     int scan = 0;
     double p_acpt = 0;
     double denergy = 0;
     int cand = 1;
     int label = 0;
     int nLabel = 0; // neighobr label
     float intensity = 0; // current pixel intensity.
     unsigned run = 0;
     float initTemp = 1;
     float finalTemp = 0.1;
     float curTemp = 1;

     ImageType3D::IndexType nidx;    
     ImageType2D::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::SizeType sampleSize =  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;
     ImageType3D::IndexType sampleDestIdx;     
     ImageType2D::IndexType imageIdx;


     ///////        Define random number generator.        ////////

     // Uniform integer generator.
     boost::uniform_int<> uni_int(1, m_numClusters); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);


     // sampling Markov Random Field. 
     scan = 0;
     sampleIdx[2] = sampleSize[2] - 1;
     nidx[2] = sampleSize[2] - 1;

     for (run = 0; run < maxRun; run ++) {
	  curTemp = initTemp * pow((finalTemp/initTemp), (float)run/(float)maxRun);
	  
     	  for (sampleIdx[0] = 0; sampleIdx[0] < sampleSize[0]; sampleIdx[0]++) {
     	       for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1]++) {
		    cand = roll_die() - 1;
		    label = samplePtr->GetPixel(sampleIdx);
     		    denergy = 0;
     		    if (sampleIdx[0] > 0) {
			 nidx[0] = sampleIdx[0] - 1;
			 nidx[1] = sampleIdx[1];
			 nLabel = samplePtr->GetPixel(nidx); // neighbor label
			 denergy = denergy + int(cand != nLabel) - int(label != nLabel);
		    }

     		    if (sampleIdx[0] < sampleSize[0]-1) {
			 nidx[0] = sampleIdx[0] + 1;
			 nidx[1] = sampleIdx[1];
			 nLabel = samplePtr->GetPixel(nidx); // neighbor label
			 denergy = denergy + int(cand != nLabel) - int(label != nLabel);
		    }

     		    if (sampleIdx[1] > 0) {
			 nidx[0] = sampleIdx[0];
			 nidx[1] = sampleIdx[1] - 1;
			 nLabel = samplePtr->GetPixel(nidx); // neighbor label
			 denergy = denergy + int(cand != nLabel) - int(label != nLabel);
		    }

     		    if (sampleIdx[1] < sampleSize[1]-1) {
			 nidx[0] = sampleIdx[0];
			 nidx[1] = sampleIdx[1] + 1;
			 nLabel = samplePtr->GetPixel(nidx); // neighbor label
			 denergy = denergy + int(cand != nLabel) - int(label != nLabel);
		    }

     		    denergy = m_beta * denergy;

		    // likelihood energy.
		    imageIdx[0] = sampleIdx[0];
		    imageIdx[1] = sampleIdx[1];
		    intensity = imagePtr->GetPixel(imageIdx);
		    denergy = denergy
			 + (intensity - m_cls[cand].mu)*(intensity - m_cls[cand].mu)/(2*m_cls[cand].sigma * m_cls[cand].sigma) + log(m_cls[cand].sigma) // Candidate's likelihood energy.
			 - (intensity - m_cls[label].mu)*(intensity - m_cls[label].mu)/(2*m_cls[label].sigma * m_cls[label].sigma) - log(m_cls[label].sigma); // Current label's likelihood energy.
		    
		    denergy = denergy/curTemp;
     		    // if energy change less than zero, just accept
     		    // candidate. otherwise accept with exp(- energy
     		    // change).
		    
     		    if (denergy <= 0) {
			 
			 samplePtr->SetPixel(sampleIdx, cand);
     		    }
     		    else {
     			 p_acpt = exp(-denergy);
     			 if (uni() < p_acpt) {
			      samplePtr->SetPixel(sampleIdx, cand);
     			 }
     		    }
     	       }
     	  }

	  if (run <= 0) {
	       printf("Annealing Scan %i...", run);
	  }
	  else {
	       printf("%i...", run);
	  }
	  
     }
     printf("\n");
     return 0;    
}

int MCModel::estimateLabelsGraphcuts(ImageType2D::Pointer imagePtr, 
				     ImageType3D::Pointer samplePtr)
{
     ImageType2D::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::SizeType sampleSize =  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;
     ImageType2D::IndexType imageIdx;
     float imageInt = 0;

     unsigned short clsIdx = 0, clsIdx2 = 0;
     double dataEnergy = 0;

     GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(imageSize[0], imageSize[1], m_numClusters);

     // set initial labels.
     sampleIdx[2] = sampleSize[2] - 1;
     for (imageIdx[0] = 0; imageIdx[0] < imageSize[0]; imageIdx[0]++) {
     	  for (imageIdx[1] = 0; imageIdx[1] < imageSize[1]; imageIdx[1]++) {
     	       sampleIdx[0] = imageIdx[0];
     	       sampleIdx[1] = imageIdx[1];
     	       gc->setLabel(imageIdx[0] * imageSize[1] + imageIdx[1], samplePtr->GetPixel(sampleIdx));
     	  }
     }

     // Specify data cost.
     sampleIdx[2] = sampleSize[2] - 1;
     for (imageIdx[0] = 0; imageIdx[0] < imageSize[0]; imageIdx[0]++) {
	  for (imageIdx[1] = 0; imageIdx[1] < imageSize[1]; imageIdx[1]++) {
	       imageInt = imagePtr->GetPixel(imageIdx);
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
		    dataEnergy = pow((imageInt - m_cls[clsIdx].mu), 2)/(2*pow(m_cls[clsIdx].sigma, 2)) + log(m_cls[clsIdx].sigma);
     		    gc->setDataCost(imageIdx[0] * imageSize[1] + imageIdx[1], clsIdx, dataEnergy);
	       }
	  }
     }

     // Sepcify smoothness cost.

     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
     	  for (clsIdx2 = 0; clsIdx2 < m_numClusters; clsIdx2++) {
     	       gc->setSmoothCost(clsIdx, clsIdx2, (clsIdx == clsIdx2)? 0: m_beta); 
     	  }
     }

     printf("estimateLabelsGraphcuts(): Before optimization. Total energy: %f, data enerty: %f, smooth energy: %f, label energy: %f.\n",gc->compute_energy(), gc->giveDataEnergy(), gc->giveSmoothEnergy(), gc->giveLabelEnergy());

     // run expansion for 2 iterations. 
     gc->expansion(1);
     gc->swap(1);
     gc->expansion(1);
     gc->swap(1);

     printf("estimateLabelsGraphcuts(): After expansion 1, Total energy: %f, data enerty: %f, smooth energy: %f, label energy: %f.\n ",gc->compute_energy(), gc->giveDataEnergy(), gc->giveSmoothEnergy(), gc->giveLabelEnergy());

     // write back.
     sampleIdx[2] = sampleSize[2] - 1;
     for (sampleIdx[0] = 0; sampleIdx[0] < sampleSize[0]; sampleIdx[0]++) {
	  for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1]++) {

	       samplePtr->SetPixel(sampleIdx, gc->whatLabel(sampleIdx[0] * sampleSize[1] + sampleIdx[1]));
	  }
     }


     delete gc;
}

int MCModel::initWithTruth(float sigma, float offset, float scale, float beta)
{

     for (unsigned short clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
	  m_cls[clsIdx].mu = clsIdx * scale + offset;
	  m_cls[clsIdx].sigma = sigma;
     }
     m_beta = beta;
}


int MCModel::nChainInit(ImageType3D::Pointer samplePtr, std::string method)
{
     
     ImageType3D::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;     
     ImageType3D::IndexType sampleDestIdx;          

     if (method.compare("ml") == 0) {
	  sampleIdx[2] = sampleSize[2] - 1;
     }
     else if (method.compare("largestll") == 0) {
	  // choose the one with largest likelihood.
	  unsigned short maxLLIdx = 0;
	  double maxLL = -10000000;
	  for (unsigned mcSampleIdx = 0; mcSampleIdx < m_numScan; mcSampleIdx++) {
	       if (m_singleMCSampleLLNew[mcSampleIdx] > maxLL) {
		    maxLL = m_singleMCSampleLLNew[mcSampleIdx];
		    maxLLIdx = mcSampleIdx;
	       }
	  }
	  sampleIdx[2] = maxLLIdx;

	  if (m_verbose >=1) {
	       printf("nChainInit(): MC sample %i is used for init value of all samples. \n", sampleIdx[2]);
	  }
     }
     for (sampleIdx[0] = 0; sampleIdx[0] < sampleSize[0]; sampleIdx[0]++) {
	  for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1]++) {
	       sampleDestIdx[0] = sampleIdx[0];
	       sampleDestIdx[1] = sampleIdx[1];
	       for (sampleDestIdx[2] = 0; sampleDestIdx[2] < m_numScan; sampleDestIdx[2] ++) {
		    samplePtr->SetPixel(sampleDestIdx, samplePtr->GetPixel(sampleIdx));
	       }
	  }
     }


}

int MCModel::nChainSampling(ImageType2D::Pointer imagePtr,
		   ImageType3D::Pointer samplePtr,
		   unsigned int burnin)
{
     int scan = 0;
     double p_acpt = 0;
     double denergy = 0;
     int cand = 1;
     int label = 0;
     int nLabel = 0; // neighobr label
     float intensity = 0; // current pixel intensity.

     ImageType3D::IndexType nidx;    
     ImageType2D::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::SizeType sampleSize =  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;
     ImageType3D::IndexType sampleDestIdx;     
     ImageType2D::IndexType imageIdx;


     //////////////////////////////////////////////////////////////
     ///////        Define random number generator.        ////////
     //////////////////////////////////////////////////////////////

     // Uniform integer generator.
     boost::uniform_int<> uni_int(1, m_numClusters); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);


     // sampling Markov Random Field. 
     scan = 0;
     nidx[2] = sampleSize[2] - 1;

     for (sampleIdx[2] = 0; sampleIdx[2] < m_numScan; sampleIdx[2]++) {
	  nidx[2] = sampleIdx[2];
	  for (scan = 0; scan < burnin + 1; scan ++) {
	       for (sampleIdx[0] = 0; sampleIdx[0] < sampleSize[0]; sampleIdx[0]++) {
		    for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1]++) {
			 cand = roll_die() - 1;
			 label = samplePtr->GetPixel(sampleIdx);
			 denergy = 0;
			 if (sampleIdx[0] > 0) {
			      nidx[0] = sampleIdx[0] - 1;
			      nidx[1] = sampleIdx[1];
			      nLabel = samplePtr->GetPixel(nidx); // neighbor label
			      denergy = denergy + int(cand != nLabel) - int(label != nLabel);
			 }

			 if (sampleIdx[0] < sampleSize[0]-1) {
			      nidx[0] = sampleIdx[0] + 1;
			      nidx[1] = sampleIdx[1];
			      nLabel = samplePtr->GetPixel(nidx); // neighbor label
			      denergy = denergy + int(cand != nLabel) - int(label != nLabel);
			 }

			 if (sampleIdx[1] > 0) {
			      nidx[0] = sampleIdx[0];
			      nidx[1] = sampleIdx[1] - 1;
			      nLabel = samplePtr->GetPixel(nidx); // neighbor label
			      denergy = denergy + int(cand != nLabel) - int(label != nLabel);
			 }

			 if (sampleIdx[1] < sampleSize[1]-1) {
			      nidx[0] = sampleIdx[0];
			      nidx[1] = sampleIdx[1] + 1;
			      nLabel = samplePtr->GetPixel(nidx); // neighbor label
			      denergy = denergy + int(cand != nLabel) - int(label != nLabel);
			 }

			 denergy = m_beta * denergy;

			 // likelihood energy.
			 imageIdx[0] = sampleIdx[0];
			 imageIdx[1] = sampleIdx[1];
			 intensity = imagePtr->GetPixel(imageIdx);
			 denergy = denergy
			      + (intensity - m_cls[cand].mu)*(intensity - m_cls[cand].mu)/(2*m_cls[cand].sigma * m_cls[cand].sigma) + log(m_cls[cand].sigma) // Candidate's likelihood energy.
			      - (intensity - m_cls[label].mu)*(intensity - m_cls[label].mu)/(2*m_cls[label].sigma * m_cls[label].sigma) - log(m_cls[label].sigma); // Current label's likelihood energy.
		    
			 // if energy change less than zero, just accept
			 // candidate. otherwise accept with exp(- energy
			 // change).
			 
			 if (denergy <= 0) {
			 
			      samplePtr->SetPixel(sampleIdx, cand);
			 }
			 else {
			      p_acpt = exp(-denergy);
			      if (uni() < p_acpt) {
				   samplePtr->SetPixel(sampleIdx, cand);
			      }
			 }
		    } // sampleIdx[1]
	       } // sampleIdx[0]
	       
	  } // scan
	  
	  if (m_verbose >= 2) {
	       printf("nChainSamping(): MC sample %i done.\n", sampleIdx[2]);
	  }
     } // mc sample

     m_oldCll = m_cll;
     m_cll = eval_cll(imagePtr, samplePtr);

     m_oldPriorll = priorll_;
     priorll_ = eval_ll(samplePtr, m_beta);

     return 0;    
}

int MCModel::updateSingleSampleLL(ImageType2D::Pointer imagePtr,
				  ImageType3D::Pointer samplePtr)
{
     unsigned clsIdx = 0;
     unsigned mcSampleIdx = 0;

     m_singleMCSampleLLNew[mcSampleIdx] = 0;
     for(mcSampleIdx = 0; mcSampleIdx < m_numScan; mcSampleIdx ++) {
	  m_singleMCSampleLLNew[mcSampleIdx] = singleMCSampleLL(imagePtr, samplePtr, mcSampleIdx, m_beta, m_cls, m_numClusters);
	  m_singleMCSampleLLOld[mcSampleIdx] = singleMCSampleLL(imagePtr, samplePtr, mcSampleIdx, m_oldBeta, m_oldCls, m_numClusters);
	  if (m_verbose >= 2) {
	       printf("singleMCSampleLLNew = %f, singleMCSampleLLOld = %f, oneSampleLambda = %f.\n", m_singleMCSampleLLNew[mcSampleIdx], m_singleMCSampleLLOld[mcSampleIdx],
		    m_singleMCSampleLLNew[mcSampleIdx] - m_singleMCSampleLLOld[mcSampleIdx]);
	  }
     }

}

double MCModel::computeASE(ImageType2D::Pointer imagePtr,
			   ImageType3D::Pointer samplePtr)
{
     unsigned clsIdx = 0;
     unsigned mcSampleIdx = 0;

     vnl_vector<float>  oneSampleLambda;
     double varLambda = 0;
     double meanLambda;

     oneSampleLambda = m_singleMCSampleLLNew - m_singleMCSampleLLOld;
     meanLambda = oneSampleLambda.mean();

     for(mcSampleIdx = 0; mcSampleIdx < m_numScan; mcSampleIdx ++) {
	  varLambda = varLambda + oneSampleLambda[mcSampleIdx] * oneSampleLambda[mcSampleIdx]/m_numScan;
     }

     varLambda = varLambda - meanLambda * meanLambda;
     if (m_verbose >= 1) {
	  printf("meanLambda = %f, varLambda = %f.\n", meanLambda, varLambda);
     }

     return sqrt(varLambda/m_numScan);
}

int MCModel::ascentAcceptEstimation(ImageType3D::Pointer samplePtr,
				    float signLevel,
				    float ASE,
				    float addSamplePercent)
{
     using namespace boost::math;

     float criticalValue = 0;
     unsigned short clsIdx = 0;
     double deltaQ = 0;
     unsigned newNumScan = 0;

     ImageType3D::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleSrcIdx;     
     ImageType3D::IndexType sampleDestIdx;          

     
     // Given signLevel, compute critical value.
     normal_distribution<> normal_dist(0, 1); // Normal distribution.
     criticalValue = quantile(normal_dist, 1 - signLevel); // i.e. z_alpha.

     deltaQ = (priorll_ + m_cll) - (m_oldPriorll + m_oldCll);

     printf("deltaQ = %f, deltaQ - criticalValue * ASE = %f.\n",  deltaQ, deltaQ - criticalValue * ASE);
     if ( deltaQ - criticalValue * ASE < 0 && m_numScan < sampleSize[2]) {
	  // Discard parameters obtained in this EM iteration. Restore
	  // old parameters.
	  m_beta = m_oldBeta;
	  for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
	       m_cls[clsIdx].mu = m_oldCls[clsIdx].mu;
	       m_cls[clsIdx].sigma = m_oldCls[clsIdx].sigma;
	       m_cls[clsIdx].numPoints = m_oldCls[clsIdx].numPoints;
	  }
	  
	  // 
	  priorll_ = m_oldPriorll;
	  m_cll = m_oldCll;
	  printf("ascentAddSamples(): DO NOT accept this estimation. Restore parameters to:\n");
	  print("normal");

	  // Init new sample with existing samples, and increase MC sample size.
	  newNumScan = vnl_math_min(unsigned (m_numScan + m_numScan * addSamplePercent), unsigned(sampleSize[2]));

	  for (sampleDestIdx[2] = m_numScan; sampleDestIdx[2] < newNumScan; sampleDestIdx[2] ++) {
	       sampleSrcIdx[2] = sampleDestIdx[2] - m_numScan;
	       for (sampleDestIdx[0] = 0; sampleDestIdx[0] < sampleSize[0]; sampleDestIdx[0]++) {
		    for (sampleDestIdx[1] = 0; sampleDestIdx[1] < sampleSize[1]; sampleDestIdx[1]++) {
			 sampleSrcIdx[0] = sampleDestIdx[0];
			 sampleSrcIdx[1] = sampleDestIdx[1];
			 samplePtr->SetPixel(sampleDestIdx, samplePtr->GetPixel(sampleSrcIdx));
		    }
	       }
	  }
	  
	  m_numScan = newNumScan;
	  m_singleMCSampleLLNew.set_size(m_numScan);
	  m_singleMCSampleLLOld.set_size(m_numScan);

	  printf("ascentAddSamples(): Increase MC sample size to: %d.\n", m_numScan);
	  return 0;
     }
     else {
	  printf("ascentAddSamples(): parameters estimation accepted. \n");
	  return 1;
     }
}

int MCModel::ascentConvergeTest(ImageType3D::Pointer samplePtr,float signLevel, float ASE, float minllinc)
{
     using namespace boost::math;
     float criticalValue = 0;
     double deltaQ = 0;

     ImageType3D::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     
     // Given signLevel, compute critical value.
     normal_distribution<> normal_dist(0, 1); // Normal distribution.
     criticalValue = quantile(normal_dist, 1 - signLevel); // i.e. z_alpha.

     deltaQ = (priorll_ + m_cll) - (m_oldPriorll + m_oldCll);

     if (m_verbose >=0) {
     printf("ascentConvergeTest(): deltaQ = %f,  deltaQ + criticalValue * ASE = %f, min log-likelihood increase = %f.\n", deltaQ, deltaQ + criticalValue * ASE, minllinc);
     }

     if ( deltaQ + criticalValue * ASE < minllinc) {
	  printf("ascentConvergeTest(): joint log-likelihood increasement not big enough. EM should be converged. \n");
	  return (1);
     }
     else {
	  return (0);
     }
}
     
double singleMCSampleLL(ImageType2D::Pointer imagePtr,
			ImageType3D::Pointer samplePtr,
			unsigned short mcSampleIdx,
			float beta,
			CompType * const cls,
			unsigned short numClusters)
{
     ImageType3D::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;     
     ImageType2D::IndexType imageIdx;     
     ImageType3D::IndexType nidx;     

     double cll = 0;
     double priorll = 0;
     int clsIdx = 0;
     
     sampleIdx[2] = mcSampleIdx;
     for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1] ++) {
	  for (sampleIdx[0] = 0; sampleIdx[0] < sampleSize[0]; sampleIdx[0] ++) {
	       imageIdx[0] = sampleIdx[0];
	       imageIdx[1] = sampleIdx[1];
	       clsIdx = samplePtr->GetPixel(sampleIdx);
	       cll = cll +  ((imagePtr->GetPixel(imageIdx) - cls[clsIdx].mu)*(imagePtr->GetPixel(imageIdx) - cls[clsIdx].mu))/(-2 * cls[clsIdx].sigma * cls[clsIdx].sigma) - log(cls[clsIdx].sigma);
	  } 
     }

     int S = 0, curS = 0; //  neighbors sum.
     double sumexp = 0;

     int thisLabel = 0;
     int curLabel = 0;
     
     nidx[2] = sampleIdx[2];
     for (sampleIdx[0] = 0; sampleIdx[0] < sampleSize[0]; sampleIdx[0] ++) {
	  for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1] ++) {
	       curLabel = samplePtr->GetPixel(sampleIdx);
	       sumexp = 0;
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
		    if (thisLabel == curLabel) {
			 curS = S;
		    }
	       }
			 
	       priorll += (-beta * curS - log(sumexp));
	  }
     }
     
     return (cll + priorll);
}
