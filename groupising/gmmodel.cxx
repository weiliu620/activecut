#include "common_gmmodel.h"
#include "gmmodel.h"
#include <boost/math/distributions/students_t.hpp>
extern twister_base_gen_type mygenerator; 

GMModel::GMModel(VnlVectorImageType::Pointer obsPtr,
		 ImageType3DVecF::Pointer grpPtr,
		 ImageType3DChar::Pointer maskPtr,
		 unsigned numClusters,
		 unsigned short verbose )
{
     unsigned subIdx = 0, clsIdx = 0;

     // init members.
     m_verbose = verbose;

     // get number of subjects (from a pixel of obsPtr).
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();
     VnlVectorImageType::IndexType obsIdx;
     obsIdx.Fill( 0 );
     vnl_vector<float> obsVector = obsPtr->GetPixel(obsIdx);
     m_numSubs = obsVector.size();
     m_numClusters = numClusters;

     // initial total number of data points (voxels) in the mask.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());
     ImageType3DChar::IndexType maskIdx;     
     m_gmm.totalPts = 0;
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       m_gmm.totalPts ++;
	  }  // maskIt > 0
     } // maskIt

     if(m_verbose >= 1) {
	  printf("total number of voxels within mask: %i.\n", m_gmm.totalPts);
     }



     // init m_gmm
     m_gmm.comp.resize(m_numClusters);
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_gmm.comp[clsIdx].prop = 1 / float (m_numClusters);
	  m_gmm.comp[clsIdx].numPoints = m_gmm.totalPts / float (m_numClusters);
	  m_gmm.comp[clsIdx].mu.set_size(m_numSubs);
	  m_gmm.comp[clsIdx].sigma.set_size(m_numSubs);
     }
     
}

GMModel::~GMModel()
{

}

int GMModel::printSelf(std::string format)
{
     unsigned int clsIdx = 0, subIdx = 0;

     if (format.compare("normal") == 0) {

	  for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
	       printf("comp[%i]: prop = %f, numPoints = %f\n", 
		      clsIdx, 
		      m_gmm.comp[clsIdx].prop, 
		      m_gmm.comp[clsIdx].numPoints);

	       printf("  mu = ");
	       for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
		    printf("%2.2f ", m_gmm.comp[clsIdx].mu[subIdx]);
	       }
	       printf("\n");
	       printf("  sigma = ");
	       for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
		    printf("%2.2f ", m_gmm.comp[clsIdx].sigma[subIdx]);
	       }
	       printf("\n");
	  } // clsIdx.
     }
     else if (format.compare("table") == 0) {
     }
     else {
	  printf("GMModel::print(): print format not recognized.\n");
     }
     fflush(stdout);
     
     return 0;
}


int GMModel::Estep(ImageType3DVecF::Pointer grpPtr,
		   VnlVectorImageType::Pointer obsPtr,
		   ImageType3DChar::Pointer maskPtr)
{
     unsigned clsIdx = 0, subIdx = 0;

     VnlVectorImageType::IndexType obsIdx;    
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();
     IteratorTypeVnlVector obsIt(obsPtr, obsPtr->GetLargestPossibleRegion() );
     
     vnl_vector <float> obsVector(m_numSubs, 0);

     //mask
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType maskIdx;     

     // group map.
     IteratorType3DVecF grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );
     itk::VariableLengthVector<float> grpVector(m_numClusters);

     double const Pi = 4 * atan(1);
     double M = 0, logGamma_znk = 0, logsumexp = 0;
     vnl_vector<float> expterm(m_numClusters, 0);

     for (grpIt.GoToBegin(), obsIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ grpIt, ++ maskIt, ++ obsIt) {
	  if (maskIt.Get() > 0) {
	       maskIdx = maskIt.GetIndex();
	       obsVector = obsIt.Get();

	       // deubg.
	       if (maskIdx[0] == 59 && maskIdx[1] == 17) {
		    printf("\n");
	       }

	       for (unsigned clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
		    expterm[clsIdx] = 0;
		    for (unsigned subIdx = 0; subIdx < m_numSubs; subIdx ++) {
			 expterm[clsIdx] +=  -log(sqrt(2*Pi)) - log( m_gmm.comp[clsIdx].sigma[subIdx]) -  pow((obsVector[subIdx]-m_gmm.comp[clsIdx].mu[subIdx]), 2) / ( 2 * pow(m_gmm.comp[clsIdx].sigma[subIdx], 2));
		    }// subIdx
	       } // clsIdx
	       M = expterm.min_value();

	       // compute logsumexp
	       logsumexp = 0;
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
		    logsumexp += m_gmm.comp[clsIdx].prop * exp (expterm[clsIdx] - M);
	       }
	       logsumexp = log (logsumexp);

	       for (unsigned clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
		    logGamma_znk = log (m_gmm.comp[clsIdx].prop) + expterm[clsIdx] - M - logsumexp;
		    grpVector[clsIdx] = exp(logGamma_znk);
	       }

	       // save back group map
	       grpIt.Set( grpVector );
	  } // maskIt > 0
     } // neiGrpIt
     return 0;
}

int GMModel::EstimateMu(ImageType3DVecF::Pointer grpPtr,
			   VnlVectorImageType::Pointer obsPtr,
			   ImageType3DChar::Pointer maskPtr)

{
     unsigned int clsIdx = 0;

     // images
     VnlVectorImageType::IndexType obsIdx;    
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();
     IteratorTypeVnlVector obsIt(obsPtr, obsPtr->GetLargestPossibleRegion() );

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DChar::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();

     // grps
     ImageType3DVecF::IndexType grpIdx;
     itk::VariableLengthVector<float> grpVector(m_numClusters);
     IteratorType3DVecF grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );

     // reset all mu to zero.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_gmm.comp[clsIdx].numPoints = 0;
	  m_gmm.comp[clsIdx].mu = 0;
     }


     for (maskIt.GoToBegin(), obsIt.GoToBegin(), grpIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++ obsIt, ++ grpIt) {
	  if (maskIt.Get() > 0) {
	       grpVector = grpIt.Get();
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
		    m_gmm.comp[clsIdx].mu += grpVector[clsIdx] * obsIt.Get();
		    m_gmm.comp[clsIdx].numPoints += grpVector[clsIdx];
	       } // clsIdx
	  } // maskIt > 0
     } //maskIt  

     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_gmm.comp[clsIdx].mu = m_gmm.comp[clsIdx].mu / m_gmm.comp[clsIdx].numPoints;
	  m_gmm.comp[clsIdx].prop = m_gmm.comp[clsIdx].numPoints / m_gmm.totalPts;
     }
     return 0;
}

int GMModel::estimateSigma(ImageType3DVecF::Pointer grpPtr,
			   VnlVectorImageType::Pointer obsPtr,
			   ImageType3DChar::Pointer maskPtr)
{
     unsigned clsIdx = 0;

     // images
     VnlVectorImageType::IndexType obsIdx;    
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();
     IteratorTypeVnlVector obsIt(obsPtr, obsPtr->GetLargestPossibleRegion() );

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DChar::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();

     // grps
     ImageType3DVecF::IndexType grpIdx;
     itk::VariableLengthVector<float> grpVector(m_numClusters);
     IteratorType3DVecF grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );

     // reset all sigma
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_gmm.comp[clsIdx].sigma = 0;
     }

     vnl_vector<float> singleVar(m_numSubs, 0);
     vnl_vector<float> variance(m_numSubs, 0);

     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  variance = 0;
	  for (maskIt.GoToBegin(), obsIt.GoToBegin(), grpIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++ obsIt, ++ grpIt) {
	       if (maskIt.Get() > 0) {
		    grpVector = grpIt.Get();
		    // the sample variance for each subject is
		    // independent. so just compute sigma for each subject
		    // saparately. This is like multivariate Gaussian with
		    // dimension equal to num of subjects, and covariance
		    // matrix is diagonal.
		    singleVar = obsIt.Get() - m_gmm.comp[clsIdx].mu;

		    for (unsigned subIdx = 0; subIdx < m_numSubs; subIdx ++) {
			 variance[subIdx]  = variance[subIdx] + grpVector[clsIdx] * singleVar[subIdx] * singleVar[subIdx];
		    } 
	       }  // maskIt > 0
	  } // maskIt

	  // normalize.
	  for (unsigned subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	  m_gmm.comp[clsIdx].sigma[subIdx] = sqrt ( variance[subIdx] / m_gmm.comp[clsIdx].numPoints);
	  } // subIdx
	  
     } // clsIdx

     return 0;
}


double GMModel::llhood(ImageType3DVecF::Pointer grpPtr,
		       VnlVectorImageType::Pointer obsPtr,
		       ImageType3DChar::Pointer maskPtr)
{
     unsigned subIdx = 0, clsIdx = 0;;
     VnlVectorImageType::IndexType obsIdx;    
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();
     IteratorTypeVnlVector obsIt(obsPtr, obsPtr->GetLargestPossibleRegion() );
     vnl_vector <float> obsVector(obsSize[3], 0);

     //mask
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType maskIdx;     

     // group map.
     IteratorType3DVecF grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );
     itk::VariableLengthVector<float> grpVector(m_numClusters);

     double const Pi = 4 * atan(1);

     // sum of sigma of all subjects.
     vnl_vector<float> sumOfLogSigma(m_numClusters, 0);
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       sumOfLogSigma[clsIdx] += log(m_gmm.comp[clsIdx].sigma[subIdx]);
	  }
     }

     double jlogll = 0;
     for (grpIt.GoToBegin(), maskIt.GoToBegin(), obsIt.GoToBegin(); 
	  !maskIt.IsAtEnd(); 
	  ++ grpIt, ++ maskIt, ++ obsIt) {
	  if (maskIt.Get() > 0) {
	       grpVector = grpIt.Get();

	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
		    float condll = - m_numSubs * log (sqrt(2*Pi)) - sumOfLogSigma[clsIdx];
			 for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
			      condll += pow((obsIt.Get()[subIdx]  - m_gmm.comp[clsIdx].mu[subIdx]), 2) / (2* pow(m_gmm.comp[clsIdx].sigma[subIdx], 2));
			 }
		    
		    jlogll += grpVector[clsIdx] * (log(m_gmm.comp[clsIdx].prop) + condll );
						   
	       } // clsIdx
	  } //maskIt > 0
     } // groupIt
	       
}
