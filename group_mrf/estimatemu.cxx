#include "commonalt.h"
#include "MCModel_init.h"

// Given group's labels, estimate mu and numPoints.
int MCModel::estimateMu(ImageType3DChar::Pointer groupPtr,
			ImageType5DFloat::Pointer imagePtr,
			ImageType3DChar::Pointer maskPtr)
{
     unsigned int clsIdx = 0, subIdx = 0;

     // reset all mu to zero.
     for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	  for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	       m_vmm[subIdx].comp[clsIdx].numPoints = 0;
	       m_vmm[subIdx].comp[clsIdx].mu = 0;
	  }
     }

     // compute mu and numPoints for all sub, all clusters.
#pragma omp parallel for
     for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	  estimateMuSub(groupPtr, imagePtr, maskPtr, subIdx);	  
	  if (m_verbose >= 2) {
	       printf("estimateMu(): estimabeMuSub() done for subject %i.\n", subIdx);
	  }
     } // for subIdx

     return 0;
}


int MCModel::estimateMuSub(ImageType3DChar::Pointer groupPtr,
			   ImageType5DFloat::Pointer imagePtr,
			   ImageType3DChar::Pointer maskPtr,
			   unsigned subIdx)
{
     unsigned int clsIdx = 0;

     // images
     ImageType5DFloat::IndexType imageIdx;    
     ImageType5DFloat::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     ImageType3DChar::RegionType maskRegion = maskPtr->GetLargestPossibleRegion();
     IteratorType3DChar maskIt(maskPtr, maskRegion);
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();

     // groups
     ImageType3DChar::IndexType groupIdx;

     vnl_vector <float> timeSeries(imageSize[3], 0);

     // reset all mu to zero.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_vmm[subIdx].comp[clsIdx].numPoints = 0;
	  m_vmm[subIdx].comp[clsIdx].mu = 0;
     }

     imageIdx[4] = subIdx;
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       maskIdx = maskIt.GetIndex();

	       // First save time series at this position into a vector.
	       imageIdx[0] = maskIdx[0];
	       imageIdx[1] = maskIdx[1];
	       imageIdx[2] = maskIdx[2];
	       for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
		    timeSeries[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
	       }

	       groupIdx = maskIdx;
	       clsIdx = groupPtr->GetPixel(groupIdx);

	       m_vmm[subIdx].comp[clsIdx].mu += timeSeries;
	       m_vmm[subIdx].comp[clsIdx].numPoints ++;
	  }  // maskIt > 0
     } // for maskIt

     // compute mean time series and meanNorm.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  // meanNorm see page 1350 of "Clustering on the Unit
	  // Hypersphere using von Mises-Fisher distributions" by
	  // Banerjee.
	  m_vmm[subIdx].comp[clsIdx].meanNorm = m_vmm[subIdx].comp[clsIdx].mu.two_norm() / m_vmm[subIdx].comp[clsIdx].numPoints;
	  m_vmm[subIdx].comp[clsIdx].mu = m_vmm[subIdx].comp[clsIdx].mu.normalize();
     }
     return 0;
}
