#include "commonalt.h"
#include "MCModel_init.h"
#include <boost/math/distributions/students_t.hpp>
#include <utilalt.h>
extern twister_base_gen_type mygenerator;

MCModel::MCModel(ImageType5DFloat::Pointer imagePtr,
		 ImageType3DChar::Pointer groupPtr,
		 unsigned numClusters,
		 unsigned short verbose = 0)


{
     unsigned subIdx = 0, clsIdx = 0;

     // init members.
     m_verbose = verbose;
     m_numClusters = numClusters;
     m_temperature = 1;

     ImageType5DFloat::SizeType imageSize = 
	  imagePtr->GetLargestPossibleRegion().GetSize();
     m_numSubs = imageSize[4];

     // init m_vmm
     m_vmm.resize(m_numSubs); 
     for (subIdx = 0; subIdx < m_numSubs; subIdx++) {
	  m_vmm[subIdx].comp.resize(m_numClusters);
	  for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	       m_vmm[subIdx].comp[clsIdx].mu.set_size(imageSize[3]);
	       m_vmm[subIdx].comp[clsIdx].mu.fill(0);
	       m_vmm[subIdx].comp[clsIdx].meanNorm = 0;
	       m_vmm[subIdx].comp[clsIdx].kappa = 0;
	       m_vmm[subIdx].comp[clsIdx].label = 0;
	       m_vmm[subIdx].comp[clsIdx].numPoints = 0;
	  }
     }
}
     
MCModel::~MCModel()
{

}

int MCModel::printSelf(std::string format)
{
     unsigned int clsIdx = 0, subIdx = 0;

     if (format.compare("normal") == 0) {

	  printf("alpha = %f, beta_g = %f, beta_z = %f, temperature = %f\n",
		 m_alpha, m_beta_g, m_beta_z, m_temperature);

	  for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       printf("subject %i: \n", subIdx +1);
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
		    printf("       cluster[%i]: # of points: %ld, meanNorm = %f, kappa = %4.2f, label = %i. ", 
			   clsIdx + 1, 
			   m_vmm[subIdx].comp[clsIdx].numPoints, 
			   m_vmm[subIdx].comp[clsIdx].meanNorm, 
			   m_vmm[subIdx].comp[clsIdx].kappa,
			   m_vmm[subIdx].comp[clsIdx].label); 
		    printf("mu = ");
		    printVnlVector(m_vmm[subIdx].comp[clsIdx].mu, 3);
	       } // clsIdx
	  } // subIdx.
     }
     else if (format.compare("table") == 0) {
     }
     else {
	  printf("MCModel::print(): print format not recognized.\n");
     }
     
     return 0;
}



// Init cluster center by K-means.
double MCModel::kmeans(ImageType5DFloat::Pointer imagePtr, 
		       ImageType3DChar::Pointer groupPtr,
		       ImageType3DChar::Pointer maskPtr,
		       float converge_eps)
{
     unsigned subIdx = 0, clsIdx = 0;
     ImageType5DFloat::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();

     // Init kmeans by kmeans++ method.
     kmeanspp(imagePtr, maskPtr);

     // Define a prevmm to save previous vmm.
     std::vector <VMM> prevmm(m_numSubs);
     for (subIdx = 0; subIdx < m_numSubs; subIdx++) {
	  prevmm[subIdx].comp.resize(m_numClusters);
	  for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	       prevmm[subIdx].comp[clsIdx].mu.set_size(imageSize[3]);
	  }
     }

     double meanSSR = 0; 
     bool converge = 1;
     do {
     	  // update label assignment.
     	  meanSSR = estLabelsByMean(imagePtr, maskPtr, groupPtr);
	  
     	  // estimate mu. Save vmm before doing that.
	  prevmm = m_vmm;
     	  estimateMu(groupPtr, imagePtr, maskPtr);


	  for (subIdx = 0; subIdx < m_numSubs; subIdx++) { 
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) { 
		    m_vmm[subIdx].comp[clsIdx].mu.set_size(imageSize[3]);
		    if ((m_vmm[subIdx].comp[clsIdx].mu - prevmm[subIdx].comp[clsIdx].mu).two_norm() > converge_eps) {
			 converge = 0;
			 break;
		    }
	       } // clsIdx
	  } // subIdx
     }
     while(!KMeansConverge(m_vmm, prevmm));
     
     return meanSSR;
}

int MCModel::KMeansConverge(std::vector<VMM> & curvmm, std::vector<VMM> & prevmm)
{
     float eps = 0.01;
     unsigned maxChangePts = 500;
     unsigned numSubs = curvmm.size();
     unsigned numClusters = curvmm[0].comp.size();

     float normDiff = 0;
     unsigned changePts = 0;
     for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	  for (unsigned clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	       // label change?
	       changePts =  vnl_math_abs(curvmm[subIdx].comp[clsIdx].numPoints - prevmm[subIdx].comp[clsIdx].numPoints);
	       if (m_verbose >=2) {
		    printf("KmeansConverge(): sub[%i], cls[%i],  numPoints change = %i.\n", subIdx, clsIdx, changePts);
	       }
	       if (changePts > maxChangePts) {  return 0;     }

	       // mean time series vector change?
	       normDiff =  (curvmm[subIdx].comp[clsIdx].mu - prevmm[subIdx].comp[clsIdx].mu).two_norm();
	       if (m_verbose >=2) {
		    printf("sub[%i], cls[%i], normDiff (mean diff) = %f.\n", subIdx, clsIdx, normDiff);
	       }
	       if (normDiff > eps ) { return 0;  }

	  } // clsIdx
     } // subIdx

     return 1;
	       
     
}
int MCModel::kmeanspp(ImageType5DFloat::Pointer imagePtr,
		      ImageType3DChar::Pointer maskPtr)
{

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);


     unsigned int clsIdx = 0;
     float randNum = 0;
     float cdfValue = 0;
     float sumOfPdf = 0;
     unsigned subIdx = 0;

     // images
     ImageType5DFloat::IndexType imageIdx;    
     ImageType5DFloat::RegionType imageRegion;
     imageRegion = imagePtr->GetLargestPossibleRegion();
     ImageType5DFloat::SizeType imageSize = imageRegion.GetSize();

     // mask image.
     ImageType3DFloat::IndexType maskIdx;     
     ImageType3DFloat::RegionType maskRegion = maskPtr->GetLargestPossibleRegion();
     IteratorType3DChar maskIt(maskPtr, maskRegion);
     ImageType3DFloat::SizeType maskSize = maskRegion.GetSize();
     
     // Allocate memory for pdf of each data point.
     ImageType3DFloat::Pointer pdfPtr = ImageType3DFloat::New();
     ImageType3DFloat::IndexType pdfStart;
     ImageType3DFloat::IndexType pdfIdx;
     pdfStart.Fill(0);
     ImageType3DFloat::SizeType pdfSize;
     pdfSize[0] = imageSize[0];
     pdfSize[1] = imageSize[1];
     pdfSize[2] = imageSize[2];

     ImageType3DFloat::RegionType pdfRegion;
     pdfRegion.SetSize(pdfSize);
     pdfRegion.SetIndex(pdfStart);
     pdfPtr->SetRegions(pdfRegion);
     pdfPtr->Allocate();
     pdfPtr->FillBuffer(0);

     vnl_vector <float> timeSeries(imageSize[3], 0);

     // Find a random bogus center.
     do {
	  maskIdx[0] = floor(uni() * maskSize[0]);
	  maskIdx[1] = floor(uni() * maskSize[1]);
	  maskIdx[2] = floor(uni() * maskSize[2]);
     }
     while(maskPtr->GetPixel(maskIdx) <= 0);
     // OK I find the point, and assign it to mu_0 (for all sub)
     imageIdx[0] = maskIdx[0];
     imageIdx[1] = maskIdx[1];
     imageIdx[2] = maskIdx[2];
     for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	  imageIdx[4] = subIdx;
	  for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) { // time points.
	       m_vmm[subIdx].comp[0].mu[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
	  }
     }

     // Find all cluster center by Furthest-First rule. This also
     // includes first cluster. At last we throw away the bogus
     // center.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  // assume we have 1:clsIdx-1 mu available.

#pragma omp parallel for
	  for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       minDistPerSub(pdfPtr, imagePtr,maskPtr, clsIdx, subIdx);
	  }


	  // compute sum of pdf.
	  sumOfPdf = 0;
	  for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	       if (maskIt.Get() > 0) {
		    maskIdx = maskIt.GetIndex();
		    pdfIdx[0] = maskIdx[0];
		    pdfIdx[1] = maskIdx[1];
		    pdfIdx[2] = maskIdx[2];
		    sumOfPdf = sumOfPdf + pdfPtr->GetPixel(pdfIdx);
	       }
	  }

	  if (m_verbose >= 2) {
	       printf("sumOfPdf = %f\n", sumOfPdf);
	  }

	  // normalize pdf.
	  for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	       if (maskIt.Get() > 0) {
		    maskIdx = maskIt.GetIndex();
		    pdfIdx[0] = maskIdx[0];
		    pdfIdx[1] = maskIdx[1];
		    pdfIdx[2] = maskIdx[2];
		    pdfPtr->SetPixel(pdfIdx, pdfPtr->GetPixel(pdfIdx) / sumOfPdf);
	       }
	  }

	  randNum = uni();
	  cdfValue = 0;
	  maskIt.GoToBegin();
	       
	  while (randNum > cdfValue) {
	       while(maskIt.Get() <= 0) {
		    ++ maskIt;
	       }
	       maskIdx = maskIt.GetIndex();
	       cdfValue += pdfPtr->GetPixel( maskIdx );
	       ++ maskIt;
	  }
	  
	  // found the data point. Assign it to current cluster for all subs.
	  imageIdx[0] = maskIdx[0];
	  imageIdx[1] = maskIdx[1];
	  imageIdx[2] = maskIdx[2];

	  for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       imageIdx[4] = subIdx;
	       for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
		    m_vmm[subIdx].comp[clsIdx].mu[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
	       }
	  }

     } // clsIdx

}

double MCModel::estLabelsByMean(ImageType5DFloat::Pointer imagePtr,
				ImageType3DChar::Pointer maskPtr,
				ImageType3DChar::Pointer groupPtr)
{
     unsigned int clsIdx = 0;
     unsigned int subIdx = 0;

     // images
     ImageType5DFloat::RegionType imageRegion;
     imageRegion = imagePtr->GetLargestPossibleRegion();
     ImageType5DFloat::SizeType imageSize = imageRegion.GetSize();

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     ImageType3DChar::RegionType maskRegion = maskPtr->GetLargestPossibleRegion();
     IteratorType3DChar maskIt(maskPtr, maskRegion);
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();

     // group level label map.
     ImageType3DChar::IndexType groupIdx;
     ImageType3DChar::RegionType groupRegion;

     ImageType3DChar::SizeType groupSize = groupPtr->GetLargestPossibleRegion().GetSize();
     groupIdx[0] = 0;
     groupIdx[1] = 0;
     groupIdx[2] = 0;


     unsigned nearestClsLabel = 0;
     double nearestDistance = 1000000;

     double meanSSR = 0;
     unsigned long numAllPoints = 0;

     vnl_vector <float> timeSeries(imageSize[3], 0);

     // make a distance array 
     ImageType4DFloat::Pointer allSubDistPtr = ImageType4DFloat::New();
     ImageType4DFloat::IndexType allSubDistIdx;
     allSubDistIdx.Fill(0);
     ImageType4DFloat::SizeType allSubDistSize;
     allSubDistSize[0] = imageSize[0]; 
     allSubDistSize[1] = imageSize[1]; 
     allSubDistSize[2] = imageSize[2]; 
     allSubDistSize[3] = m_numClusters;

     ImageType4DFloat::RegionType allSubDistRegion;
     allSubDistRegion.SetSize(allSubDistSize);
     allSubDistRegion.SetIndex(allSubDistIdx);
     allSubDistPtr->SetRegions(allSubDistRegion);
     allSubDistPtr->Allocate();
     allSubDistPtr->FillBuffer(0);
     
     // 
#pragma omp parallel for
     for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	  estLabelsByMeanSub(imagePtr, maskPtr, allSubDistPtr, subIdx);
	  if (m_verbose >= 2 ) {
	       printf("estlabelByMeanSub() done for sub %i\n", subIdx);
	  }
     }

     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() > 0) {	  
	       numAllPoints ++;
	       maskIdx = maskIt.GetIndex();

	       allSubDistIdx[0] = maskIdx[0];
	       allSubDistIdx[1] = maskIdx[1];
	       allSubDistIdx[2] = maskIdx[2];

	       nearestDistance = 1000000;
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
		    allSubDistIdx[3] = clsIdx;
		    if (allSubDistPtr->GetPixel(allSubDistIdx) < nearestDistance) {
			 nearestClsLabel = clsIdx;
			 nearestDistance = allSubDistPtr->GetPixel(allSubDistIdx);
		    }
	       }

	       // Found the nearest cluster label.
	       groupIdx = maskIdx;
	       groupPtr->SetPixel(groupIdx, nearestClsLabel);
	       meanSSR += nearestDistance;
	  }
     }
     meanSSR = meanSSR / numAllPoints;
     return meanSSR;
}

int MCModel::minDistPerSub(ImageType3DFloat::Pointer pdfPtr, 
		       ImageType5DFloat::Pointer imagePtr,
		       ImageType3DChar::Pointer maskPtr,
		       unsigned int clsIdx,
		       unsigned int subIdx)
{
     // images
     ImageType5DFloat::IndexType imageIdx;    
     ImageType5DFloat::RegionType imageRegion;
     imageRegion = imagePtr->GetLargestPossibleRegion();
     ImageType5DFloat::SizeType imageSize = imageRegion.GetSize();

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     ImageType3DChar::RegionType maskRegion = maskPtr->GetLargestPossibleRegion();
     IteratorType3DChar maskIt(maskPtr, maskRegion);
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();

     ImageType3DFloat::IndexType pdfIdx;          
     vnl_vector <float> timeSeries(imageSize[3], 0);

     double minDistance = 0;
     double thisDistance = 0;
     int prevIdx = 0;
     imageIdx[4] = subIdx;

     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() > 0) {

	       // First save time series at this position into a vector.
	       maskIdx = maskIt.GetIndex();
	       imageIdx[0] = maskIdx[0];
	       imageIdx[1] = maskIdx[1];
	       imageIdx[2] = maskIdx[2];	  
	       for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
		    timeSeries[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
	       }
	       
	       minDistance = 1e10;
	       if (clsIdx == 0) {
		    minDistance = (timeSeries - m_vmm[subIdx].comp[0].mu).two_norm();
	       }
	       else {
		    for (prevIdx = 0; prevIdx < clsIdx; prevIdx ++) {
			 thisDistance = (timeSeries - m_vmm[subIdx].comp[prevIdx].mu).two_norm();
			 if (thisDistance < minDistance) {
			      minDistance = thisDistance;
			 }
		    }
	       }
		    
	       // Got the minDistance for current point. Save it in pdfPtr.
	       pdfIdx = maskIdx;

	       // reduction for pdf.
#pragma omp critical(updatepdf)
	       {
		    pdfPtr->SetPixel(pdfIdx, pdfPtr->GetPixel(pdfIdx) + minDistance * minDistance);
	       } 
	  } // in mask.
     } // maskIt.
}
		       
double MCModel::estLabelsByMeanSub(ImageType5DFloat::Pointer imagePtr,
				   ImageType3DChar::Pointer maskPtr,
				   ImageType4DFloat::Pointer allSubDistPtr,
				   unsigned int subIdx)
{
     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     ImageType3DChar::RegionType maskRegion = maskPtr->GetLargestPossibleRegion();
     IteratorType3DChar maskIt(maskPtr, maskRegion);
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();
     ImageType4DFloat::IndexType allSubDistIdx;
     allSubDistIdx[3] = subIdx;

     // images
     ImageType5DFloat::IndexType imageIdx;    
     ImageType5DFloat::RegionType imageRegion;
     imageRegion = imagePtr->GetLargestPossibleRegion();
     ImageType5DFloat::SizeType imageSize = imageRegion.GetSize();

     vnl_vector <float> timeSeries(imageSize[3], 0);
     float dist = 0;
     unsigned clsIdx = 0;

     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       maskIdx = maskIt.GetIndex();
	       imageIdx[0] = maskIdx[0];
	       imageIdx[1] = maskIdx[1];
	       imageIdx[2] = maskIdx[2];
	       imageIdx[4] = subIdx;
	       for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) { // time points.
		    timeSeries[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
	       }

	       allSubDistIdx[0] = maskIdx[0];
	       allSubDistIdx[1] = maskIdx[1];
	       allSubDistIdx[2] = maskIdx[2];
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
		    dist = (timeSeries - m_vmm[subIdx].comp[clsIdx].mu).two_norm();
		    allSubDistIdx[3] = clsIdx;
#pragma omp critical(updateAllSubDist)
		    {
			 allSubDistPtr->SetPixel(allSubDistIdx, allSubDistPtr->GetPixel(allSubDistIdx) + dist * dist);
		    }
	       } // end for clsIdx
	  } // maskIt > 0
     } // for maskIt
}

int MCModel::SetVMM(std::vector <VMM> & vmm)
{
     m_vmm = vmm;
     return 0;
}

int MCModel::GetVMM(std::vector <VMM> & vmm)
{
     vmm = m_vmm;
     return 0;
}
