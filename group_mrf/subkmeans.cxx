#include "commonalt.h"
#include "MCModel_init.h"
#include <utilalt.h>
struct Parm
{
     unsigned numClusters;
     float converge_eps;
     std::vector<vnl_vector<float> > mu;
     std::vector<unsigned> numPoints;
     unsigned verbose;
};

double kmeans(ImageType4DFloat::Pointer fmriPtr, 
	      ImageType3DChar::Pointer maskPtr,
	      ImageType3DChar::Pointer outPtr,
	      Parm & parm);

int kmeanspp(ImageType4DFloat::Pointer fmriPtr,
	     ImageType3DChar::Pointer maskPtr,
	     Parm & parm);

int ComputeDistance(ImageType3DFloat::Pointer pdfPtr, 
		    ImageType4DFloat::Pointer imagePtr,
		    ImageType3DChar::Pointer maskPtr,
		    Parm & parm,
		    unsigned int clsIdx);

int EstimateMu(ImageType3DChar::Pointer outPtr,
	       ImageType4DFloat::Pointer imagePtr,
	       ImageType3DChar::Pointer maskPtr,
	       Parm & parm);

double EstimateLabels(ImageType4DFloat::Pointer imagePtr,
		      ImageType3DChar::Pointer maskPtr,
		      Parm & parm,
		      ImageType3DChar::Pointer outPtr);


twister_base_gen_type mygenerator(42u);

int main(int argc, char* argv[])
{
     unsigned short kmeansIter = 0;
     unsigned seed = 0;
     int numClusters = 4;
     unsigned short verbose = 0;
     std::string fmrifile, labelmapfile, maskimage;
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "kmeans clustering for single subject fmri image.")
	  ("kmeansiter", po::value<unsigned short>(&kmeansIter)->default_value(20),
	   "iteration # for kmeans. ")
	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("numClusters,k", po::value<int>(&numClusters)->default_value(4),
	       "Number of labels. Default is 5.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")

	  ("fmri,f", po::value<std::string>(&fmrifile)->default_value("fmri.nii.gz"),
	   "fmri image file")
	  ("labelmap,o", po::value<std::string>(&labelmapfile)->default_value("labelmap.nii.gz"),
	   "output labeled lebel image.")
	  ("mask,m", po::value<std::string>(&maskimage)->default_value("maskimage.nii"), "input maskimage file");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: generateimage [options]\n";
	       std::cout << mydesc << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // init random generator.
     mygenerator.seed(static_cast<unsigned int>(seed));

     // read in mask image
     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(maskimage);
     maskReader->Update();
     ImageType3DChar::Pointer maskPtr = maskReader->GetOutput();
     ImageType3DChar::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();

     // Allocate memory for group level label map.
     ImageType3DChar::Pointer outPtr = ImageType3DChar::New();
     ImageType3DChar::IndexType outIdx;
     outIdx.Fill(0);

     ImageType3DChar::RegionType outRegion;
     outRegion.SetSize(maskSize);
     outRegion.SetIndex(outIdx);
     outPtr->SetRegions(outRegion);
     outPtr->Allocate();
     outPtr->FillBuffer(-1);

     outPtr->SetOrigin( maskPtr->GetOrigin() );
     outPtr->SetSpacing( maskPtr->GetSpacing() ); 
     outPtr->SetDirection( maskPtr->GetDirection () );

 
     // read in fmri image
     ReaderType4DFloat::Pointer fmriReader = ReaderType4DFloat::New();
     fmriReader->SetFileName(fmrifile);
     fmriReader->Update();
     ImageType4DFloat::Pointer fmriPtr = fmriReader->GetOutput();
     ImageType4DFloat::SizeType fmriSize = fmriPtr->GetLargestPossibleRegion().GetSize();

     // Define parameters, including mu.
     Parm parm;
     parm.numClusters = numClusters;
     parm.converge_eps = 1e-4;
     parm.mu.resize(parm.numClusters);
     parm.numPoints.resize(parm.numClusters);
     parm.verbose = verbose;

     std::vector<vnl_vector<float> > bestmu(numClusters);
     unsigned clsIdx = 0;
     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  parm.mu[clsIdx].set_size(fmriSize[3]);
	  parm.mu[clsIdx].fill( 0 );
	  bestmu[clsIdx].set_size(fmriSize[3]);
	  bestmu[clsIdx].fill( 0 );
     }


     unsigned repeatIdx = 0; // repeat index of K-means.
     double minMeanSSR = 1e10, meanSSR = 0;

     for (repeatIdx = 0; repeatIdx < kmeansIter; repeatIdx ++) {
     	  printf("kmeans for %s, run  %i begin:\n", fmrifile.c_str(), repeatIdx+1);

     	  meanSSR = kmeans(fmriPtr, maskPtr, outPtr, parm);

	  if (verbose >= 1) {
	       printf("minMeanSSR = %f, current MeanSSR = %f.\n", minMeanSSR, meanSSR);
	  }

     	  if (meanSSR < minMeanSSR) {
     	       minMeanSSR = meanSSR;
     	       // Save best mu in to bestmu
	       bestmu = parm.mu;
     	  }
     }

     // Found the best mu. Restore the best mu into parm.
     parm.mu = bestmu;

     // Given best mu, estimate labels. 
     EstimateLabels(fmriPtr, maskPtr, parm, outPtr);

     printf("subkmeans(): %s, Final minMeanSSR = %f.\n", fmrifile.c_str(), minMeanSSR);
     save3dcharInc(outPtr, labelmapfile);
}

double kmeans(ImageType4DFloat::Pointer fmriPtr, 
	      ImageType3DChar::Pointer maskPtr,
	      ImageType3DChar::Pointer outPtr,
	      Parm & parm)
{
     float meanSSR = 0, meanSSR_old = 0;
     ImageType4DFloat::SizeType fmriSize = fmriPtr->GetLargestPossibleRegion().GetSize();
     // Init kmeans by kmeans++ method.
     kmeanspp(fmriPtr, maskPtr, parm);

     do {
     	  // update label assignment.
	  meanSSR_old = meanSSR;
     	  meanSSR = EstimateLabels(fmriPtr, maskPtr, parm, outPtr);

     	  EstimateMu(outPtr, fmriPtr, maskPtr, parm);
	  if (parm.verbose >= 1) {
	       printf("    kemans(), meanSSR = %.4f, meanSSR_old = %.4f\n", meanSSR, meanSSR_old);
	  }
     } while (fabs((meanSSR - meanSSR_old) / meanSSR_old) > 1e-4);

     return meanSSR;
}

int kmeanspp(ImageType4DFloat::Pointer fmriPtr,
	     ImageType3DChar::Pointer maskPtr,
	     Parm & parm)
{
     float randNum = 0, cdfValue = 0;
     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     // fmris
     ImageType4DFloat::IndexType fmriIdx;    
     ImageType4DFloat::SizeType fmriSize = fmriPtr->GetLargestPossibleRegion().GetSize();

     // mask image.
     ImageType3DFloat::IndexType maskIdx;     
     ImageType3DFloat::RegionType maskRegion = maskPtr->GetLargestPossibleRegion();
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DFloat::SizeType maskSize = maskRegion.GetSize();

     // Allocate memory for pdf of each data point.
     ImageType3DFloat::Pointer pdfPtr = ImageType3DFloat::New();
     ImageType3DFloat::IndexType pdfStart;
     ImageType3DFloat::IndexType pdfIdx;
     pdfStart.Fill(0);
     ImageType3DFloat::SizeType pdfSize;
     pdfSize[0] = fmriSize[0];
     pdfSize[1] = fmriSize[1];
     pdfSize[2] = fmriSize[2];

     ImageType3DFloat::RegionType pdfRegion;
     pdfRegion.SetSize(pdfSize);
     pdfRegion.SetIndex(pdfStart);
     pdfPtr->SetRegions(pdfRegion);
     pdfPtr->Allocate();
     pdfPtr->FillBuffer(0);

     vnl_vector <float> timeSeries(fmriSize[3], 0);

     // Find a random bogus center.
     do {
	  maskIdx[0] = floor(uni() * maskSize[0]);
	  maskIdx[1] = floor(uni() * maskSize[1]);
	  maskIdx[2] = floor(uni() * maskSize[2]);
     }
     while(maskPtr->GetPixel(maskIdx) <= 0);

     fmriIdx[0] = maskIdx[0];
     fmriIdx[1] = maskIdx[1];
     fmriIdx[2] = maskIdx[2];

     // OK I find the point, and assign it to mu_0 (for all sub)

     for (fmriIdx[3] = 0; fmriIdx[3] < fmriSize[3]; fmriIdx[3] ++) {
	  parm.mu[0][fmriIdx[3]] = fmriPtr->GetPixel(fmriIdx);
     }
     unsigned clsIdx = 0;
     for (clsIdx = 0; clsIdx < parm.numClusters; clsIdx ++) {
	  // compute the distance to exisint cluster centers.
	  ComputeDistance(pdfPtr, fmriPtr, maskPtr, parm, clsIdx);
	  // compute sum of pdf.

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

	  if (parm.verbose >= 1) {
	       printf("kmeanspp(): cls[%i]: maskIdx[%i][%i][%i]. \n", clsIdx, maskIdx[0], maskIdx[1], maskIdx[2]);
	  }

	  // found the data point. Assign it to current cluster for all subs.
	  fmriIdx[0] = maskIdx[0];
	  fmriIdx[1] = maskIdx[1];
	  fmriIdx[2] = maskIdx[2];

	  for (fmriIdx[3] = 0; fmriIdx[3] < fmriSize[3]; fmriIdx[3] ++) {
	       parm.mu[clsIdx][fmriIdx[3]] = fmriPtr->GetPixel(fmriIdx);
	  }
     } // clsIdx

}

int ComputeDistance(ImageType3DFloat::Pointer pdfPtr, 
		    ImageType4DFloat::Pointer imagePtr,
		    ImageType3DChar::Pointer maskPtr,
		    Parm & parm,
		    unsigned int clsIdx)
{
     // images
     ImageType4DFloat::IndexType imageIdx;    
     ImageType4DFloat::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     ImageType3DChar::RegionType maskRegion = maskPtr->GetLargestPossibleRegion();
     IteratorType3DChar maskIt(maskPtr, maskRegion);
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();

     ImageType3DFloat::IndexType pdfIdx;          
     IteratorType3DFloat pdfIt(pdfPtr, pdfPtr->GetLargestPossibleRegion().GetSize());
     vnl_vector <float> timeSeries(imageSize[3], 0);

     double minDistance = 0;
     double thisDistance = 0;
     int prevIdx = 0;

     float sumPDF = 0;
     for (maskIt.GoToBegin(), pdfIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++ pdfIt) {
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
		    minDistance = (timeSeries - parm.mu[0]).two_norm();
	       }
	       else {
		    for (prevIdx = 0; prevIdx < clsIdx; prevIdx ++) {
			 thisDistance = (timeSeries - parm.mu[clsIdx]).two_norm();
			 if (thisDistance < minDistance) {
			      minDistance = thisDistance;
			 }
		    }
	       }
	       // Got the minDistance for current point. Save it in pdfPtr.
	       pdfIt.Set(minDistance * minDistance);
	       sumPDF += minDistance * minDistance;
	  } // in mask
     } // maskIT
     
     // normalize the pdf.
     for (pdfIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++ pdfIt) {
	  if (maskIt.Get() > 0) {
	       pdfIt.Set( pdfIt.Get() / sumPDF );
	  } // in mask
     } // maskIT

}


double EstimateLabels(ImageType4DFloat::Pointer imagePtr,
		       ImageType3DChar::Pointer maskPtr,
		       Parm & parm,
		       ImageType3DChar::Pointer outPtr)
{
     unsigned int clsIdx = 0;

     // images
     ImageType4DFloat::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();
     ImageType4DFloat::IndexType imageIdx;    

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     ImageType3DChar::RegionType maskRegion = maskPtr->GetLargestPossibleRegion();
     IteratorType3DChar maskIt(maskPtr, maskRegion);
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();

     // output label map.
     ImageType3DChar::IndexType outIdx;
     IteratorType3DChar outIt(outPtr, outPtr->GetLargestPossibleRegion());

     ImageType3DChar::SizeType outSize = outPtr->GetLargestPossibleRegion().GetSize();
     outIdx[0] = 0;
     outIdx[1] = 0;
     outIdx[2] = 0;

     unsigned nearestClsLabel = 0;
     double nearestDistance = 1000000;

     double meanSSR = 0;
     unsigned long numAllPoints = 0;

     vnl_vector <float> timeSeries(imageSize[3], 0);

     float thisDist = 0;
     for (maskIt.GoToBegin(), outIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++ outIt) {
	  if (maskIt.Get() > 0) {	  
	       numAllPoints ++;
	       maskIdx = maskIt.GetIndex();
	       imageIdx[0] = maskIdx[0];
	       imageIdx[1] = maskIdx[1];
	       imageIdx[2] = maskIdx[2];	  
	       for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
		    timeSeries[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
	       }

	       nearestDistance = 1000000;
	       for (clsIdx = 0; clsIdx < parm.numClusters; clsIdx++) {
		    thisDist = (timeSeries - parm.mu[clsIdx]).two_norm();
		    if (thisDist < nearestDistance) {
			 nearestClsLabel = clsIdx;
			 nearestDistance = thisDist;
		    }
	       }

	       // Found the nearest cluster label.
	       outIt.Set(nearestClsLabel);
	       meanSSR += nearestDistance;
	  }
     }
     meanSSR = meanSSR / numAllPoints;
     return meanSSR;
}

     
int EstimateMu(ImageType3DChar::Pointer outPtr,
	       ImageType4DFloat::Pointer imagePtr,
	       ImageType3DChar::Pointer maskPtr,
	       Parm & parm)

{
     unsigned int clsIdx = 0;

     // images
     ImageType4DFloat::IndexType imageIdx;    
     ImageType4DFloat::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     ImageType3DChar::RegionType maskRegion = maskPtr->GetLargestPossibleRegion();
     IteratorType3DChar maskIt(maskPtr, maskRegion);
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();

     // outs
     ImageType3DChar::IndexType outIdx;
     IteratorType3DChar outIt(outPtr, outPtr->GetLargestPossibleRegion());

     vnl_vector <float> timeSeries(imageSize[3], 0);

     // reset all mu to zero.
     for (clsIdx = 0; clsIdx < parm.numClusters; clsIdx ++) {
	  parm.mu[clsIdx] = 0;
	  parm.numPoints[clsIdx] = 0;
     }

     for (maskIt.GoToBegin(), outIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++ outIt) {
	  if (maskIt.Get() > 0) {
	       maskIdx = maskIt.GetIndex();

	       // First save time series at this position into a vector.
	       imageIdx[0] = maskIdx[0];
	       imageIdx[1] = maskIdx[1];
	       imageIdx[2] = maskIdx[2];
	       for (imageIdx[3] = 0; imageIdx[3] < imageSize[3]; imageIdx[3] ++) {
		    timeSeries[imageIdx[3]] = imagePtr->GetPixel(imageIdx);
	       }

	       clsIdx = outIt.Get();

	       parm.mu[clsIdx] += timeSeries;
	       parm.numPoints[clsIdx] ++;
	  }  // maskIt > 0
     } // for maskIt

     // compute mean time series and meanNorm.
     float meanNorm = 0;
     for (clsIdx = 0; clsIdx < parm.numClusters; clsIdx ++) {
	  // meanNorm see page 1350 of "Clustering on the Unit
	  // Hypersphere using von Mises-Fisher distributions" by
	  // Banerjee.
	  parm.mu[clsIdx] = parm.mu[clsIdx].normalize();
     }
     return 0;
}
