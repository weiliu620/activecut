#include "common_gmmodel.h"
twister_base_gen_type mygenerator(42u);

int MergeObs(ImageType4DFloat::Pointer inPtr,      
	     VnlVectorImageType::Pointer & obsVec);

unsigned ReadObs(std::string obspath,
	    VnlVectorImageType::Pointer & obsPtr,
	    ImageType3DChar::Pointer maskPtr);

double kmeans(VnlVectorImageType::Pointer obsPtr,
	      ImageType3DChar::Pointer maskPtr,
	      ImageType3DChar::Pointer groupPtr,
	      std::vector< vnl_vector<float> > & cc,
	      float converge_eps);

int kmeanspp(VnlVectorImageType::Pointer obsPtr,
	     ImageType3DChar::Pointer maskPtr,
	     std::vector< vnl_vector<float> > & cc);

int minDist(VnlVectorImageType::Pointer obsPtr,
	    ImageType3DChar::Pointer maskPtr,
	    ImageType3DFloat::Pointer pdfPtr,
	    std::vector< vnl_vector<float> > & cc);

int EstimateMu(VnlVectorImageType::Pointer obsPtr,
	       ImageType3DChar::Pointer maskPtr,
	       ImageType3DChar::Pointer groupPtr,
	       std::vector< vnl_vector<float> > & cc);

double EstLabelsByMean(VnlVectorImageType::Pointer obsPtr,
		       ImageType3DChar::Pointer maskPtr,
		       ImageType3DChar::Pointer groupPtr,
		       std::vector< vnl_vector<float> > & cc);

int PrintClusterCenter( std::vector< vnl_vector<float> > & cc);

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     unsigned int kmeansIter = 1;
     int numClusters = 4;
     std::string obspath, outGrpLabel, maskimage;
     unsigned seed = 0;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Kmeans clustering on multiple subject images.")
	  ("numclusters,k", po::value<int>(&numClusters)->default_value(6),
	       "Number of clusters. Default is 6.")
	  ("kmeansiter,i", po::value<unsigned int>(&kmeansIter)->default_value(20),
	   "iteration # for kmeans. ")
	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("mask,m", po::value<std::string>(&maskimage)->default_value("maskimage.nii"), "input maskimage file")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")
	  ("obspath,f", po::value<std::string>(&obspath)->default_value("."),
	   "noised image file")
	  ("grouplabel,g", po::value<std::string>(&outGrpLabel)->default_value("outgrplabel.nii"), "output group label file. nii format.");

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
     unsigned numSubs;

     // read in mask image
     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(maskimage);
     maskReader->Update();
     ImageType3DChar::Pointer maskPtr = maskReader->GetOutput();

     ImageType3DChar::RegionType maskRegion;
     maskRegion = maskReader->GetOutput()->GetLargestPossibleRegion();
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();

     ImageType3DChar::PointType maskOrigin = maskPtr->GetOrigin();
     ImageType3DChar::SpacingType maskSpacing = maskPtr->GetSpacing();
     ImageType3DChar::DirectionType maskDirection = maskPtr->GetDirection();

	   
     // Allocate memory for saving vectorized images.
     VnlVectorImageType::Pointer  obsPtr = VnlVectorImageType::New();
     numSubs = ReadObs(obspath, obsPtr, maskPtr);
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();

     if (maskSize[0] != obsSize[0] 
     	 || maskSize[1] != obsSize[1] 
     	 || maskSize[2] != obsSize[2]) {
	  std::cout << "mask image and true label image have different size. Need double check before masking. Exit. " << std::endl;
	  exit(1);
     }


     // Allocate memory for group level label map.
     ImageType3DChar::Pointer grpPtr = ImageType3DChar::New();
     ImageType3DChar::IndexType grpIdx;
     grpIdx.Fill(0);

     ImageType3DChar::SizeType grpSize;
     grpSize = maskSize;

     ImageType3DChar::RegionType grpRegion;
     grpRegion.SetSize(grpSize);
     grpRegion.SetIndex(grpIdx);
     grpPtr->SetRegions(grpRegion);
     grpPtr->Allocate();
     grpPtr->FillBuffer(-1);

     grpPtr->SetOrigin( grpPtr->GetOrigin() );
     grpPtr->SetSpacing( grpPtr->GetSpacing() ); 
     grpPtr->SetDirection( grpPtr->GetDirection () );

     // define a cc and bestcc.
     std::vector< vnl_vector<float> > cc(numClusters);
     std::vector< vnl_vector<float> > bestcc(numClusters);
     unsigned clsIdx = 0;
     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  cc[clsIdx].set_size(numSubs);
	  bestcc[clsIdx].set_size(numSubs);
     }
     
     unsigned repeatIdx = 0; // repeat index of K-means.
     double minMeanSSR = 1e10, meanSSR = 0;

     for (repeatIdx = 0; repeatIdx < kmeansIter; repeatIdx ++) {
     	  printf("kmeans run  %i begin:\n", repeatIdx+1);
     	  meanSSR = kmeans(obsPtr, grpPtr, maskPtr, cc, 1e-6);

	  if (verbose >= 1) {
	       printf("minMeanSSR = %f, current MeanSSR = %f.\n", minMeanSSR, meanSSR);
	  }

     	  if (meanSSR < minMeanSSR) {
     	       minMeanSSR = meanSSR;
     	       // Save best mu in to bestvmm
	       bestcc = cc;
     	  }
     }

     // Found the bestcc. Restore the bestcc to cc.
     cc = bestcc;

     // Given bestcc, estimate labels.
     EstLabelsByMean(obsPtr, grpPtr, maskPtr, cc);     
     save3dchar(grpPtr, outGrpLabel);
}

double kmeans(VnlVectorImageType::Pointer obsPtr,
	      ImageType3DChar::Pointer groupPtr,
	      ImageType3DChar::Pointer maskPtr,
	      std::vector< vnl_vector<float> > & cc,
	      float converge_eps)
{
     unsigned subIdx = 0, clsIdx = 0;
     unsigned numClusters = cc.size();
     unsigned numSubs = cc[0].size();
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();

     // Init kmeans by kmeans++ method.
     kmeanspp(obsPtr, maskPtr, cc);
     printf("after kmans++, cluster centers:\n");
     PrintClusterCenter(cc);

     // Define a prevcc to save previous cc.
     std::vector< vnl_vector<float> > prevcc(numClusters);

     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  cc[clsIdx].set_size(numSubs);
	  prevcc[clsIdx].set_size(numSubs);
     }

     double meanSSR = 0; 
     bool converge = 1;
     do {
     	  // update label assignment.
     	  meanSSR = EstLabelsByMean(obsPtr, maskPtr, groupPtr, cc);

	  
     	  // estimate mu. Save vmm before doing that.
	  prevcc = cc;
     	  EstimateMu(obsPtr, maskPtr, groupPtr, cc);

	  printf("kmeans(): afterEstimateMu().\n");
	  PrintClusterCenter(cc);
	  converge = 1;
	  for (clsIdx = 0; clsIdx < numClusters; clsIdx++) { 
	       if ((cc[clsIdx] - prevcc[clsIdx]).two_norm() > converge_eps) {
		    converge = 0;
		    break;
	       }
	  } // clsIdx
     }
     while(!converge);
     
     return meanSSR;
}

int kmeanspp(VnlVectorImageType::Pointer obsPtr,
	     ImageType3DChar::Pointer maskPtr,
	     std::vector< vnl_vector<float> > & cc)
{

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);


     unsigned int clsIdx = 0;
     float randNum = 0;
     float cdfValue = 0;
     float sumOfPdf = 0;
     unsigned subIdx = 0;
     unsigned numClusters = cc.size();
     unsigned numSubs = cc[0].size();

     // images.
     VnlVectorImageType::IndexType obsIdx;    
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DFloat::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     
     // Allocate memory for pdf of each data point.
     ImageType3DFloat::Pointer pdfPtr = ImageType3DFloat::New();
     ImageType3DFloat::IndexType pdfStart;
     ImageType3DFloat::IndexType pdfIdx;
     pdfStart.Fill(0);
     ImageType3DFloat::SizeType pdfSize = obsSize;

     ImageType3DFloat::RegionType pdfRegion;
     pdfRegion.SetSize(pdfSize);
     pdfRegion.SetIndex(pdfStart);
     pdfPtr->SetRegions(pdfRegion);
     pdfPtr->Allocate();
     pdfPtr->FillBuffer(0);

     vnl_vector <float> obsVector(numSubs, 0);

     // Find a random bogus center.
     do {
	  maskIdx[0] = floor(uni() * maskSize[0]);
	  maskIdx[1] = floor(uni() * maskSize[1]);
	  maskIdx[2] = floor(uni() * maskSize[2]);
     }
     while(maskPtr->GetPixel(maskIdx) <= 0);
     // OK I find the point, and assign it to mu_0 (for all sub)
     obsIdx = maskIdx;
     cc[0] =  obsPtr->GetPixel(obsIdx);

     // Look for all cluster center by kmeans++. This also
     // includes first cluster. At last we throw away the bogus
     // center.
     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  // assume we have 1:clsIdx-1 mu available.
	  minDist(obsPtr, maskPtr, pdfPtr, cc);
	       
	  // compute sum of pdf.
	  sumOfPdf = 0;
	  for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	       if (maskIt.Get() > 0) {
		    maskIdx = maskIt.GetIndex();
		    pdfIdx = maskIdx;
		    sumOfPdf = sumOfPdf + pdfPtr->GetPixel(pdfIdx);
	       }
	  }

	  // normalize pdf.
	  for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	       if (maskIt.Get() > 0) {
		    maskIdx = maskIt.GetIndex();
		    pdfIdx = maskIdx;
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
	  obsIdx  = maskIdx;
	  cc[clsIdx] = obsPtr->GetPixel(obsIdx);
     }

}

int minDist(VnlVectorImageType::Pointer obsPtr,
	    ImageType3DChar::Pointer maskPtr,
	    ImageType3DFloat::Pointer pdfPtr,
	     std::vector< vnl_vector<float> > & cc)
{
     unsigned numClusters = cc.size();
     unsigned numSubs = cc[0].size();

     // images.
     VnlVectorImageType::IndexType obsIdx;    
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();
     IteratorTypeVnlVector obsIt(obsPtr, obsPtr->GetLargestPossibleRegion() );

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DFloat::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();

     ImageType3DFloat::IndexType pdfIdx;          
     IteratorType3DFloat pdfIt(pdfPtr, pdfPtr->GetLargestPossibleRegion() );
     double minDistance = 0;
     double thisDistance = 0;
     int prevIdx = 0;
     vnl_vector <float> obsVector(numSubs, 0);
     unsigned clsIdx = 0;

     for (pdfIt.GoToBegin(), obsIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++obsIt, ++pdfIt) {
	  if (maskIt.Get() > 0) {
	       obsVector = obsIt.Get();
	       minDistance = 1e10;
	       if (clsIdx == 0) {
		    minDistance = (obsVector - cc[0]).squared_magnitude();
	       }
	       else {
		    for (prevIdx = 0; prevIdx < clsIdx; prevIdx ++) {
			 thisDistance = (obsVector - cc[prevIdx]).two_norm();
			 if (thisDistance < minDistance) {
			      minDistance = thisDistance;
			 }
		    }
	       }
		    
	       // Got the minDistance for current point. Save it in pdfPtr.
	       pdfIt.Set(minDistance);
	  } // in mask.
     } // maskIt.     
}


double EstLabelsByMean(VnlVectorImageType::Pointer obsPtr,
		       ImageType3DChar::Pointer maskPtr,
		       ImageType3DChar::Pointer grpPtr,
		       std::vector< vnl_vector<float> > & cc)

{
     unsigned int clsIdx = 0;
     unsigned int subIdx = 0;
     unsigned numClusters = cc.size();
     unsigned numSubs = cc[0].size();

     // images.
     VnlVectorImageType::IndexType obsIdx;    
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();
     IteratorTypeVnlVector obsIt(obsPtr, obsPtr->GetLargestPossibleRegion() );

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DFloat::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     

     // group level label map.
     ImageType3DChar::IndexType grpIdx;
     ImageType3DChar::RegionType grpRegion;
     ImageType3DChar::SizeType grpSize = grpPtr->GetLargestPossibleRegion().GetSize();
     grpIdx.Fill( 0 );
     IteratorType3DChar grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );

     unsigned nearestClsLabel = 0;
     double nearestDistance = 1000000;
     float thisDistance = 0;

     double meanSSR = 0;
     unsigned long numAllPoints = 0;


     vnl_vector <float> obsVector(numSubs, 0);

     for (grpIt.GoToBegin(), maskIt.GoToBegin(), obsIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++ obsIt, ++grpIt) {
	  if (maskIt.Get() > 0) {	  
	       numAllPoints ++;
	       maskIdx = maskIt.GetIndex();
	       obsVector = obsIt.Get();

	       nearestDistance = 1000000;
	       for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
		    thisDistance = (obsVector - cc[clsIdx]).two_norm();
		    if (thisDistance < nearestDistance) {
			 nearestDistance = thisDistance;
			 nearestClsLabel = clsIdx;
		    }
	       }
	       // Found the nearest cluster label.
	       grpIt.Set( nearestClsLabel );
	       meanSSR += nearestDistance;
	  } // maskIt > 0
     } // for
     meanSSR = meanSSR / numAllPoints;
     return meanSSR;
}

int EstimateMu(VnlVectorImageType::Pointer obsPtr,
	       ImageType3DChar::Pointer maskPtr,
	       ImageType3DChar::Pointer grpPtr,
	       std::vector< vnl_vector<float> > & cc)

{
     unsigned int clsIdx = 0;
     unsigned int subIdx = 0;
     unsigned numClusters = cc.size();
     unsigned numSubs = cc[0].size();

     // images.
     VnlVectorImageType::IndexType obsIdx;    
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();
     IteratorTypeVnlVector obsIt(obsPtr, obsPtr->GetLargestPossibleRegion() );

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DFloat::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     

     // group level label map.
     ImageType3DChar::IndexType grpIdx;
     ImageType3DChar::RegionType grpRegion;
     ImageType3DChar::SizeType grpSize = grpPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DChar grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );


     vnl_vector<unsigned> numPoints(numClusters, 0);

     // clear to zero for all cc.
     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  cc[clsIdx].fill (0);
     }
     for (grpIt.GoToBegin(), maskIt.GoToBegin(), obsIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++ obsIt, ++grpIt) {
	  if (maskIt.Get() > 0) {	  
	       clsIdx = grpIt.Get();
	       cc[clsIdx] += obsIt.Get();
	       numPoints[clsIdx] ++;
	  }
     }

     // devide by total number of pts in this component.
     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  cc[clsIdx] = cc[clsIdx] / numPoints[clsIdx];
     }
     
}

int PrintClusterCenter( std::vector< vnl_vector<float> > & cc)
{
     unsigned numClusters = cc.size();
     unsigned numSubs = cc[0].size();

     for (unsigned clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  printf("cls[%i]: ", clsIdx);
	  for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	       printf("%2.2f ", cc[clsIdx][subIdx]);
	  }
	  printf("\n");
     }
     printf("\n");
}
