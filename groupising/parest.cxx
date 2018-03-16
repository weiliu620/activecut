#include "common_mcem.h"


typedef itk::ConstantBoundaryCondition< ImageType3DChar >  MyBoundCondType;
typedef itk::NeighborhoodIterator< ImageType3DChar, MyBoundCondType > MyNeiItType;

// model parameters.
struct ParType 
{
     unsigned numClusters;
     unsigned numSubs;
     float alpha;
     float betag;
     float betaz;
     unsigned short verbose;
};
	  
unsigned ReadSubToVec(std::string subdir,
		      std::vector<ImageType3DChar::Pointer> & subVec,
		      ImageType3DChar::Pointer maskPtr);

float GrpBetaDrv(ImageType3DChar::Pointer & grpPtr,
	       std::vector<ImageType3DChar::Pointer>  & subVec,
	       ImageType3DVecUS::Pointer histPtr,
	       ImageType3DChar::Pointer maskPtr,
	       double & drv1, double & drv2,
	       ParType & par);

float SubBetaDrv(ImageType3DChar::Pointer &  grpPtr,
	       std::vector<ImageType3DChar::Pointer>  & subVec,
	       ImageType3DChar::Pointer maskPtr,
	       double & drv1, double & drv2,
	       ParType & par);

float GrpAlphaDrv(ImageType3DChar::Pointer & grpPtr,
	       std::vector<ImageType3DChar::Pointer>  & subVec,
	       ImageType3DVecUS::Pointer histPtr,
	       ImageType3DChar::Pointer maskPtr,
	       double & drv1, double & drv2,
		ParType & par);

float SubAlphaDrv(ImageType3DChar::Pointer &  grpPtr,
	       std::vector<ImageType3DChar::Pointer>  & subVec,
	       ImageType3DChar::Pointer maskPtr,
	       double & drv1, double & drv2,
		ParType & par);

int main (int argc, char* argv[])
{
     ParType par;
     float tol;
     std::string trueimage, maskimage;
     std::string grpLabel, subdir;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "create group label map.")
	  ("betag", po::value<float>(&par.betag)->default_value(0.5),
	   "Initial value of group level pairwise interation term")

	  ("betaz", po::value<float>(&par.betaz)->default_value(0.5),
	   "Initial value of subject level pairwise interation term")

	  ("alpha", po::value<float>(&par.alpha)->default_value(0.5),
	   "Initial value of connection between group lebel and individual subject level.")

	  ("numClusters,k", po::value<unsigned>(&par.numClusters)->default_value(6),
	       "Number of clusters.")

	  ("grouplabel,g", po::value<std::string>(&grpLabel)->default_value("outgrplabel.nii"), "group label file. nii format.")

	  ("subdir,s", po::value<std::string>(&subdir)->default_value("subdir"), "subject label maps directory")

	  ("tol", po::value<float>(&tol)->default_value(1e-4),
	   "Change of joint log likelihood before the algorithm stops.")

	  ("verbose,v", po::value<unsigned short>(&par.verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");



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

     // Read group label map as mask. Decrease 1 and get group label map for
     // internal use.
     ReaderType3DChar::Pointer grpReader = ReaderType3DChar::New();
     grpReader->SetFileName(grpLabel);
     ImageType3DChar::Pointer maskPtr = grpReader->GetOutput();
     maskPtr->Update();


     AddConstantToImageFilterType::Pointer myfilter = AddConstantToImageFilterType::New();

     // keep original input image unchanged.
     myfilter->InPlaceOff();
     myfilter->SetInput( grpReader->GetOutput());
     myfilter->SetConstant(-1);

     ImageType3DChar::Pointer grpPtr = myfilter->GetOutput();
     grpPtr->Update();


     // Read subjects images.
     std::vector<ImageType3DChar::Pointer>  subVec;
     par.numSubs = ReadSubToVec(subdir, subVec, maskPtr);


     // create histogram data block.
     itk::VariableLengthVector<unsigned short> histVector(par.numClusters);

     ImageType3DVecUS::Pointer histPtr = ImageType3DVecUS::New();
     ImageType3DVecUS::RegionType histRegion = grpPtr->GetLargestPossibleRegion();
     histPtr->SetRegions(histRegion);
     histPtr->SetNumberOfComponentsPerPixel(par.numClusters);
     histPtr->SetVectorLength(par.numClusters);
     histPtr->Allocate();

     histVector.Fill(0);
     histPtr->FillBuffer(histVector);
     IteratorType3DVecUS histIt(histPtr, histPtr->GetLargestPossibleRegion() );

     // compute histogram
     for (unsigned short subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  ConstIteratorType3DChar sampleIt(subVec[subIdx], subVec[subIdx]->GetLargestPossibleRegion() );

	  for (sampleIt.GoToBegin(), histIt.GoToBegin(); 
	       !histIt.IsAtEnd(); ++ histIt, ++sampleIt) {
	       if (sampleIt.Get() >= 0) {
		    histVector = histIt.Get();
		    histVector[sampleIt.Get()] ++;
		    histIt.Set(histVector);
	       }
	  } // for sampelIt.
     } // for subIdx



     double beta_g_old = par.betag, beta_z_old = par.betaz, alpha_old = par.alpha;
     double drv1 = 0, drv2 = 0, Q1 = 0, Q2 = 0, old_Q1 = 0, old_Q2 = 0;
     double drv1g = 0, drv1z = 0, drv2g = 0, drv2z = 0;

     

     // compute Q1 and Q2.

     old_Q2 = SubBetaDrv(grpPtr, subVec, maskPtr, drv1, drv2, par);
     do {

	  ////// Debug ///////
	  if (par.verbose >= 4) {
	       beta_g_old = par.betag;
	       for (par.betag = beta_g_old - 0.1; par.betag < beta_g_old + 0.1; par.betag += 0.01) {
		    Q1 = GrpBetaDrv(grpPtr, subVec, histPtr, maskPtr, drv1, drv2, par);
		    printf("        beta_g: %5.4f,  Q1: %.3f,  drv1 = %.2f, drv2 = %.2f\n", par.betag, Q1, drv1, drv2);
	       }
	       par.betag = beta_g_old;
	  }
	  ///// end of debug ////

	  // beta_g
	  old_Q1 = GrpBetaDrv(grpPtr, subVec, histPtr, maskPtr, drv1, drv2, par);
	  beta_g_old = par.betag;
	  par.betag = beta_g_old - drv1/drv2;
	  
	  // update Q1
	  Q1 = GrpBetaDrv(grpPtr, subVec, histPtr, maskPtr, drv1, drv2, par);

	  if (par.verbose >= 2) {
	       printf("    mbeta_g: %5.4f -> %5.4f. Q1: %.3f -> %.3f. drv1 = %.4f, drv2 = %.4f\n", beta_g_old, par.betag, old_Q1, Q1, drv1, drv2);
	  }

	  // beta_z
	  old_Q2 = SubBetaDrv(grpPtr, subVec, maskPtr, drv1, drv2, par);
	  beta_z_old = par.betaz;
	  par.betaz = beta_z_old - drv1/drv2;
	  
	  if (par.betaz < 0) {
	       par.betaz = (beta_z_old)/2;
	  }

	  // update Q2
	  Q2 = SubBetaDrv(grpPtr, subVec, maskPtr, drv1, drv2, par);

	  if (par.verbose >= 2) {
	       printf("    mbeta_z: %5.4f -> %5.4f. Q2: %.3f -> %.3f. drv1 = %.4f, drv2 = %.4f\n", beta_z_old, par.betaz, old_Q2, Q2, drv1, drv2);
	  }

	  // alpha
	  old_Q1 = GrpAlphaDrv(grpPtr, subVec, histPtr, maskPtr, drv1g, drv2g, par);
	  old_Q2 = SubAlphaDrv(grpPtr, subVec, maskPtr, drv1z, drv2z, par);
	  alpha_old = par.alpha;

	  // the derivative of (q1+q2) is equal to the sume of their
	  // derivatives.
	  par.alpha = alpha_old - (drv1g + drv1z)/(drv2g + drv2z);

	  // update Q1 and Q2
	  Q1 = GrpAlphaDrv(grpPtr, subVec, histPtr, maskPtr, drv1g, drv2g, par);
	  Q2 = SubAlphaDrv(grpPtr, subVec, maskPtr, drv1z, drv2z, par);


	  if (par.verbose >= 3) {
	       printf("    alpha: %5.4f -> %5.4f. (Q1+Q2=Q): (%.3f + %.3f = %.3f) -> (%.3f + %.3f = %.3f)\n", alpha_old, par.alpha, old_Q1, old_Q2, old_Q1 + old_Q2, Q1, Q2, Q1+Q2);
	       printf("           (drv1g+z = drv1): (%.4f + %.4f = %.4f), (drv2g+z = drv2) = (%.4f + %.4f  = %.4f)\n", drv1g, drv1z, drv1g+drv1z, drv2g, drv2z, drv2g+drv2z);
	  }

	  if (par.verbose >= 1) {
	       printf("beta_g = %5.4f, beta_z = %5.4f, alpha = %5.4f, Q = %f\n", par.betag, par.betaz, par.alpha, (Q1+Q2));
	  }
     } while ((fabs(par.betag - beta_g_old) + fabs(par.betaz - beta_z_old) + fabs(par.alpha - alpha_old)) > tol );

     printf("TB: beta_g = %5.4f, beta_z = %5.4f, alpha = %5.4f, Q = %f\n", par.betag, par.betaz, par.alpha, (Q1+Q2));


     // Q1(grpPtr, maskPtr, par);


     return 0;
}


unsigned ReadSubToVec(std::string subdir,
		      std::vector<ImageType3DChar::Pointer> & subVec,
		      ImageType3DChar::Pointer maskPtr)
{
     // some of the code below refers to boost filesystem tutorial tut04.cpp.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     // default constructed iterator acts as a end iterator.
     boost::filesystem::path subdirVar(subdir);
     boost::filesystem::directory_iterator subdirEnd;

     ReaderType3DChar::Pointer subReader = ReaderType3DChar::New();

     
     // prepare for sorting directory entries.
     typedef std::vector<boost::filesystem::path> PathVec;
     PathVec  sortedEntries;
     copy(boost::filesystem::directory_iterator(subdir), boost::filesystem::directory_iterator(), std::back_inserter(sortedEntries) );
     sort(sortedEntries.begin(), sortedEntries.end() );

     // Get number of subjects.
     unsigned numSubs = 0;

     numSubs = sortedEntries.size();
     subVec.resize(numSubs);

     ImageType3DChar::Pointer singleSubPtr = ImageType3DChar::New();

     // allocate memory for subPtr.
     ImageType3DChar::IndexType subPixelIdx;
     subPixelIdx.Fill ( 0 );
     ImageType3DChar::SizeType subSize;

     // get image size.
     boost::filesystem::directory_iterator subdirIt(subdirVar);
     subReader->SetFileName( (*subdirIt).path().string() );
     subReader->Update();
     singleSubPtr = subReader->GetOutput();
     subSize = singleSubPtr->GetLargestPossibleRegion().GetSize();

     // region
     ImageType3DChar::RegionType subRegion;
     subRegion.SetSize (subSize);
     subRegion.SetIndex( subPixelIdx );


     unsigned subIdx = 0;
     for (PathVec::const_iterator subdirIt(sortedEntries.begin() ); subdirIt != sortedEntries.end(); ++ subdirIt) {

	  // Destination.
	  subVec[subIdx] = ImageType3DChar::New();
	  subVec[subIdx]->SetRegions( subRegion );
	  subVec[subIdx]->Allocate();
	  subVec[subIdx]->FillBuffer ( -1 );
	  IteratorType3DChar subIt(subVec[subIdx], subVec[subIdx]->GetLargestPossibleRegion() );

	  // source.
	  subReader->SetFileName( (*subdirIt).string() );
	  subReader->Update();
	  singleSubPtr = subReader->GetOutput();
	  IteratorType3DChar singleSubIt(singleSubPtr, singleSubPtr->GetLargestPossibleRegion() );
	  std::cout <<  "add " << (*subdirIt).string() << "\n";

	  // read in data for this subject. 1-based convert to 0-based for
	  // internal use.

	  for (maskIt.GoToBegin(), subIt.GoToBegin(), singleSubIt.GoToBegin(); !subIt.IsAtEnd(); ++ subIt, ++maskIt, ++singleSubIt){
	       if (maskIt.Get() > 0) {
		    subIt.Set( singleSubIt.Get() -1);
	       } // maskIt > 0
	  } // maskIt
	  
	  // done with read data for this sub. 
	  subIdx ++;
     }
     return numSubs;
}


float GrpBetaDrv(ImageType3DChar::Pointer & grpPtr,
	       std::vector<ImageType3DChar::Pointer>  & subVec,
	       ImageType3DVecUS::Pointer histPtr,
	       ImageType3DChar::Pointer maskPtr,
	       double & drv1, double & drv2,
	       ParType & par)
{
     // mask image.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     // Define neighborhood iterator
     MyBoundCondType constCondition;
     constCondition.SetConstant(-1);
     MyNeiItType::RadiusType radius;
     radius.Fill(1);

     MyNeiItType::OffsetType xplus = {{1,0, 0}};
     MyNeiItType::OffsetType xminus = {{-1, 0, 0}};
     MyNeiItType::OffsetType yplus = {{0, 1, 0}};
     MyNeiItType::OffsetType yminus = {{0, -1, 0}};
     MyNeiItType::OffsetType zplus = {{0, 0, 1}};
     MyNeiItType::OffsetType zminus = {{0, 0, -1}};

     double a = 0, b = 0, M0 = 0, M1 = 0, M2 = 0, bcur = 0, acur = 0;
     unsigned curLabel = 0;
     unsigned clsIdx = 0;
     drv1 = 0, drv2 = 0;
     
     MyNeiItType neiGrpIt(radius, grpPtr, grpPtr->GetLargestPossibleRegion() );
     neiGrpIt.OverrideBoundaryCondition(&constCondition);

     IteratorType3DVecUS histIt(histPtr, histPtr->GetLargestPossibleRegion() );

     float Q1 = 0;
     for (neiGrpIt.GoToBegin(), maskIt.GoToBegin(), histIt.GoToBegin(); !maskIt.IsAtEnd(); ++ neiGrpIt, ++maskIt, ++ histIt) {
	  if (maskIt.Get() > 0) {
	       curLabel = neiGrpIt.GetCenterPixel();
	       M0 = 0, M1 = 0, M2 = 0;
	       for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
		    a = (double)(histIt.Get()[clsIdx]) - par.numSubs;

		    // debug ///
		    itk::VariableLengthVector<unsigned short> histVector(par.numClusters);
		    histVector = histIt.Get();
		    /// end of debug ///
		    b = - (int(neiGrpIt.GetPixel(xminus) >= 0 && clsIdx != neiGrpIt.GetPixel(xminus) )
			   + int(neiGrpIt.GetPixel(xplus) >= 0 && clsIdx != neiGrpIt.GetPixel(xplus))
			   + int(neiGrpIt.GetPixel(yminus) >= 0 && clsIdx != neiGrpIt.GetPixel(yminus))
			   + int(neiGrpIt.GetPixel(yplus) >= 0 && clsIdx != neiGrpIt.GetPixel(yplus))
			   + int(neiGrpIt.GetPixel(zminus) >= 0 && clsIdx != neiGrpIt.GetPixel(zminus))
			   + int(neiGrpIt.GetPixel(zplus) >= 0 && clsIdx != neiGrpIt.GetPixel(zplus)) );

		    if (clsIdx == curLabel) {
			 acur = a;
			 bcur = b;
		    }
		    M0 += exp (a * par.alpha + b * par.betag);
		    M1 += b * exp(a * par.alpha + b * par.betag);
		    M2 += b * b * exp(a * par.alpha + b * par.betag);
	       } // for clsIdx
	       drv1 = drv1 - bcur + M1 / M0;
	       drv2 += (M2 * M0 - M1 * M1)/ (M0 * M0);
	       Q1 = Q1 + (-acur * par.alpha - bcur * par.betag + log(M0));
	  }  // maskIt > 0
     } // for maskIt

     return Q1;
}



float SubBetaDrv(ImageType3DChar::Pointer &  grpPtr,
	       std::vector<ImageType3DChar::Pointer>  & subVec,
	       ImageType3DChar::Pointer maskPtr,
	       double & drv1, double & drv2,
	       ParType & par)
{
     // mask image.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );


     // Define neighborhood iterator
     MyBoundCondType constCondition;
     constCondition.SetConstant(-1);
     MyNeiItType::RadiusType radius;
     radius.Fill(1);

     MyNeiItType::OffsetType xplus = {{1,0, 0}};
     MyNeiItType::OffsetType xminus = {{-1, 0, 0}};
     MyNeiItType::OffsetType yplus = {{0, 1, 0}};
     MyNeiItType::OffsetType yminus = {{0, -1, 0}};
     MyNeiItType::OffsetType zplus = {{0, 0, 1}};
     MyNeiItType::OffsetType zminus = {{0, 0, -1}};

     double a = 0, b = 0, M0 = 0, M1 = 0, M2 = 0, bcur = 0, acur = 0, Q2 = 0;
     unsigned curSubLabel = 0, grpLabel = 0;
     unsigned clsIdx = 0;
     drv1 = 0, drv2 = 0;
     
     IteratorType3DChar grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );

     for (unsigned short subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  MyNeiItType neiSampleIt( radius, subVec[subIdx], subVec[subIdx]->GetLargestPossibleRegion() );
	  neiSampleIt.OverrideBoundaryCondition(&constCondition);
	  
	  for (neiSampleIt.GoToBegin(), grpIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ neiSampleIt, ++ grpIt, ++ maskIt) {
	       if (maskIt.Get() > 0) {
		    curSubLabel = neiSampleIt.GetCenterPixel();

		    M0 = 0, M1 = 0, M2 = 0;
		    for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
			 a = - int (grpIt.Get() != clsIdx);
			 b = - (int(neiSampleIt.GetPixel(xminus) >= 0 && clsIdx != neiSampleIt.GetPixel(xminus) )
				+ int(neiSampleIt.GetPixel(xplus) >= 0 && clsIdx != neiSampleIt.GetPixel(xplus))
				+ int(neiSampleIt.GetPixel(yminus) >= 0 && clsIdx != neiSampleIt.GetPixel(yminus))
				+ int(neiSampleIt.GetPixel(yplus) >= 0 && clsIdx != neiSampleIt.GetPixel(yplus))
				+ int(neiSampleIt.GetPixel(zminus) >= 0 && clsIdx != neiSampleIt.GetPixel(zminus))
				+ int(neiSampleIt.GetPixel(zplus) >= 0 && clsIdx != neiSampleIt.GetPixel(zplus)) );
			 if (clsIdx == curSubLabel) {
			      bcur = b;
			      acur = a;
			 }
			 M0 += exp (a * par.alpha + b * par.betaz);
			 M1 += b * exp(a * par.alpha + b * par.betaz);
			 M2 += b * b * exp(a * par.alpha + b * par.betaz);
		    } // for clsIdx
		    drv1 += (- bcur + M1 / M0);
		    drv2 += (M2 * M0 - M1 * M1)/ (M0 * M0);
		    Q2 = Q2 + (- acur * par.alpha - bcur * par.betaz + log (M0) );
	       }  // maskIt > 0
	  } // for maskIt
     } // subIdx

     return Q2;
}


float GrpAlphaDrv(ImageType3DChar::Pointer & grpPtr,
	       std::vector<ImageType3DChar::Pointer>  & subVec,
	       ImageType3DVecUS::Pointer histPtr,
	       ImageType3DChar::Pointer maskPtr,
	       double & drv1, double & drv2,
	       ParType & par)
{
     // mask image.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     // Define neighborhood iterator
     MyBoundCondType constCondition;
     constCondition.SetConstant(-1);
     MyNeiItType::RadiusType radius;
     radius.Fill(1);

     MyNeiItType::OffsetType xplus = {{1,0, 0}};
     MyNeiItType::OffsetType xminus = {{-1, 0, 0}};
     MyNeiItType::OffsetType yplus = {{0, 1, 0}};
     MyNeiItType::OffsetType yminus = {{0, -1, 0}};
     MyNeiItType::OffsetType zplus = {{0, 0, 1}};
     MyNeiItType::OffsetType zminus = {{0, 0, -1}};

     double a = 0, b = 0, M0 = 0, M1 = 0, M2 = 0, bcur = 0, acur = 0;
     unsigned curLabel = 0;
     unsigned clsIdx = 0;
     drv1 = 0, drv2 = 0;
     
     MyNeiItType neiGrpIt(radius, grpPtr, grpPtr->GetLargestPossibleRegion() );
     neiGrpIt.OverrideBoundaryCondition(&constCondition);

     IteratorType3DVecUS histIt(histPtr, histPtr->GetLargestPossibleRegion() );
     float Q1 = 0;
     for (neiGrpIt.GoToBegin(), maskIt.GoToBegin(), histIt.GoToBegin(); !maskIt.IsAtEnd(); ++ neiGrpIt, ++maskIt, ++ histIt) {
	  if (maskIt.Get() > 0) {
	       curLabel = neiGrpIt.GetCenterPixel();
	       M0 = 0, M1 = 0, M2 = 0;
	       for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
		    a = (double)(histIt.Get()[clsIdx]) - par.numSubs;

		    // debug ///
		    itk::VariableLengthVector<unsigned short> histVector(par.numClusters);
		    histVector = histIt.Get();
		    /// end of debug ///
		    b = - (int(neiGrpIt.GetPixel(xminus) >= 0 && clsIdx != neiGrpIt.GetPixel(xminus) )
			   + int(neiGrpIt.GetPixel(xplus) >= 0 && clsIdx != neiGrpIt.GetPixel(xplus))
			   + int(neiGrpIt.GetPixel(yminus) >= 0 && clsIdx != neiGrpIt.GetPixel(yminus))
			   + int(neiGrpIt.GetPixel(yplus) >= 0 && clsIdx != neiGrpIt.GetPixel(yplus))
			   + int(neiGrpIt.GetPixel(zminus) >= 0 && clsIdx != neiGrpIt.GetPixel(zminus))
			   + int(neiGrpIt.GetPixel(zplus) >= 0 && clsIdx != neiGrpIt.GetPixel(zplus)) );

		    if (clsIdx == curLabel) {
			 acur = a;
			 bcur = b;
		    }
		    M0 += exp (a * par.alpha + b * par.betag);
		    M1 += a * exp(a * par.alpha + b * par.betag);
		    M2 += a * a * exp(a * par.alpha + b * par.betag);
	       } // for clsIdx
	       drv1 = drv1 - acur + M1 / M0;
	       drv2 += (M2 * M0 - M1 * M1)/ (M0 * M0);
	       Q1 = Q1 + (-acur * par.alpha - bcur * par.betag + log(M0));
	  }  // maskIt > 0
     } // for maskIt

     return Q1;
}

float SubAlphaDrv(ImageType3DChar::Pointer &  grpPtr,
	       std::vector<ImageType3DChar::Pointer>  & subVec,
	       ImageType3DChar::Pointer maskPtr,
	       double & drv1, double & drv2,
	       ParType & par)
{
     // mask image.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );


     // Define neighborhood iterator
     MyBoundCondType constCondition;
     constCondition.SetConstant(-1);
     MyNeiItType::RadiusType radius;
     radius.Fill(1);

     MyNeiItType::OffsetType xplus = {{1,0, 0}};
     MyNeiItType::OffsetType xminus = {{-1, 0, 0}};
     MyNeiItType::OffsetType yplus = {{0, 1, 0}};
     MyNeiItType::OffsetType yminus = {{0, -1, 0}};
     MyNeiItType::OffsetType zplus = {{0, 0, 1}};
     MyNeiItType::OffsetType zminus = {{0, 0, -1}};

     double a = 0, b = 0, M0 = 0, M1 = 0, M2 = 0, bcur = 0, acur = 0, Q2 = 0;
     unsigned curSubLabel = 0, grpLabel = 0;
     unsigned clsIdx = 0;
     drv1 = 0, drv2 = 0;
     
     IteratorType3DChar grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );

     for (unsigned short subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  MyNeiItType neiSampleIt( radius, subVec[subIdx], subVec[subIdx]->GetLargestPossibleRegion() );
	  neiSampleIt.OverrideBoundaryCondition(&constCondition);
	  
	  for (neiSampleIt.GoToBegin(), grpIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ neiSampleIt, ++ grpIt, ++ maskIt) {
	       if (maskIt.Get() > 0) {
		    curSubLabel = neiSampleIt.GetCenterPixel();

		    M0 = 0, M1 = 0, M2 = 0;
		    for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
			 a = - int (grpIt.Get() != clsIdx);
			 b = - (int(neiSampleIt.GetPixel(xminus) >= 0 && clsIdx != neiSampleIt.GetPixel(xminus) )
				+ int(neiSampleIt.GetPixel(xplus) >= 0 && clsIdx != neiSampleIt.GetPixel(xplus))
				+ int(neiSampleIt.GetPixel(yminus) >= 0 && clsIdx != neiSampleIt.GetPixel(yminus))
				+ int(neiSampleIt.GetPixel(yplus) >= 0 && clsIdx != neiSampleIt.GetPixel(yplus))
				+ int(neiSampleIt.GetPixel(zminus) >= 0 && clsIdx != neiSampleIt.GetPixel(zminus))
				+ int(neiSampleIt.GetPixel(zplus) >= 0 && clsIdx != neiSampleIt.GetPixel(zplus)) );
			 if (clsIdx == curSubLabel) {
			      acur = a;
			      bcur = b;
			 }
			 M0 += exp (a * par.alpha + b * par.betaz);
			 M1 += a * exp(a * par.alpha + b * par.betaz);
			 M2 += a * b * exp(a * par.alpha + b * par.betaz);
		    } // for clsIdx
		    drv1 += (- acur + M1 / M0);
		    drv2 += (M2 * M0 - M1 * M1)/ (M0 * M0);
		    Q2 = Q2 + (-acur * par.alpha - bcur * par.betaz + log (M0) );
	       }  // maskIt > 0
	  } // for maskIt
     } // subIdx

     return Q2;
}
