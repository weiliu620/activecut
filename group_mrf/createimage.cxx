#include <commonalt.h>
#include <utilalt.h>
using std::cout;
using std::string;

int applyMask(ImageType3DChar::Pointer imagePtr,
	      ImageType3DChar::Pointer maskPtr);

// Init random generator as a global object.
twister_base_gen_type mygenerator(42u);
int main (int argc, char* argv[])
{
     unsigned verbose = 0;
     unsigned seed = 0;
     float beta = 0.7;
     int scan = 0;
     int num_scan;
     int numClusters = 1;
     float denergy = 0; // energy change = candidate energy - current energy.
     int cand;
     double p_acpt = 0; // probability to accept the candidate sample.
     signed int currentLabel = 0;
     std::string observedimage, trueimage, maskimage;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "create group label map.")
	  ("beta,b", po::value<float>(&beta)->default_value(2), 
	   "Set MRF parameter beta.")
	  ("numclusters,k", po::value<int>(&numClusters)->default_value(4),
	       "Number of clusters.")
	  ("scan,s", po::value<int>(&num_scan)->default_value(500),
	   "Number of scan on MRF.")
	  ("true,t", po::value<std::string>(&trueimage)->default_value("trueimage.nii"), "output true image file")
	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("mask,m", po::value<std::string>(&maskimage)->default_value("maskimage.nii"), "input maskimage file")
	  ("verbose,v", po::value<unsigned int>(&verbose)->default_value(0), 
	   "Verbose level.");

     
     po::options_description cmdline_options;
     cmdline_options.add(generic);

     po::positional_options_description p;
     p.add("in", -1);

     po::variables_map vm;        
     po::store(po::command_line_parser(argc, argv).
	       options(cmdline_options).positional(p).run(), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") || argc == 1) {
	       std::cout << "Create group label map.\n";
	       cout << cmdline_options << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  cout << e.what() << "\n";
	  return 1;
     }    


     // read in mask image
     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(maskimage);
     maskReader->Update();
     ImageType3DChar::Pointer maskPtr = maskReader->GetOutput();

     ImageType3DChar::RegionType maskRegion;
     maskRegion = maskReader->GetOutput()->GetLargestPossibleRegion();
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();
     IteratorType3DChar maskIt(maskPtr, maskRegion);
     ImageType3DChar::IndexType maskIdx;     



     // Allocate memory for images.
     ImageType3DChar::Pointer labelImagePtr = ImageType3DChar::New();
     ImageType3DChar::IndexType start;
     start[0] = 0;
     start[1] = 0;
     start[2] = 0;

     ImageType3DChar::SizeType imageSize; 
     // imageSize[0] = imageSize0;
     // imageSize[1] = imageSize1;
     // imageSize[2] = imageSize2;

     imageSize = maskSize;

     ImageType3DChar::RegionType region;
     region.SetSize(imageSize);
     region.SetIndex(start);
     labelImagePtr->SetRegions(region);
     labelImagePtr->Allocate();
     
     // mygenerator.seed(static_cast<unsigned int>(seed));
     mygenerator.seed(static_cast<unsigned int>(std::time(0)));

     // Uniform integer generator.
     boost::uniform_int<> uni_int(0, numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // normal distribution generator.
     boost::normal_distribution<> normal_dist(0, 1); // Normal distribution.
     boost::variate_generator<twister_base_gen_type&, boost::normal_distribution<> > nor(mygenerator, normal_dist);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     // Define neighborhood iterator
     typedef itk::ConstantBoundaryCondition< ImageType3DChar >  BoundaryConditionType;
     BoundaryConditionType constCondition;
     constCondition.SetConstant(-1);
     typedef itk::NeighborhoodIterator< ImageType3DChar, BoundaryConditionType > NeighborhoodIteratorType;
     
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);

     NeighborhoodIteratorType neighborIt( radius, labelImagePtr, labelImagePtr->GetRequestedRegion() );
     neighborIt.OverrideBoundaryCondition(&constCondition);

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1}};

     typedef itk::ImageRegionIterator< ImageType3DChar> IteratorType;
     IteratorType labelImageIt(labelImagePtr, labelImagePtr->GetRequestedRegion());

     // Init image.
     for (labelImageIt.GoToBegin(); !labelImageIt.IsAtEnd(); ++labelImageIt) {
     	  labelImageIt.Set(roll_die());
     }

     // ******* mask image so voxels outside of mask will be negative. ****** //
     applyMask(labelImagePtr, maskPtr);

     ImageType3DChar::IndexType imageIdx;     
     // sampling Markov Random Field.
     for (scan = 0; scan < num_scan; scan++) {
     	  for (neighborIt.GoToBegin(); !neighborIt.IsAtEnd(); ++neighborIt) {
	       currentLabel = neighborIt.GetCenterPixel();
	       if (currentLabel >= 0) {
		    cand = roll_die();
		    imageIdx = neighborIt.GetIndex();

		    denergy 
			 = int(cand != neighborIt.GetPixel(xminus)) - int(currentLabel != neighborIt.GetPixel(xminus))
			 + int(cand != neighborIt.GetPixel(xplus)) - int(currentLabel != neighborIt.GetPixel(xplus))
			 + int(cand != neighborIt.GetPixel(yminus)) - int(currentLabel != neighborIt.GetPixel(yminus))
			 + int(cand != neighborIt.GetPixel(yplus)) - int(currentLabel != neighborIt.GetPixel(yplus))
			 + int(cand != neighborIt.GetPixel(zminus)) - int(currentLabel != neighborIt.GetPixel(zminus))
			 + int(cand != neighborIt.GetPixel(zplus)) - int(currentLabel != neighborIt.GetPixel(zplus));

		    denergy = beta * denergy;
		    // if energy change less than zero, just accept
		    // candidate. otherwise accept with exp(- energy
		    // change).
		    if (denergy <= 0) {
			 neighborIt.SetCenterPixel(cand);
		    }
		    else {
			 p_acpt = exp(-denergy);
			 if (uni() < p_acpt) {
			      neighborIt.SetCenterPixel(cand);
			 }
		    }
	       } // in mask
	  } // iterator.

     	  if (scan%10 == 0 ) { 
	       printf("scan %d. \n", scan); 
	  }
     } // end for (scan)
     save3dcharInc(labelImagePtr, trueimage);

}

int applyMask(ImageType3DChar::Pointer imagePtr,
	      ImageType3DChar::Pointer maskPtr)
{
     
     // images
     ImageType3DChar::IndexType imageIdx;
     ImageType3DChar::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();

     // mask.
     ImageType3DChar::RegionType maskRegion;
     maskRegion = maskPtr->GetLargestPossibleRegion();
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();
     IteratorType3DCharIdx maskIt(maskPtr, maskRegion);
     ImageType3DChar::IndexType maskIdx;    


     if (maskSize[0] != imageSize[0] 
	 || maskSize[1] != imageSize[1] 
	 || maskSize[2] != imageSize[2]) {
	  std::cout << "mask image and label image have different size. Need double check before masking. Exit. " << std::endl;
     }

     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  maskIdx = maskIt.GetIndex();
	  imageIdx[0] = maskIdx[0];
	  imageIdx[1] = maskIdx[1];
	  imageIdx[2] = maskIdx[2];

	  // For pixels outside of the mask, we set it to -1. So later
	  // routine will know this is outside of mask. No mask needed later!
	  if (maskIt.Get() <= 0) {
	       imagePtr->SetPixel(imageIdx, -1);
	  }
     }
     return (0);
}
