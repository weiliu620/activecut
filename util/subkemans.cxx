#include "common_gmmodel.h"

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     unsigned int kmeansIter = 1;
     int numClusters = 4;
     std::string obs, outLabel, maskimage;
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
	  ("obs,f", po::value<std::string>(&obsimage)->default_value("image.nii.gz"),
	   "noised image file")
	  ("outlabel,o", po::value<std::string>(&outLabel)->default_value("outlabel.nii"), "output label file. nii.gz format.");

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

     // read in observed image.
     ReaderType3DChar::Pointer obsReader = ReaderType3DChar::New();
     obsReader->SetFileName(obsimage);
     obsReader->Update();
     ImageType3DChar::Pointer obsPtr = obsReader->GetOutput();


     // create image buffer for out label.
     ImageType3DChar::IndexType start;
     start.Fill(0);
     ImageType3DChar::RegionType outRegion;
     outRegion.SetSize(maskSize);
     outRegion.SetIndex(start);
     ImageType3DChar::Pointer outPtr = ImageType3DChar::New();
     outPtr->SetRegions(outRegion);
     outPtr->Allocate();
     outPtr->FillBuffer(-1);

     // define a cluster center and bestcc.
     std::vector<float> cc(numClusters);
     std::vector<float> bestcc(numClusters);

     unsigned repeatIdx = 0; // repeat index of K-means.
     double minMeanSSR = 1e10, meanSSR = 0;

     for (repeatIdx = 0; repeatIdx < kmeansIter; repeatIdx ++) {
     	  printf("kmeans run  %i begin:\n", repeatIdx+1);
     	  meanSSR = kmeans(obsPtr, outPtr, maskPtr, cc, 1e-6);

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
     


double kmeans(ImageType3DChar::Pointer obsPtr,
	      ImageType3DChar::Pointer groupPtr,
	      ImageType3DChar::Pointer maskPtr,
	      std::vector< float > & cc,
	      float converge_eps)

{
     unsigned clsIdx = 0;
     unsigned numClusters = cc.size();
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();

     // Init kmeans by kmeans++ method.
     kmeanspp(obsPtr, maskPtr, cc);
     printf("after kmans++, cluster centers:\n");
     PrintClusterCenter(cc);

     // Define a prevcc to save previous cc.
     std::vector< float > prevcc(numClusters);

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

     }
     while(!converge);
}
