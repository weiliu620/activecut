#include "common_gmmodel.h"
#include "gmmodel.h"

using std::cout;
using std::endl;

twister_base_gen_type mygenerator(42u);

unsigned ReadObs(std::string obspath,
		 VnlVectorImageType::Pointer & obsPtr,
		 ImageType3DChar::Pointer maskPtr);

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     unsigned int EMIter = 1;
     int numClusters = 4;
     std::string obspath, initGrplabel, outGrpLabel, outGrpProb, outSubBase;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Inference on group mrf.")
	  ("emiter", po::value<unsigned int>(&EMIter)->default_value(30),
	   "initial number of Monte Carlo samples. ")
	  ("numClusters,k", po::value<int>(&numClusters)->default_value(6),
	       "Number of clusters. Default is 6.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")

	  ("obspath,f", po::value<std::string>(&obspath)->default_value("."),
	   "noised image file")

	  ("grouplabel,g", po::value<std::string>(&outGrpLabel)->default_value("outgrplabel.nii"), "output group label file. nii format.")

	   ("initgrouplabel,i", po::value<std::string>(&initGrplabel)->default_value("grouplabel.nii"), "Initial group level label map. Also used as mask file.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: generateimage [options]\n";
	       cout << mydesc << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  cout << e.what() << "\n";
	  return 1;
     }    



     // read in initial grp label. Also used as mask.
     ReaderType3DChar::Pointer initGrpReader = ReaderType3DChar::New();
     initGrpReader->SetFileName(initGrplabel);
     ImageType3DChar::Pointer initGrpPtr = initGrpReader->GetOutput();
     initGrpPtr->Update();

     // original 1-based group map only used for mask.
     ImageType3DChar::Pointer maskPtr = initGrpReader->GetOutput();


     // Allocate memory for saving vectorized images.
     VnlVectorImageType::Pointer  obsPtr = VnlVectorImageType::New();
     ReadObs(obspath, obsPtr, maskPtr);
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();

     // Allocate memory for group label map. (vector image)
     ImageType3DVecF::Pointer grpPtr = ImageType3DVecF::New();
     ImageType3DVecF::IndexType grpIdx;
     grpIdx.Fill(0);

     ImageType3DVecF::SizeType grpSize;
     grpSize = obsSize;

     ImageType3DVecF::RegionType grpRegion;
     grpRegion.SetSize(grpSize);
     grpRegion.SetIndex(grpIdx);
     grpPtr->SetRegions(grpRegion);
     grpPtr->SetVectorLength(numClusters);
     grpPtr->SetNumberOfComponentsPerPixel(numClusters);
     grpPtr->Allocate();

     // copy initital group map to group map vector image.
     ReadInitGrpLabel(grpPtr, initGrpPtr);


     GMModel gmmodel(obsPtr,
		     grpPtr,
		     maskPtr,
		     numClusters,
     		     verbose);

     if (verbose >= 1) {
	  SaveGrpLabel(grpPtr, maskPtr, outGrpLabel);
     }

     gmmodel.EstimateMu(grpPtr, obsPtr, maskPtr);
     gmmodel.estimateSigma(grpPtr, obsPtr, maskPtr);
     gmmodel.printSelf("normal");
     
     for (unsigned short emIterIdx = 0; emIterIdx < EMIter; emIterIdx ++) {

	  if (verbose >= 1) {
	       printf("EM iteration %i begin:\n", emIterIdx + 1);
	  }

	  // sampling both group and individual subject label map.
	  gmmodel.Estep(grpPtr, obsPtr, maskPtr);

	  if (verbose >= 1) {
	       printf("EM iteration %i, E step done. \n", emIterIdx + 1);
	       SaveGrpLabel(grpPtr, maskPtr, outGrpLabel);
	       SaveGrpProb(grpPtr, maskPtr, "rawgroupmap.nii");
	  }

	  // estimate vMF parameters mu, kappa.
	  gmmodel.EstimateMu(grpPtr, obsPtr, maskPtr);
	  gmmodel.estimateSigma(grpPtr, obsPtr, maskPtr);

	  printf("EM iteration %i, parameter estimation done. \n", emIterIdx + 1 );
	  gmmodel.printSelf("normal");
	  fflush(stdout);
	  
	  gmmodel.llhood(grpPtr, obsPtr, maskPtr);

	  if (verbose >= 1) {
	       // Group label map.
	       SaveGrpLabel(grpPtr, maskPtr, outGrpLabel);
	  }
     }
     
     // Group label map.
     SaveGrpLabel(grpPtr, maskPtr, outGrpLabel);
}


