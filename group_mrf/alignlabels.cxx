#include "commonalt.h"
#include <algorithm>
namespace po = boost::program_options;
int main (int argc, char* argv[])
{
     unsigned verbose = 0;
     unsigned numClusters = 0;
     std::string refimage, targetimage, outimage;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Given reference label map, permute target labelmap such that it has best consistency with reference label map.")
	  ("target,t", po::value<std::string>(&targetimage)->default_value("targetimage.nii.gz"), "recovered label map that need to compare with reference label map.")
	  ("ref,r", po::value<std::string>(&refimage)->default_value("refimage.nii"), "reference label map.")
	  ("out,o", po::value<std::string>(&outimage)->default_value("outimage.nii"), "aligned label map.")
	  ("numClusters,k", po::value<unsigned>(&numClusters)->default_value(6),
	       "Number of clusters. Default is 6.")
	  ("verbose,v", po::value<unsigned int>(&verbose)->default_value(0), 
	   "Verbose level.");

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

     // read in true label image
     ReaderType3DChar::Pointer refimageReader = ReaderType3DChar::New();
     refimageReader->SetFileName(refimage);
     refimageReader->Update();
     ImageType3DChar::Pointer refimagePtr = refimageReader->GetOutput();
     IteratorType3DChar refimageIt(refimagePtr, refimageReader->GetOutput()->GetLargestPossibleRegion() );

     // read recovered label image.
     ReaderType3DChar::Pointer recoveredReader = ReaderType3DChar::New();
     recoveredReader->SetFileName(targetimage);
     recoveredReader->Update();
     ImageType3DChar::Pointer recoveredPtr = recoveredReader->GetOutput();
     IteratorType3DChar recoveredIt(recoveredPtr, recoveredReader->GetOutput()->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType recoveredIdx;


     std::vector<unsigned short> permumap(numClusters);
     std::vector<unsigned short> best_permumap(numClusters);

     // use 1-based cluster index.
     for (unsigned clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  permumap[clsIdx] = clsIdx;
     }
     
     sort (permumap.begin(), permumap.begin() + numClusters);
     
     unsigned totalNumber = 0, correctedNumber = 0;
     unsigned predictedLabel = 0;
     float predRate = 0, thisPredRate = 0;
     do {
	  // begin compare images.
	  totalNumber = 0, correctedNumber = 0;
	  for (recoveredIt.GoToBegin(), refimageIt.GoToBegin(); !recoveredIt.IsAtEnd(); ++ recoveredIt, ++ refimageIt) {
	       if (recoveredIt.Get() > 0) { // no maskIt here.
		    totalNumber ++;
		    // 0-based --> 1-based --> 0-based.
		    predictedLabel = permumap[recoveredIt.Get() - 1] + 1;
		    if (predictedLabel == refimageIt.Get() ) {
			 correctedNumber ++;
		    }
	       }
	  }
	  
	  // compare this permutation's prediction rate with other
	  // permutations'.
	  thisPredRate = float(correctedNumber) / float(totalNumber);

	  if (verbose >= 1) {
	       printf("current permutation: ");
	       for (unsigned clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
		    printf("%i -> %i  ", clsIdx+1, permumap[clsIdx]+1);
	       }
	       std::cout << "thisPredRate = " << thisPredRate << std::endl;
	  }

	  if (thisPredRate > predRate) {
	       predRate = thisPredRate;
	       best_permumap = permumap;
	  }
     } while ( next_permutation (permumap.begin(), permumap.begin() + numClusters) );


     printf("best permutation: ");
     for (unsigned clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  printf("%i->%i, ", clsIdx+1, best_permumap[clsIdx]+1);
     }
     printf("\n");
     printf("Matching rate: %.2f\n",predRate);
     
     // apply the permutation so the label map get aligned.
     for (recoveredIt.GoToBegin(); !recoveredIt.IsAtEnd(); ++ recoveredIt) {
	  if (recoveredIt.Get() > 0) { // no maskIt here.
	       recoveredIt.Set( best_permumap[recoveredIt.Get() - 1] + 1 );
	  }
     }
     
     WriterType3DChar::Pointer writer = WriterType3DChar::New();
     
     writer->SetInput(recoveredPtr);
     writer->SetFileName(outimage);
     try 
     { 
	  writer->Update(); 
     } 
     catch( itk::ExceptionObject & err ) 
     { 
	  std::cerr << "ExceptionObject caught !" << std::endl; 
	  std::cerr << err << std::endl; 
	  return EXIT_FAILURE;
     } 

     std::cout << "alignlabels(): File " << outimage << " saved.\n";

     return 0;
}

     

