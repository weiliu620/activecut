#include "common_mcem.h"
#include <algorithm>
using std::cout;
using std::cin;
using std::string;

int main (int argc, char* argv[])
{
     unsigned verbose = 0;
     unsigned numClusters = 0;
     std::string trueimage, recoveredimage;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "given label map and mean time series, generate time series in von Mises-Fisher distribution.")
	  ("recovered,r", po::value<std::string>(&recoveredimage)->default_value("recoveredimage.nii.gz"), "recovered label map that need to compare with true label map.")
	  ("true,t", po::value<std::string>(&trueimage)->default_value("trueimage.nii"), "true label map.")
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
	       cout << mydesc << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  cout << e.what() << "\n";
	  return 1;
     }    

     // read in true label image
     ReaderType3DChar::Pointer trueImageReader = ReaderType3DChar::New();
     trueImageReader->SetFileName(trueimage);
     trueImageReader->Update();
     ImageType3DChar::Pointer trueImagePtr = trueImageReader->GetOutput();
     IteratorType3DChar trueImageIt(trueImagePtr, trueImageReader->GetOutput()->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType trueImageIdx;

     // read recovered label image.
     ReaderType3DChar::Pointer recoveredReader = ReaderType3DChar::New();
     recoveredReader->SetFileName(recoveredimage);
     recoveredReader->Update();
     ImageType3DChar::Pointer recoveredPtr = recoveredReader->GetOutput();
     IteratorType3DChar recoveredIt(recoveredPtr, recoveredReader->GetOutput()->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType recoveredIdx;


     std::vector<unsigned short> permumap(numClusters);
     for (unsigned clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  permumap[clsIdx] = clsIdx;
     }

     
     sort (permumap.begin(), permumap.begin() + numClusters);

     
     unsigned totalNumber = 0, correctedNumber = 0;
     unsigned predictedLabel = 0;
     float predRate = 0, thisPredRate = 0;
     do {

	  if (verbose >= 1) {
	       for (unsigned clsIdx = 1; clsIdx < numClusters + 1; clsIdx ++) {
		    std::cout << permumap[clsIdx] << " ";
	       }
	       std::cout << std::endl;
	  }

	  // begin compare images.
	  totalNumber = 0, correctedNumber = 0;
	  for (recoveredIt.GoToBegin(), trueImageIt.GoToBegin(); !recoveredIt.IsAtEnd(); ++ recoveredIt, ++ trueImageIt) {
	       if (recoveredIt.Get() > 0) {
		    totalNumber ++;
		    // 0-based --> 1-based --> 0-based.
		    predictedLabel = permumap[recoveredIt.Get() - 1] + 1;
		    if (predictedLabel == trueImageIt.Get() ) {
			 correctedNumber ++;
		    }
	       }
	  }
	  
	  // compare this permutation's prediction rate with other
	  // permutations'.
	  thisPredRate = float(correctedNumber) / float(totalNumber);
	  if (verbose >= 1) {
	       std::cout << "thisPredRate = " << thisPredRate << std::endl;
	  }

	  if (thisPredRate > predRate) {
	       predRate = thisPredRate;
	  }
     } while ( next_permutation (permumap.begin(), permumap.begin() + numClusters) );

     std::cout << "corrected clustering rate (in percentage): " << predRate << std::endl;
}

     

