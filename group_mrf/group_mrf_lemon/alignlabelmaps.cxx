#include <common.h>
#include <algorithm>
int main (int argc, char* argv[])
{
     unsigned verbose = 0;
     std::string inImage, refImage, outImage;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Given reference label map, permute target labelmap such that it has best consistency with reference label map.")
	  ("input,i", po::value<std::string>(&inImage)->default_value("inimage.nii.gz"), "Input image file that need to be aligned to the reference image.")
	  ("ref,r", po::value<std::string>(&refImage)->default_value("trueimage.nii"), "reference label map.")
	  ("out,o", po::value<std::string>(&outImage)->default_value("outimage.nii"), "aligned label map.")
	  ("verbose,v", po::value<unsigned int>(&verbose)->default_value(0), 
	   "Verbose level.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: alignlabelmaps [options]\n";
	       std::cout << mydesc << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    



     ReaderType3DShort::Pointer refReader = ReaderType3DShort::New();
     refReader->SetFileName(refImage);
     ImageType3DShort::Pointer refPtr = refReader->GetOutput();
     refPtr->Update();

     ReaderType3DShort::Pointer inputImageReader = ReaderType3DShort::New();
     inputImageReader->SetFileName(inImage);
     ImageType3DShort::Pointer inputImagePtr = inputImageReader->GetOutput();
     inputImagePtr->Update();

     // Compute the maximum intensity value of input image.
     typedef itk::MinimumMaximumImageCalculator <ImageType3DShort>  ImageCalculatorFilterType;
     ImageCalculatorFilterType::Pointer inputImageCalculatorFilter = ImageCalculatorFilterType::New ();
     inputImageCalculatorFilter->SetImage(inputImagePtr);
     inputImageCalculatorFilter->Compute();
     unsigned numClustersInputImage = inputImageCalculatorFilter->GetMaximum();

     ImageCalculatorFilterType::Pointer refCalculatorFilter = ImageCalculatorFilterType::New ();
     refCalculatorFilter->SetImage(refPtr);
     refCalculatorFilter->Compute();
     unsigned numClustersRef = refCalculatorFilter->GetMaximum();

     // init contingency table.
     unsigned U = numClustersInputImage;
     unsigned V = numClustersRef;
     vnl_matrix <double> conTable(U, V, 0);
     IteratorType3DShort inputIt(inputImagePtr, inputImagePtr->GetLargestPossibleRegion() );
     IteratorType3DShort refIt(refPtr, refPtr->GetLargestPossibleRegion() );

     // compute contingency table.
     unsigned N = 0; 
     for (inputIt.GoToBegin(), refIt.GoToBegin(); !inputIt.IsAtEnd(); ++ inputIt, ++ refIt) {
	  if (inputIt.Get() > 0 && refIt.Get() > 0) {
	       conTable[inputIt.Get()-1][refIt.Get()-1] ++;
	       N ++;
	  }
     }


     // Now permute the contingency table, which is equivalent to permute the
     // iput label map.
     std::vector<unsigned short> permumap(numClustersInputImage);
     std::vector<unsigned short> best_permumap(numClustersInputImage);

     // use 1-based cluster index.
     for (unsigned clsIdx = 0; clsIdx < numClustersInputImage; clsIdx ++) {
	  permumap[clsIdx] = clsIdx;
     }

     unsigned numMatches = 0, bestNumMatches = 0;
     unsigned minNumClusters = numClustersInputImage<numClustersRef?numClustersInputImage:numClustersRef;
     do {

	  if (verbose >= 1) {
	       printf("working on: ");
	       for (unsigned clsIdx = 0; clsIdx < numClustersInputImage; clsIdx ++) {
		    printf("%i->%i, ", clsIdx, permumap[clsIdx]);
	       }
	       printf("\n");
	  }

	  numMatches = 0;
	  for (unsigned diagIdx = 0; diagIdx < minNumClusters; diagIdx ++) {
	      numMatches += conTable(permumap[diagIdx], diagIdx);
	  }

	  if (numMatches > bestNumMatches) {
	       bestNumMatches = numMatches;
	       best_permumap = permumap;

	       if (verbose >= 1) {
		    printf(" best permumation map: ");
		    for (unsigned clsIdx = 0; clsIdx < numClustersInputImage; clsIdx ++) {
			 printf("%i->%i, ", clsIdx + 1, best_permumap[clsIdx] + 1);
		    }
		    printf("Percent of voxels that matches: %2.2f\n", (double)(bestNumMatches) / (double)(N));
	       } // verbose >= 1    
	  }
	  
     } while (next_permutation (permumap.begin(), permumap.begin() + numClustersInputImage) );

     printf("Finally, best permutation: ");
     for (unsigned clsIdx = 0; clsIdx < numClustersInputImage; clsIdx ++) {
	  printf("%i -> %i, ", clsIdx+1, best_permumap[clsIdx]+1);
     }
     
     printf("Percent of voxels that matches: %2.2f\n", (double)(bestNumMatches) / (double)(N));

     // Apply the permutation to ininput image.
     for (inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++ inputIt) {
	  if (inputIt.Get() > 0) { // no maskIt here.
	       inputIt.Set( best_permumap[inputIt.Get() - 1] + 1 );
	  }
     }
     
     WriterType3DShort::Pointer writer = WriterType3DShort::New();
     
     writer->SetInput(inputImagePtr);
     writer->SetFileName(outImage);
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

     std::cout << "alignlabelmaps: File " << outImage << " saved.\n";

	      

     return 0;
}
