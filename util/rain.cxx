#include <common.h>
using namespace std;

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     std::string image1, image2, mask;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "output histograme of image after being masked.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")
	  ("ref,r", po::value<std::string>(&image1)->default_value("reference.nii.gz"), 
	   "reference image.")
	  ("image,i", po::value<std::string>(&image2)->default_value("image.nii.gz"), 
	   "test image.")
	  ("mask,m", po::value<std::string>(&mask)->default_value("mask.nii.gz"), 
	   "mask binary image.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: rain [options]\n";
	       cout << mydesc << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  cout << e.what() << "\n";
	  return 1;
     }    


     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(mask);
     ImageType3DChar::Pointer maskPtr = maskReader->GetOutput();
     maskPtr->Update();

     ReaderType3DChar::Pointer image1Reader = ReaderType3DChar::New();
     image1Reader->SetFileName(image1);
     ImageType3DChar::Pointer image1Ptr = image1Reader->GetOutput();
     image1Ptr->Update();

     ReaderType3DChar::Pointer image2Reader = ReaderType3DChar::New();
     image2Reader->SetFileName(image2);
     ImageType3DChar::Pointer image2Ptr = image2Reader->GetOutput();
     image2Ptr->Update();

     IteratorType3DChar image1ItA(image1Ptr, image1Ptr->GetLargestPossibleRegion());
     IteratorType3DChar image1ItB(image1Ptr, image1Ptr->GetLargestPossibleRegion());

     IteratorType3DChar image2ItA(image2Ptr, image2Ptr->GetLargestPossibleRegion());
     IteratorType3DChar image2ItB(image2Ptr, image2Ptr->GetLargestPossibleRegion());

     IteratorType3DChar maskItA(maskPtr, maskPtr->GetLargestPossibleRegion());
     IteratorType3DChar maskItB(maskPtr, maskPtr->GetLargestPossibleRegion());


     unsigned totalPts = 0;
     unsigned truePos = 0;
     unsigned trueNeg = 0;
     for (maskItA.GoToBegin(); !maskItA.IsAtEnd(); ++ maskItA) {
     	  if (maskItA.Get() > 0) {
     	       totalPts ++;
     	  }  // maskIt > 0
     } // maskIt

     ImageType3DChar::IndexType posIndex;
     for (maskItA.GoToBegin(), image1ItA.GoToBegin(), image2ItA.GoToBegin(); 
	  !maskItA.IsAtEnd(); 
	  ++ maskItA, ++image1ItA, ++ image2ItA) {
	  if (maskItA.Get() > 0) {
	       posIndex = maskItA.GetIndex();
	       maskItB.SetIndex(posIndex), ++ maskItB;
	       image1ItB.SetIndex(posIndex), ++ image1ItB;
	       image2ItB.SetIndex(posIndex), ++ image2ItB;
	       for (; !maskItB.IsAtEnd(); ++ maskItB, ++ image1ItB, ++ image2ItB) {
	  	    if (maskItB.Get() > 0) {
	  		 if ((image1ItA.Get() == image1ItB.Get()) && (image2ItA.Get() == image2ItB.Get()) ) {
	  		      // true positive. This voxel pair in same set for
	  		      // both iamges.
	  		      truePos ++;
			      
	  		 } // if
	  		 else if ((image1ItA.Get() != image1ItB.Get()) && (image2ItA.Get() != image2ItB.Get()) ) {
	  		      // true negative. This pair is in different set
	  		      // for both iamges.
	  		      trueNeg ++;

	  		 } // else if
	  		 else {
	  		 }
	  	    } // maskItB > 0
	       } // maskItB 

	  } // maskItA > 0
     } // maskItA

     unsigned long allComb = 0;
     allComb = boost::math::binomial_coefficient <double> (totalPts, 2);

     // printf("allComb = %i, truePos = %i, trueNeg = %i\n", allComb, truePos, trueNeg);
     
     // printf("Rand index (unjusted): %f. Percentage: %f. \n", float(truePos + trueNeg) / allComb,  100 * float(truePos + trueNeg) / allComb);

     printf("%s: Rand index (unjusted) in Percentage: %6.3f. \n", image2.c_str(), 100 * float(truePos + trueNeg) / allComb);

}


