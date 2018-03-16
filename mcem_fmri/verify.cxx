#include "commonalt.h"
#include "MCModel.h"

int main(int argc, char* argv[])
{
     int _DBG = 0;
     std::string trueimage;
     std::string labeled;
     
     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help", "produce help message")
	  ("truth,t", po::value<std::string>(&trueimage)->default_value("trueimage.nii"), "true labeled image file. ")
	  ("test,e", po::value<std::string>(&labeled)->default_value("labeled.nii"), "labeled image file by segmentation algorithm. ");

     po::options_description cmdline_options;
     cmdline_options.add(generic);

     po::positional_options_description p;
     p.add("in", -1);

     po::variables_map vm;        
     po::store(po::command_line_parser(argc, argv).
	       options(cmdline_options).positional(p).run(), vm);
     po::notify(vm);    

     try {

	  if (vm.count("help")) {
	       std::cout << "Usage: verifyimage [options]\n";
	       std::cout << cmdline_options << "\n";
	       return 0;
	  }
     
	  if (vm.count("version")) {
	       std::cout << " version 1.0.\n";
	       return 0;
	  }
     }

     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // read in image data.
     ReaderType::Pointer trueimageReader = ReaderType::New();
     ReaderType::Pointer labeledReader = ReaderType::New();
     trueimageReader->SetFileName(trueimage);
     labeledReader->SetFileName(labeled);

     trueimageReader->Update();
     labeledReader->Update();

     ImageType2D::Pointer trueimagePtr = trueimageReader->GetOutput();
     ImageType2D::Pointer labeledPtr = labeledReader->GetOutput();

     ImageType2D::RegionType trueimageRegion = trueimagePtr->GetLargestPossibleRegion();
     ImageType2D::RegionType labeledRegion = labeledPtr->GetLargestPossibleRegion();

     ImageType2D::SizeType size =
          trueimageReader->GetOutput()->GetLargestPossibleRegion().GetSize();

     ConstIteratorType2D trueimageIt(trueimagePtr, trueimageRegion);
     ConstIteratorType2D labeledIt(labeledPtr, labeledRegion);
     trueimageIt.GoToBegin();
     labeledIt.GoToBegin();
     unsigned numError = 0;
     while(!trueimageIt.IsAtEnd()) {
	  if (trueimageIt.Get() != labeledIt.Get()) {
	       numError++;
	  }
	  ++trueimageIt;
	  ++labeledIt;
     }

     printf("verify(): # of misclassified = %i, error rate = %.2f\% \n", numError, 100 * double(numError)/double(size[0]*size[1]));

     return 0;
}


