#include "commonalt.h"
using std::cout;
using std::endl;
using std::string;

int main (int argc, char* argv[])
{
     unsigned verbose = 0;
     short outsideValue = 0;
     std::string trueimage, maskimage, outimage;
     float threshhold = 0;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "produce help message")
	  ("true,t", po::value<std::string>(&trueimage)->default_value("trueimage.nii"), "input true image file")
	  ("mask,m", po::value<std::string>(&maskimage)->default_value("maskimage.nii"), "input maskimage file")
	  ("out,o", po::value<std::string>(&outimage)->default_value("trueimage_m.nii"), "output masked image")
	  ("verbose,v", po::value<unsigned int>(&verbose)->default_value(0), 
	   "verbose level. 0 is no verbose.")
	  ("threshhold,t", po::value<float>(&threshhold)->default_value(0), 
	   "Threshhold for mask image. Only pixels greater than threshhold will be regared as image, otherwise will be outside.")
	  ("outsidevalue,u", po::value<short>(&outsideValue)->default_value(-1), 
	   "pixel value outside of the mask.");

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
	       std::cout << "Usage: generateimage [options]\n";
	       cout << cmdline_options << "\n";
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

     ImageType3DChar::RegionType trueImageRegion;
     trueImageRegion = trueImageReader->GetOutput()->GetLargestPossibleRegion();
     ImageType3DChar::SizeType trueImageSize = trueImageRegion.GetSize();

     IteratorType3DChar trueImageIt(trueImagePtr, trueImageRegion);
     ImageType3DChar::IndexType trueImageIdx;     

     // read in mask image
     ReaderType3DFloat::Pointer maskReader = ReaderType3DFloat::New();
     maskReader->SetFileName(maskimage);
     maskReader->Update();
     ImageType3DFloat::Pointer maskPtr = maskReader->GetOutput();

     ImageType3DFloat::RegionType maskRegion;
     maskRegion = maskReader->GetOutput()->GetLargestPossibleRegion();
     ImageType3DFloat::SizeType maskSize = maskRegion.GetSize();
     IteratorType3DFloat maskIt(maskPtr, maskRegion);
     ImageType3DFloat::IndexType maskIdx;     


     if (maskSize[0] != trueImageSize[0] 
	 || maskSize[1] != trueImageSize[1] 
	 || maskSize[2] != trueImageSize[2]) {
	  cout << "mask image and true label image have different size. Need double check before masking. Exit. " << endl;
     }
     
     for (maskIt.GoToBegin(), trueImageIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt, ++trueImageIt)
     {
	  trueImageIdx = trueImageIt.GetIndex();
	  maskIdx = maskIt.GetIndex();
	  if (maskIt.Get() <= threshhold) {
	       trueImageIt.Set(outsideValue);
	  }
     }

     WriterType3DChar::Pointer writer = WriterType3DChar::New();
     writer->SetInput(trueImagePtr);
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

     return 0;
}
