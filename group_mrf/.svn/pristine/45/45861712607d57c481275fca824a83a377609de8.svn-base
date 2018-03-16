#include "commonalt.h"
#include "MCModel.h"

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     std::string maskPath, maskimage, outMask;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "Given all subjects' mask file (--maskpath), compute average mask file (--out).")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")
	  ("maskpath", po::value<std::string>(&maskPath)->default_value("."),
	   "noised image file")
	  ("out", po::value<std::string>(&outMask)->default_value("meanmask"),
	   "output clean mask image.");

     po::options_description cmdline_options;
     cmdline_options.add(generic);

     po::positional_options_description p;
     p.add("in", -1);

     po::variables_map vm;        
     po::store(po::command_line_parser(argc, argv).
	       options(cmdline_options).positional(p).run(), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: cleanmask [options]\n";
	       std::cout << cmdline_options << "\n";
	       return 0;
	  }
     }

     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    


     // read mask files.
     boost::filesystem::path maskPathVar(maskPath);
     boost::filesystem::directory_iterator maskPathEnd;

     ImageType4DFloat::IndexType maskIdx;          

     SeriesReaderType4DFloat::Pointer maskReader = SeriesReaderType4DFloat::New();

     for (boost::filesystem::directory_iterator maskPathIt(maskPathVar); maskPathIt != maskPathEnd; ++maskPathIt) {

     	  maskReader->AddFileName( (*maskPathIt).path().string());	  
     	  cout <<  "add " << (*maskPathIt).path().string() << "\n";
     }

     maskReader->Update();
     ImageType4DFloat::Pointer maskPtr = maskReader->GetOutput();
     ImageType4DFloat::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();;
     ImageType4DFloat::PointType maskOrigin = maskPtr->GetOrigin();
     ImageType4DFloat::SpacingType maskSpacing = maskPtr->GetSpacing();
     ImageType4DFloat::DirectionType maskDirection = maskPtr->GetDirection();
     

     // allocate memory for saving output combined mask.
     ImageType3DFloat::Pointer outPtr = ImageType3DFloat::New();
     ImageType3DFloat::IndexType outIdx;
     outIdx.Fill(0);

     ImageType3DFloat::SizeType outSize;
     outSize[0] = maskSize[0];
     outSize[1] = maskSize[1];
     outSize[2] = maskSize[2];
     ImageType3DFloat::RegionType outRegion;
     outRegion.SetSize(outSize);
     outRegion.SetIndex(outIdx);
     outPtr->SetRegions(outRegion);
     outPtr->Allocate();
     outPtr->FillBuffer(0);

     ImageType3DFloat::PointType outOrigin;
     ImageType3DFloat::SpacingType outSpacing;
     ImageType3DFloat::DirectionType outDirection;
     outOrigin[0] = maskOrigin[0];
     outOrigin[1] = maskOrigin[1];
     outOrigin[2] = maskOrigin[2];
     outSpacing[0] = maskSpacing[0];
     outSpacing[1] = maskSpacing[1];
     outSpacing[2] = maskSpacing[2];

     outDirection(0,0) = maskDirection(0,0);
     outDirection(0,1) = maskDirection(0,1);
     outDirection(0,2) = maskDirection(0,2);

     outDirection(1,0) = maskDirection(1,0);
     outDirection(1,1) = maskDirection(1,1);
     outDirection(1,2) = maskDirection(1,2);

     outDirection(2,0) = maskDirection(2,0);
     outDirection(2,1) = maskDirection(2,1);
     outDirection(2,2) = maskDirection(2,2);

     outPtr->SetOrigin(outOrigin);
     outPtr->SetSpacing(outSpacing);
     outPtr->SetDirection(outDirection);

     IteratorType3DFloat outIt(outPtr, outRegion);

     float meanValue = 0;
     for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++ outIt) {
	  outIdx = outIt.GetIndex();
	  maskIdx[0] = outIdx[0];
	  maskIdx[1] = outIdx[1];
	  maskIdx[2] = outIdx[2];

	  meanValue = 0;
	  for (maskIdx[3] = 0; maskIdx[3] < maskSize[3]; maskIdx[3] ++) {
	       meanValue += maskPtr->GetPixel(maskIdx);
	  }
	  meanValue = meanValue / maskSize[3];
	  outIt.Set(meanValue);
     }
     
     WriterType3DFloat::Pointer writer = WriterType3DFloat::New();
     
     writer->SetInput(outPtr);
     writer->SetFileName(outMask);
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

     std::cout << "File " << outMask<< " saved.\n";
     return 0;     

	       
	  
     

}

