#include "common_mcem.h"
using std::cout;
using std::string;

int main (int argc, char* argv[])
{

     
     int imageSize0 = 0;
     int imageSize1 = 0;
     int imageSize2 = 0;
     
     unsigned short radius1 = 0;
     unsigned short radius2 = 0;
     

     unsigned short verbose = 0;
     std::string maskImage;  

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "create binary mask image. a sphere.")

	  ("imagesizex,X", po::value<int>(&imageSize0)->default_value(64),
	   "Image size along x axis.")
	  ("imagesizey,Y", po::value<int>(&imageSize1)->default_value(64),
	   "Image size along y axis")
	  ("imagesizez,Z", po::value<int>(&imageSize2)->default_value(64),
	   "Image size along z axis.")
	  ("radius1,R", po::value<unsigned short>(&radius1)->default_value(64),
	   "mask sphere's radius.")
	  ("radius2,r", po::value<unsigned short>(&radius2)->default_value(0),
	   "mask sphere's inner radius. voxels in this radius will be masked out.")
	  ("maskimage,m", po::value<std::string>(&maskImage)->default_value("maskimage.nii"), "output mask image file")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
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


     // Allocate memory for saving mask image.
     ImageType3DChar::Pointer maskImagePtr = ImageType3DChar::New();
     ImageType3DChar::IndexType start;
     start[0] = 0;
     start[1] = 0;
     start[2] = 0;

     ImageType3DChar::SizeType imageSize; 
     imageSize[0] = imageSize0;
     imageSize[1] = imageSize1;
     imageSize[2] = imageSize2;

     ImageType3DChar::RegionType region;
     region.SetSize(imageSize);
     region.SetIndex(start);
     maskImagePtr->SetRegions(region);
     maskImagePtr->Allocate();

     ImageType3DChar::IndexType center;     
     ImageType3DChar::IndexType maskIdx;     
     center[0] = imageSize[0]/2;
     center[1] = imageSize[1]/2;
     center[2] = imageSize[2]/2;

     printf("center[2] = %i\n", center[2]);
     printf("radius2 = %i\n", radius2);
     
     IteratorType3DChar maskIt(maskImagePtr, region);
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt) {
	  maskIdx = maskIt.GetIndex();
	  if ((maskIdx[0] - center[0]) * (maskIdx[0] - center[0])
	      + (maskIdx[1] - center[1]) * (maskIdx[1] - center[1])
	      + (maskIdx[2] - center[2]) * (maskIdx[2] - center[2])
	      < radius1 * radius1
	       && (maskIdx[0] - center[0]) * (maskIdx[0] - center[0])
	      + (maskIdx[1] - center[1]) * (maskIdx[1] - center[1])
	      + (maskIdx[2] - center[2]) * (maskIdx[2] - center[2])
	      >= radius2 * radius2){
	       // inside radius1 and outside of raius2.
	       maskIt.Set(1);
	  }
	  else {
	       maskIt.Set(0);
	  }
     }

     WriterType3DChar::Pointer writer = WriterType3DChar::New();
     writer->SetInput(maskImagePtr);
     writer->SetFileName(maskImage);

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

     std::cout << "save3dchar(): File " << maskImage << " saved.\n";
     return 0;
}
