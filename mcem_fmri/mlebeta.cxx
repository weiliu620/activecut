#include "commonalt.h"
#include "MCModel.h"

int main(int argc, char* argv[])
{
     unsigned verbose = 0;
     float beta = 0.1;
     unsigned numScan = 1;

     int nlabels = 5;

     std::string observedimage, trueimage;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "produce help message")
	  ("beta,b", po::value<float>(&beta)->default_value(0.1), 
	   "MRF smoothness prior. When cheat = no, used as initnial value of EM. When cheat = yes, used as true value.")
	  ("numclusters,l", po::value<int>(&nlabels)->default_value(5),
	       "Number of labels. Default is 5.")
	  ("obs", po::value<std::string>(&observedimage)->default_value("observedimage.nii"),
	   "noised image file")
	  ("true", po::value<std::string>(&trueimage)->default_value("trueimage.nii"),
	   "true image file");

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
	       std::cout << cmdline_options << "\n";
	       return 0;
	  }
	  
     }

     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // read in image data.
     ReaderType::Pointer reader = ReaderType::New();
     reader->SetFileName(observedimage);
     reader->Update();
     ImageType2D::Pointer imagePtr = reader->GetOutput();
     ImageType2D::SizeType size =
          reader->GetOutput()->GetLargestPossibleRegion().GetSize();

     // read in true label data.
     ReaderType::Pointer trueReader = ReaderType::New();
     trueReader->SetFileName(trueimage);
     trueReader->Update();
     ImageType2D::Pointer truePtr = trueReader->GetOutput();
     ImageType2D::SizeType trueSize =
          trueReader->GetOutput()->GetLargestPossibleRegion().GetSize();

     // Allocate memory for saving samples of MC.
     ImageType3D::Pointer samplePtr = ImageType3D::New();
     ImageType3D::IndexType start;
     start[0] = 0;
     start[1] = 0;
     start[2] = 0;
     ImageType3D::SizeType sampleSize;
     sampleSize[0] = size[0];
     sampleSize[1] = size[1];


     sampleSize[2] = 1;
     numScan = 1;

     ImageType3D::RegionType region;
     region.SetSize(sampleSize);
     region.SetIndex(start);
     samplePtr->SetRegions(region);
     samplePtr->Allocate();

     MCModel mcmodel(imagePtr,
		     samplePtr,
		     nlabels,
		     numScan,
		     beta,
		     "icm",
		     verbose);

     // assign true labels to last slice of samplePtr.
     ImageType3D::IndexType sampleIdx;     
     ImageType2D::IndexType trueIdx;     
     sampleIdx[2] = sampleSize[2] - 1;
     for (sampleIdx[0] = 0; sampleIdx[0] < sampleSize[0]; sampleIdx[0] ++) {
	  for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1] ++) {
	       trueIdx[0] = sampleIdx[0];
	       trueIdx[1] = sampleIdx[1];
	       samplePtr->SetPixel(sampleIdx, truePtr->GetPixel(trueIdx));
	  }
     }

     printf("samplePtr assigned to true label. \n");
     mcmodel.estimatePriorParameter(0.1, 8, samplePtr);
     mcmodel.print("normal");
     mcmodel.estPriorPramNewton(samplePtr, 5);
     mcmodel.print("normal");

}
