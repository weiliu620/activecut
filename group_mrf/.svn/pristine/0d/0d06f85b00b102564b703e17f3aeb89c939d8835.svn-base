#include "commonalt.h"
#include "MCModel.h"

std::vector<ImageType4DFloat::Pointer>  ReadGroupImage(std::string fmriPath);

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     float maskThresh = 0;

     std::string maskPath, maskimage, outMask, fmriPath;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "Update input mask image (-m)  and remove points in the mask where any fmri time series  (--fmripath) are zero. And apply threshold on mask to get a binary mask image (--outmask).")
	  ("mask,m", po::value<std::string>(&maskimage)->default_value("maskimage.nii"), "input maskimage file")
	  ("fmripath", po::value<std::string>(&fmriPath)->default_value("."),
	   "fmri image file path.")
	  ("outmask", po::value<std::string>(&outMask)->default_value("cleanmask.nii"),
	   "output clean mask image.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");


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


     // read fmri files.
     boost::filesystem::path fmriPathVar(fmriPath);
     boost::filesystem::directory_iterator fmriPathEnd;

     ImageType5DFloat::IndexType fmriIdx;          

     SeriesReaderType5DFloat::Pointer fmriReader = SeriesReaderType5DFloat::New();

     for (boost::filesystem::directory_iterator fmriPathIt(fmriPathVar); fmriPathIt != fmriPathEnd; ++fmriPathIt) {

     	  fmriReader->AddFileName( (*fmriPathIt).path().string());	  
     	  cout <<  "add " << (*fmriPathIt).path().string() << "\n";
     }

     fmriReader->Update();
     ImageType5DFloat::Pointer fmriPtr = fmriReader->GetOutput();
     ImageType5DFloat::SizeType fmriSize = fmriPtr->GetLargestPossibleRegion().GetSize();

     // read in mask image
     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(maskimage);
     maskReader->Update();
     ImageType3DChar::Pointer maskPtr = maskReader->GetOutput();

     ImageType3DChar::RegionType maskRegion;
     maskRegion = maskReader->GetOutput()->GetLargestPossibleRegion();
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();
     ImageType3DChar::IndexType maskIdx;     
     IteratorType3DChar maskIt(maskPtr, maskRegion);

     if (maskSize[0] != fmriSize[0] 
     	 || maskSize[1] != fmriSize[1] 
     	 || maskSize[2] != fmriSize[2]) {
     	  cout << "mask image and true label image have different size. Need double check before masking. Exit. " << endl;
	  exit(1);
     }
     
     vnl_vector <float> timeSeries(fmriSize[3], 0);

     unsigned numSupress = 0;
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() == 1) {
	       maskIdx = maskIt.GetIndex();
	       fmriIdx[0] = maskIdx[0];
	       fmriIdx[1] = maskIdx[1];
	       fmriIdx[2] = maskIdx[2];

	       fmriIdx[4] = 0;  // subIdx
	       do {
		    for (fmriIdx[3] = 0; fmriIdx[3] < fmriSize[3]; fmriIdx[3] ++) {
			 timeSeries[fmriIdx[3]] = fmriPtr->GetPixel(fmriIdx);
		    }
		    fmriIdx[4] = fmriIdx[4] + 1;
	       }
	       while (timeSeries.two_norm() == 0 && fmriIdx[4] < fmriSize[4] );
	  
	       if (timeSeries.two_norm() == 0) {
		    maskIt.Set(0);
		    numSupress ++;
	       }
	  }
     }

     printf("number of voxels supressed to zero in original mask file: %i.\n", numSupress);
     
     WriterType3DChar::Pointer writer = WriterType3DChar::New();
     
     writer->SetInput(maskReader->GetOutput() );
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

