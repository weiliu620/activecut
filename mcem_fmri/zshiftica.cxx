
#include "commonalt.h"

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;

     int numClusters = 5;
     unsigned short targetLabel = 0;
     float maxDotProd = 0.3;
     bool doAllLabels = 0;
     std::string fmriFile, inputFile, maskFile, outFile;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", " Normalize input 3D volume image so the voxel intensities have zero mean and standard deviation. The (smaple) mean and variance are computed from the image data.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")
	  ("input,i", po::value<std::string>(&inputFile)->default_value("observedimage.nii"),
	   "input 3D volume image file.")
	  ("mask,m", po::value<std::string>(&maskFile)->default_value("mask.nii"),
	   "mask file. Only voxels inside the mask are counted. Mask file need to be binary. mask > 0 is assumed to be in the mask.")
	  ("out,o", po::value<std::string>(&outFile)->default_value("normalized.nii"),
	   "output normalized image file");

     po::options_description cmdline_options;
     cmdline_options.add(generic);

     po::positional_options_description p;
     p.add("in", -1);

     po::variables_map vm;        
     po::store(po::command_line_parser(argc, argv).
	       options(cmdline_options).positional(p).run(), vm);
     po::notify(vm);    

     try {

	  if (vm.count("help") | argc == 1) {
	       std::cout << "Usage: " << argv[0] << " [options]\n";
	       std::cout << cmdline_options << "\n";
	       return 0;
	  }
	  
	  if (vm.count("label")) {
	       doAllLabels = 0;
	  }
	  else {
	       doAllLabels = 1;
	  }
     }

     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // read in file.
     ReaderType3DFloat::Pointer inputReader = ReaderType3DFloat::New();
     inputReader->SetFileName(inputFile);
     inputReader->Update();
     ImageType3DFloat::Pointer inputPtr = inputReader->GetOutput();
     ImageType3DFloat::SizeType inputSize = inputPtr->GetLargestPossibleRegion().GetSize();
     ImageType3DFloat::IndexType inputIdx;
     IteratorType3DFloat inputIt(inputPtr, inputPtr->GetLargestPossibleRegion() );

     // read mask file.
     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(maskFile);
     ImageType3DChar::Pointer maskPtr = maskReader->GetOutput();
     maskPtr->Update();
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType maskIdx;


     double mean = 0;
     double numPoints = 0;
     double std = 0;
     for (maskIt.GoToBegin(), inputIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++ inputIt) {
	  maskIdx = maskIt.GetIndex();
	  // inputIdx[0] = maskIdx[0];
	  // inputIdx[1] = maskIdx[1];
	  // inputIdx[2] = maskIdx[2];
	  
	  if (maskIt.Get() > 0) {
	       numPoints ++;
	       mean += inputIt.Get();
	  }
     }
     
     // compute mean.
     mean = mean / numPoints;
     

     for (maskIt.GoToBegin(), inputIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++ inputIt) {
	  if (maskIt.Get() > 0) {
	       std = std + (inputIt.Get() - mean) * (inputIt.Get() - mean);
	  }
     }
     std = sqrt ( std / (numPoints - 1));

     for (maskIt.GoToBegin(), inputIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++ inputIt) {
	  if (maskIt.Get() > 0) {
	       inputIt.Set( (inputIt.Get() - mean ) / std );
	  }
     }
     

     // write back.
     WriterType3DFloat::Pointer writer = WriterType3DFloat::New();
     writer->SetInput(inputPtr);
     writer->SetFileName(outFile);
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

     std::cout << "File " << outFile << " saved.\n";



}
