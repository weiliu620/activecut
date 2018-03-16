#include "commonalt.h"

// #include "vnl/algo/vnl_qr.h"
// #include "vnl/algo/vnl_matrix_inverse.h"
// #include "vnl/vnl_inverse.h"

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;

     std::string infile, maskfile, outfile;
     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "project time series onto a sphere.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")

	  ("in,i", po::value<std::string>(&infile)->default_value("observedimage.nii.gz"),
	   "original input nii file.")
	  ("mask,m", po::value<std::string>(&maskfile)->default_value("maskimage.nii.gz"),
	   "original input nii file.")
	  ("out,o", po::value<std::string>(&outfile)->default_value("fmrionsphere.nii"),
	   "output file. Time series vector are subtracted mean and normalized to unit length.");


     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: projectfmri [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // read in image data.
     ReaderType4DFloat::Pointer imageReader = ReaderType4DFloat::New();
     imageReader->SetFileName(infile);
     imageReader->Update();
     ImageType4DFloat::Pointer imagePtr = imageReader->GetOutput();
     ImageType4DFloat::SizeType imageSize = imagePtr->GetLargestPossibleRegion().GetSize();
     unsigned tsLength = imageSize[3];
     ImageType4DFloat::IndexType imageIdx;    

     // Allocate memory for saving out image, and std image.
     ImageType4DFloat::Pointer outPtr = ImageType4DFloat::New();
     ImageType4DFloat::IndexType outStart;
     ImageType4DFloat::IndexType outIdx;
     outStart.Fill(0);

     ImageType4DFloat::RegionType outRegion;
     outRegion.SetSize(imageSize);
     outRegion.SetIndex(outStart);
     outPtr->SetRegions(outRegion);
     outPtr->Allocate();
     outPtr->FillBuffer(0);
     IteratorType4DFloat outIt(outPtr, outRegion);

     // Allocate memory for saving middle image, and std image.
     ImageType3DFloat::Pointer middlePtr = ImageType3DFloat::New();
     ImageType3DFloat::IndexType middleStart;
     ImageType3DFloat::IndexType middleIdx;
     middleStart.Fill(0);

     ImageType3DFloat::SizeType middleSize;
     middleSize[0] = imageSize[0];
     middleSize[1] = imageSize[1];
     middleSize[2] = imageSize[2];

     ImageType3DFloat::RegionType middleRegion;
     middleRegion.SetSize(middleSize);
     middleRegion.SetIndex(middleStart);
     middlePtr->SetRegions(middleRegion);
     middlePtr->Allocate();
     middlePtr->FillBuffer(0);
     IteratorType3DFloat middleIt(middlePtr, middleRegion);

     double thisMean = 0;
     for (middleIt.GoToBegin(); !middleIt.IsAtEnd(); ++ middleIt)
     {

	  middleIdx = middleIt.GetIndex();
	  imageIdx[0] = middleIdx[0];
	  imageIdx[1] = middleIdx[1];
	  imageIdx[2] = middleIdx[2];

	  thisMean = 0;
	  for (unsigned tsIdx = 0; tsIdx < tsLength; tsIdx ++) {
	       imageIdx[3] = tsIdx;
	       thisMean += imagePtr->GetPixel(imageIdx);
	  }
	  thisMean = thisMean / tsLength;
	  middleIt.Set(thisMean);
     }

     // substract mean.
     for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++ outIt)
     {
	  outIdx = outIt.GetIndex();
	  middleIdx[0] = outIdx[0];
	  middleIdx[1] = outIdx[1];
	  middleIdx[2] = outIdx[2];
	  outIt.Set(imagePtr->GetPixel(outIdx) - middlePtr->GetPixel(middleIdx) );
     }


     // sum of square normalization.
     double ss = 0; // sum of square
     middlePtr->FillBuffer( 0 );
     for (middleIt.GoToBegin(); !middleIt.IsAtEnd(); ++ middleIt)
     {
	  middleIdx = middleIt.GetIndex();
	  imageIdx[0] = middleIdx[0];
	  imageIdx[1] = middleIdx[1];
	  imageIdx[2] = middleIdx[2];

	  ss = 0;
	  for (unsigned tsIdx = 0; tsIdx < tsLength; tsIdx ++) {
	       imageIdx[3] = tsIdx;
	       ss += pow(imagePtr->GetPixel(imageIdx), 2);
	  }
	  middleIt.Set(ss);
     }

     for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++ outIt)
     {
	  outIdx = outIt.GetIndex();
	  middleIdx[0] = outIdx[0];
	  middleIdx[1] = outIdx[1];
	  middleIdx[2] = outIdx[2];
	  if (middlePtr->GetPixel(middleIdx) > 0 ) {
	       outIt.Set(outIt.Get() / sqrt(middlePtr->GetPixel(middleIdx) ) );
	  }
	  else {
	       // must be outside of the mask.
	       outIt.Set( 0 );
	  }
     }


     // write back.
     WriterType4DFloat::Pointer writer = WriterType4DFloat::New();
     
     writer->SetInput(outPtr);
     writer->SetFileName(outfile);
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

     std::cout << "projectfmri(): File " << outfile << " saved.\n";
     return 0;
     }




