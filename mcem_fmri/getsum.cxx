#include "commonalt.h"
#include "MCModel.h"

#include "itkImageSeriesReader.h"
typedef itk::ImageSeriesReader< ImageType4DFloat >  SeriesReaderType4DFloat;

int main (int argc, char* argv[])
{
     unsigned verbose = 0;
     std::string imageFile;
     std::string imagepath, varimage;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Read all nifti files (binary label) in a dir, and compute variance of the labels at each voxel, and save it to a file.")
	  ("image,i", po::value<std::string>(&imageFile)->default_value("image.nii.gz"), "image to compute the sum of all voxels.")
	   ("verbose,v", po::value<unsigned int>(&verbose)->default_value(0), 
	    "Verbose level.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: generateimage [options]\n";
	       std::cout << mydesc << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }

     // read in image
     ReaderType3DFloat::Pointer imageReader = ReaderType3DFloat::New();
     imageReader->SetFileName(imageFile);
     imageReader->Update();
     ImageType3DFloat::Pointer imagePtr = imageReader->GetOutput();

     ImageType3DFloat::RegionType imageRegion;
     imageRegion = imageReader->GetOutput()->GetLargestPossibleRegion();
     ImageType3DFloat::SizeType imageSize = imageRegion.GetSize();
     IteratorType3DFloat imageIt(imagePtr, imagePtr->GetLargestPossibleRegion() );
     
     float sum = 0;
     for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++ imageIt) {
	  sum += imageIt.Get();
	  
	  if (verbose > 0) {
	       printf("voxel is: %f\n", imageIt.Get() );
	  }
     }
     
     printf("sum of voxels: %f", sum);
     return 0;
}
     

     
     
