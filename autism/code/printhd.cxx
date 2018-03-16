#include <common.h>
#include <utility.h>
namespace po = boost::program_options;
struct CoType {
     double x;
     double y;
     double z;
};

int main(int argc, char* argv[])
{
     std::string input_image_file, out_file, mask_file, seeds_file;
     unsigned verbose = 0;
     double radius = 0;
     
     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "print nifti file header.")
	  ("data,d", po::value<std::string>(&input_image_file),
	   "input fmri images.")
	  ;

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: printhd [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // read data images.
     ReaderType3DU::Pointer dataReader = ReaderType3DU::New();
     dataReader->SetFileName(input_image_file);
     dataReader->Update();
     ImageType3DU::Pointer dataPtr = dataReader->GetOutput();
     ImageType3DU::SizeType dataSize = dataPtr->GetLargestPossibleRegion().GetSize();
     
     ImageType3DU::PointType mniIdx;
     mniIdx = dataPtr->GetOrigin();
     ImageType3DU::DirectionType dir = dataPtr->GetDirection();
     ImageType3DU::SpacingType spac = dataPtr->GetSpacing();
     
     printf("origin: [%.1f %.1f %.1f]\n", mniIdx[0], mniIdx[1], mniIdx[2]);
     printf("direction:\n  [%f %f %f]\n  [%f %f %f]\n  [%f %f %f]\n", dir(0,0), dir(0,1), dir(0,2), dir(1, 0), dir(1,1), dir(1,2), dir(2,0), dir(2,1), dir(2,2));
     printf("spacing: [%f %f %f]\n", spac[0], spac[1], spac[2]);

     itk::OrientImageFilter<ImageType3DU,ImageType3DU>::Pointer orienter =
	  itk::OrientImageFilter<ImageType3DU,ImageType3DU>::New();

     orienter->UseImageDirectionOn();
     orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAS);
     orienter->SetInput(dataPtr);
     orienter->Update();
     ImageType3DU::Pointer outPtr = orienter->GetOutput();

     mniIdx = outPtr->GetOrigin();
     dir = outPtr->GetDirection();
     spac = outPtr->GetSpacing();
     
     printf("origin: [%.1f %.1f %.1f]\n", mniIdx[0], mniIdx[1], mniIdx[2]);
     printf("direction:\n  [%f %f %f]\n  [%f %f %f]\n  [%f %f %f]\n", dir(0,0), dir(0,1), dir(0,2), dir(1, 0), dir(1,1), dir(1,2), dir(2,0), dir(2,1), dir(2,2));
     printf("spacing: [%f %f %f]\n", spac[0], spac[1], spac[2]);

}



