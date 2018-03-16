#include "commonalt.h"
#include "MCModel.h"

#include "itkImageSeriesReader.h"
#include <boost/filesystem.hpp>
typedef itk::ImageSeriesReader< ImageType4DFloat >  SeriesReaderType4DFloat;

int main (int argc, char* argv[])
{
     unsigned verbose = 0;
     std::string observedimage, trueimage, maskimage, meantsfile;
     std::string imagepath, varimage;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Read all nifti files (binary label) in a dir, and compute variance of the labels at each voxel, and save it to a file.")
	  ("labelimagepath,p", po::value<std::string>(&imagepath)->default_value("."),
	   "binary label image file path.")
	  ("varimage,a", po::value<std::string>(&varimage)->default_value("variance.nii.gz"), "variance of the image.")
	   ("mask,m", po::value<std::string>(&maskimage)->default_value("mask.nii.gz"),
	    "binary mask image. Only voxels within mask are computed.")
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


     boost::filesystem::path labelPathVar(imagepath);
     boost::filesystem::directory_iterator labelPathEnd;


     SeriesReaderType4DFloat::Pointer labelReader = SeriesReaderType4DFloat::New();

     for (boost::filesystem::directory_iterator labelPathIt(labelPathVar); labelPathIt != labelPathEnd; ++ labelPathIt) {

     	  labelReader->AddFileName( (*labelPathIt).path().string());	  
	  std::cout <<  "add " << (*labelPathIt).path().string() << "\n";
     }

     labelReader->Update();
     ImageType4DFloat::Pointer labelPtr = labelReader->GetOutput();
     ImageType4DFloat::SizeType labelSize = labelReader->GetOutput()->GetLargestPossibleRegion().GetSize();   
     ImageType4DFloat::IndexType labelIdx; 
     unsigned numSubs = labelSize[3];

     // read in mask image
     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(maskimage);
     maskReader->Update();
     ImageType3DChar::Pointer maskPtr = maskReader->GetOutput();

     ImageType3DChar::RegionType maskRegion;
     maskRegion = maskReader->GetOutput()->GetLargestPossibleRegion();
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     
     // create mean and variance image.
     ImageType3DFloat::Pointer varPtr = ImageType3DFloat::New();
     ImageType3DFloat::SizeType varSize;
     varSize[0] = labelSize[0];
     varSize[1] = labelSize[1];
     varSize[2] = labelSize[2];
     ImageType3DFloat::IndexType varIdx;
     varIdx.Fill(0);
     ImageType3DFloat::RegionType varRegion;
     varRegion.SetSize(varSize);
     varRegion.SetIndex(varIdx);
     varPtr->SetRegions( varRegion );
     varPtr->Allocate();
     varPtr->FillBuffer( 0 );

     IteratorType3DFloat varIt( varPtr, varPtr->GetLargestPossibleRegion() );
     float mean = 0, var = 0;
     unsigned numPoints = 0;
     float aveVar = 0;
     for (maskIt.GoToBegin(), varIt.GoToBegin(); !varIt.IsAtEnd(); ++ varIt, ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       varIdx = varIt.GetIndex();
	       labelIdx[0] = varIdx[0];
	       labelIdx[1] = varIdx[1];
	       labelIdx[2] = varIdx[2];

	       mean = 0;
	       for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
		    labelIdx[3] = subIdx;
		    mean += labelPtr->GetPixel(labelIdx) / numSubs;
	       }

	       var = 0;
	       for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
		    labelIdx[3] = subIdx;
		    var += pow (( labelPtr->GetPixel(labelIdx) - mean ), 2) / numSubs;
	       }
	       varPtr->SetPixel(varIdx, var);

	       aveVar += var;
	       numPoints ++;
	  } // maskIt
     }
     
     aveVar = aveVar / numPoints;


     // write var into a file
     WriterType3DFloat::Pointer writer = WriterType3DFloat::New();
     writer->SetInput(varPtr);
     writer->SetFileName(varimage);
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

     printf("Average variance: %f\n", aveVar);

     return 0;
}
     

     
     
