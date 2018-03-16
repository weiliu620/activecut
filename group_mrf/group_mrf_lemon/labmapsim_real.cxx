#include <common.h>
#include <utility.h>
#include "itkDiscreteGaussianImageFilter.h"

twister_base_gen_type mygenerator;

int main(int argc, char* argv[])
{
     unsigned seed = 0, n_cls;
     std::string out_file, mask_file;
     unsigned short verbose = 0;
     double variance = 0;

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Generate radom map with real number in (0,1).")

	  ("seed,s", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")

	   ("mask,m", po::value<std::string>(&mask_file)->default_value("mask.nii.gz"), 
	    "mask file.")

	  ("var", po::value<double>(&variance)->default_value(1),
	   "Variance of the Gaussian filter")

	   ("out,o", po::value<std::string>(&out_file)->default_value("outlabel.nii.gz"), 
	    "Output label file. Outside mask is zero, inside mask has value [1, K]")

	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");


     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: genrandom [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // init random generator.
     mygenerator.seed(static_cast<unsigned int>(seed));

     // read mask file.
     typedef itk::ImageFileReader< ImageType3DU >  ReaderType3DU;
     ReaderType3DU::Pointer maskReader = ReaderType3DU::New();
     maskReader->SetFileName(mask_file);
     maskReader->Update();
     ImageType3DU::Pointer maskPtr = maskReader->GetOutput();
     IteratorType3DU maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());

     // create label file.
     ImageType3DFloat::RegionType labelRegion = maskPtr->GetLargestPossibleRegion();
     ImageType3DFloat::Pointer labelPtr = ImageType3DFloat::New();
     labelPtr->SetRegions(labelRegion);
     labelPtr->Allocate();
     labelPtr->FillBuffer(0);

     boost::random::uniform_real_distribution<> uni_generator(0, 1);

     IteratorType3DFloat labelIt(labelPtr, labelPtr->GetLargestPossibleRegion());
     for (labelIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++labelIt) {
	  if (maskIt.Get() > 0) {
	       labelIt.Set(uni_generator(mygenerator));
	  }
     }

     typedef itk::DiscreteGaussianImageFilter<
	  ImageType3DFloat, ImageType3DFloat >  filterType;
     filterType::Pointer gaussianFilter = filterType::New();

     gaussianFilter->SetInput( labelPtr);
     gaussianFilter->SetVariance(variance);


     // save the label map. 
     typedef itk::ImageFileWriter< ImageType3DFloat >  WriterType3DFloat;
     WriterType3DFloat::Pointer writer = WriterType3DFloat::New();
     writer->SetInput( gaussianFilter->GetOutput() );
     writer->SetFileName(out_file);


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
     std::cout << "labelmapsim(): File  " << out_file << " saved.\n";

}


	       



