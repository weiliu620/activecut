#include <common.h>

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
     std::string mask_file, input_image_file, out_file;
     double out_min = 0, out_max = 0;
     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "SLIC superpixel segmentation.")


	  ("outmin,i", po::value<double>(&out_min)->default_value(0),
	   "output image minimal value")

	  ("outmax,x", po::value<double>(&out_max)->default_value(255),
	   "output image max value")

	  ("data,d", po::value<std::string>(&input_image_file),
	   "Input  file. A 3D gipl or nii or nii.gz file.")

	  ("mask,m", po::value<std::string>(&mask_file),
	   "mask file")

	  ("out,o", po::value<std::string>(&out_file),
	   "Output rescaled file")
	  ;

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: grabcut [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }    


     // read data images.
     ReaderType3DF::Pointer dataReader = ReaderType3DF::New();
     dataReader->SetFileName(input_image_file);

     // // read mask images.
     // ReaderType3DI::Pointer maskReader = ReaderType3DI::New();
     // maskReader->SetFileName(mask_file);
     // maskReader->Update();
     // ImageType3DI::Pointer maskPtr = maskReader->GetOutput();

     // ImageType3DI::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     // IteratorType3DI maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());
     // ImageType3DI::IndexType maskIdx;


     typedef itk::RescaleIntensityImageFilter< ImageType3DF, ImageType3DF > RescaleFilterType;
     RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
     rescaleFilter->SetInput(dataReader->GetOutput());
     rescaleFilter->SetOutputMinimum(out_min);
     rescaleFilter->SetOutputMaximum(out_max);
     rescaleFilter->Update();
     std::cout << "in max: " << rescaleFilter->GetInputMaximum() << "\n";
     std::cout << "out max: " << rescaleFilter->GetOutputMaximum() << "\n";


     WriterType3DF::Pointer writer = WriterType3DF::New();
	  
     writer->SetInput(rescaleFilter->GetOutput() );
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

     std::cout << "rescale(): File " << out_file << " saved.\n";

     return 0;

     
}
