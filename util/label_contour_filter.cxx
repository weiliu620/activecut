#include <common.h>
#include <util.h>
#include <itkLabelContourImageFilter.h>

int main(int argc, char* argv[])
{
     std::string label_file, out_file;
     unsigned short verbose = 0;
     bool full_connected = false;
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Apply the LabelContourImageFilter to the input label image and keep only the pixels on the outline of the regions.")
	  ("label,l", po::value<std::string>(&label_file)->default_value("label.nii.gz"), 
	   "label image.")
	  ("out,o", po::value<std::string>(&out_file)->default_value("out.nii.gz"), 
	   "output 2d overlaied map.")
	  ("full,f", po::value<bool>(&full_connected)->default_value(false),
	   "if the connected components are full connected.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: label_contour_filter [options]\n";
	       std::cout << mydesc << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // read in gray intensity and label file. 
     ReaderType3DS::Pointer labelReader = ReaderType3DS::New();
     labelReader->SetFileName(label_file);
     ImageType3DS::Pointer labelPtr = labelReader->GetOutput();
     labelPtr->Update();

     typedef itk::LabelContourImageFilter<ImageType3DS, ImageType3DS> LabelContourImageFilterType;
     LabelContourImageFilterType::Pointer labelContourImageFilter =
	  LabelContourImageFilterType::New();
     labelContourImageFilter->SetInput(labelPtr);
     labelContourImageFilter->SetFullyConnected(full_connected);
     labelContourImageFilter->Update();

     WriterType3DS::Pointer writer = WriterType3DS::New();
     writer->SetInput(labelContourImageFilter->GetOutput());
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

     return 0;
}
