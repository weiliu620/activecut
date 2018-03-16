#include <common.h>
#include <utility.h>
namespace po = boost::program_options;
int main(int argc, char* argv[])
{
     std::string srcFile, inFile, outFile;
     
     
     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "copy one file's header (voxel dimension, orientation) and overwrite on another file.")
	  ("source,s", po::value<std::string>(&srcFile)->default_value("srcfile.nii"),
	   "source file whose header information will be extracted.")
	  ("in,i", po::value<std::string>(&inFile)->default_value("infile.nii"),
	   "input file.")
	  ("out,o", po::value<std::string>(&outFile)->default_value("outfile.nii"),
	   "output file.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: gropumrf [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     ReaderType3DU::Pointer srcreader = ReaderType3DU::New();
     srcreader->SetFileName(srcFile);
     srcreader->Update();
     ImageType3DU::Pointer srcPtr = srcreader->GetOutput();
     
     ReaderType3DU::Pointer reader = ReaderType3DU::New();
     reader->SetFileName(inFile);
     reader->Update();
     ImageType3DU::Pointer inPtr = reader->GetOutput();

     inPtr->SetOrigin( srcPtr->GetOrigin() );
     inPtr->SetSpacing(srcPtr->GetSpacing() );
     inPtr->SetDirection(srcPtr->GetDirection() );

     WriterType3DU::Pointer writer = WriterType3DU::New();
	  
     writer->SetInput(inPtr);
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

     std::cout << "copyheader(): File " << outFile << " saved.\n";

     return 0;
}
