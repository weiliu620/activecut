#include <common.h>
#include <utility.h>
twister_base_gen_type mygenerator;

int main(int argc, char* argv[])
{
     unsigned seed = 0, n_cls;
     std::string out_file, mask_file;
     unsigned short verbose = 0;

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Generate radom sampled label map.")

	  ("seed,s", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")

	  ("clusters,k", po::value<unsigned int>(&n_cls)->default_value(4),
	   "Number of clusters/labels.")

	   ("mask,m", po::value<std::string>(&mask_file)->default_value("mask.nii.gz"), 
	    "mask file.")

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
     ImageType3DU::RegionType labelRegion = maskPtr->GetLargestPossibleRegion();
     ImageType3DU::Pointer labelPtr = ImageType3DU::New();
     labelPtr->SetRegions(labelRegion);
     labelPtr->Allocate();
     labelPtr->FillBuffer(0);

     boost::random::uniform_int_distribution<> uni_generator(1, n_cls);

     IteratorType3DU labelIt(labelPtr, labelPtr->GetLargestPossibleRegion());
     for (labelIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++labelIt) {
	  if (maskIt.Get() > 0) {
	       // just use 1-based labels.
	       labelIt.Set(uni_generator(mygenerator));
	  }
     }

     // save the label map. 
     typedef itk::ImageFileWriter< ImageType3DU >  WriterType3DU;
     WriterType3DU::Pointer writer = WriterType3DU::New();
     writer->SetInput( labelPtr );
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


	       



