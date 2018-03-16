#include <common.h>
#include <utility.h>
twister_base_gen_type mygenerator;

int main(int argc, char* argv[])
{
     unsigned seed = 0;
     double percent = 0;
     std::string init_file, mask_file, out_file;
     unsigned short verbose = 0;

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Given a unsmoothed label map, flip some labels. Used for generating subject label map given group map. ")

	  ("seed,s", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")

	  ("percent,p", po::value<double>(&percent)->default_value(0.1),
	   "Real number for the percentage of voxels being changed. 0 for no change. 1 for all change. ")

	   ("init,i", po::value<std::string>(&init_file)->default_value("init.nii.gz"), 
	    "initial label file.")

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



     // read in initial label file. Voxel value in [1, K].  Also act as mask
     // since outside brain is zero.
     typedef itk::ImageFileReader< ImageType3DU >  ReaderType3DU;
     ReaderType3DU::Pointer labelReader = ReaderType3DU::New();
     labelReader->SetFileName(init_file);
     labelReader->Update();
     ImageType3DU::Pointer labelPtr = labelReader->GetOutput();

     // get the number of clusters.
     typedef itk::MinimumMaximumImageCalculator <ImageType3DU> ImageCalculatorFilterType;
 
     ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New ();
     imageCalculatorFilter->SetImage(labelReader->GetOutput());
     imageCalculatorFilter->Compute();
     unsigned n_clusters = imageCalculatorFilter->GetMaximum();

     boost::random::uniform_real_distribution<> flip_generator(0, 1);
     boost::random::uniform_int_distribution<> uni_generator(1, n_clusters);

     IteratorType3DU labelIt(labelPtr, labelPtr->GetLargestPossibleRegion());
     for (labelIt.GoToBegin(); !labelIt.IsAtEnd(); ++ labelIt) {
	  if (labelIt.Get() > 0) {
	       // if the uniform sampled number is on the right side of the
	       // percent threshold, need to change label.
	       if (flip_generator(mygenerator) < percent) {

		    // just use 1-based labels.
		    labelIt.Set(uni_generator(mygenerator));
	       } // flip
	  } // > 0
     } // It

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


	       



