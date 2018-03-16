#include <common.h>
#include <utility.h>

using namespace lemon;
extern twister_base_gen_type mygenerator;

int SamplingSub(ImageType3DChar::Pointer grpPtr,
		     ImageType3DChar::Pointer subPtr,
		     ImageType3DChar::Pointer maskPtr,
		     ParStruct & par);

int main(int argc, char* argv[])
{
     ParStruct par;
     unsigned seed = 0;
     std::string grpLabelFile, allSampleFile, initLabelFile;

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Take group label map as input, generate only one sample of sub label map. The initial sub label map can be either random map or a given image. Without the '--initlabel' argument, the sub label map will be init'd with random label. With the argument, it will be init'd with an input image. ")

	  ("burnin,b", po::value<unsigned>(&par.burnin)->default_value(20),
	   "number of scans for burn-in period. ")
	  ("inittemp", po::value<float>(&par.initTemp)->default_value(1),
	   "Initial temperature for annealing.")
	  ("finaltemp", po::value<float>(&par.finalTemp)->default_value(1),
	   "Final temperature for annealing to stop")

	  ("alpha", po::value<float>(&par.alpha)->default_value(0.7),
	   "connection between group lebel and individual subject level.")

	  ("beta", po::value<float>(&par.beta)->default_value(0.5),
	   "pairwise interation term")

	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")

	  ("numClusters,k", po::value<unsigned>(&par.numClusters)->default_value(7),
	   "Number of clusters. Default is 6.")

	   ("grouplabel,i", po::value<std::string>(&grpLabelFile)->default_value("grouplabel.nii.gz"), 
	    "Group level label map (Intensity value 1-K. Also used as mask file.")

	   ("initlabel,t", po::value<std::string>(&initLabelFile)->default_value("initlabel.nii.gz"), 
	    "Initial subject label map.")

	   ("out,o", po::value<std::string>(&allSampleFile)->default_value("testsamples.nii.gz"), 
	    "Output all test subject's sample label map, saved in a 4d file.")

	  ("verbose,v", po::value<unsigned short>(&par.verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: gentestlabel [options]\n";
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


     // read in initial grp label. Also used as mask. 1-based label converted to
     // 0-based label.
     ReaderType3DChar::Pointer grpReader = ReaderType3DChar::New();
     grpReader->SetFileName(grpLabelFile);
     grpReader->Update();
     // original 1-based group map only used for mask.
     ImageType3DChar::Pointer maskPtr = grpReader->GetOutput();
     ImageType3DChar::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());

     typedef itk::AddImageFilter <ImageType3DChar, ImageType3DChar, ImageType3DChar> AddImageFilterType;
     AddImageFilterType::Pointer addImageFilter = AddImageFilterType::New();
     addImageFilter->SetInput( maskPtr );
     addImageFilter->SetConstant2(-1);
     addImageFilter->Update();

     ImageType3DChar::Pointer grpPtr = addImageFilter->GetOutput();
     ImageType3DChar::SizeType grpSize = grpPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DChar grpIt(grpPtr, grpPtr->GetLargestPossibleRegion());

     // define a running sample.
     ImageType3DChar::IndexType start;
     start.Fill(0);
     ImageType3DChar::RegionType subRegion;
     subRegion.SetSize(grpSize);
     subRegion.SetIndex(start);
     ImageType3DChar::Pointer subPtr = ImageType3DChar::New();
     subPtr->SetRegions(subRegion);
     subPtr->Allocate();
     IteratorType3DChar subIt(subPtr, subPtr->GetLargestPossibleRegion());

     // init the running sample to the group map.
     if (vm["initlabel"].defaulted() ) {
	  // Uniform integer generator.
	  boost::uniform_int<> uni_int(0, par.numClusters - 1); // Uniform distribution of integers.
	  boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);
	  for (grpIt.GoToBegin(), subIt.GoToBegin(), maskIt.GoToBegin(); !grpIt.IsAtEnd(); ++ grpIt, ++ subIt, ++ maskIt) {
	       if (maskIt.Get() > 0) {
		    subIt.Set(roll_die() );
	       }
	       else {
		    subIt.Set( -1 );
	       }
	  }
     }
     else {
	  // specifically init sub map with a input image. 
	  ReaderType3DChar::Pointer initReader = ReaderType3DChar::New();
	  initReader->SetFileName(initLabelFile);
	  addImageFilter->SetInput( initReader->GetOutput()  );
	  addImageFilter->SetConstant2(-1);
	  addImageFilter->Update();
	  ImageType3DChar::Pointer initPtr = addImageFilter->GetOutput();
	  IteratorType3DChar initIt(initPtr, initPtr->GetLargestPossibleRegion());
	  for (initIt.GoToBegin(), subIt.GoToBegin(); !initIt.IsAtEnd(); ++ initIt, ++ subIt) {
	       subIt.Set( initIt.Get() );
	  }
     }
	  
     // sampling.
     for (unsigned bid = 0; bid < par.burnin; bid ++) {
     	  SamplingSub(grpPtr, subPtr, maskPtr, par);
	  if (par.verbose >= 1 && bid%500 == 0) {
	       printf("burnin sampling %i\n", bid);
	  }
     }
     
     save3dcharInc(subPtr, allSampleFile);
}

