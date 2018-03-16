#include <common.h>
#include <utility.h>

using namespace lemon;
extern twister_base_gen_type mygenerator;

int SaveRunningSample(ImageType3DChar::Pointer subPtr,
		      ImageType4DS::Pointer samplePtr,
		      ImageType3DChar::Pointer maskPtr,
			   unsigned mcid);

int SamplingSub(ImageType3DChar::Pointer grpPtr,
		     ImageType3DChar::Pointer subPtr,
		     ImageType3DChar::Pointer maskPtr,
		     ParStruct & par);

int main(int argc, char* argv[])
{
     ParStruct par;
     unsigned seed = 0;
     unsigned interval = 0;

     std::string grpLabelFile, allSampleFile;


     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "This routine take the group label map as input, and sample from P(test subject label | group label).")

	  ("burnin,b", po::value<unsigned>(&par.burnin)->default_value(20),
	   "number of scans for burn-in period. ")
	  ("skip,s", po::value<unsigned>(&interval)->default_value(1),
	   "number of samples skipped between save samples. 0 means samples will be saved consecutively.")
	  ("inittemp", po::value<float>(&par.initTemp)->default_value(1),
	   "Initial temperature for annealing.")
	  ("finaltemp", po::value<float>(&par.finalTemp)->default_value(1),
	   "Final temperature for annealing to stop")
	  ("numSamples,n", po::value<unsigned int>(&par.numSamples)->default_value(20),
	   "Number of Monte Carlo samples. ")

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
     for (grpIt.GoToBegin(), subIt.GoToBegin(); !grpIt.IsAtEnd(); ++ grpIt, ++ subIt) {
	  subIt.Set( grpIt.Get() );
     }


     // define sample that save all MC samples.
     ImageType4DS::IndexType sampleIdx;
     sampleIdx.Fill( 0 );
     ImageType4DS::SizeType sampleSize;
     sampleSize[0] = maskSize[0];
     sampleSize[1] = maskSize[1];
     sampleSize[2] = maskSize[2];
     sampleSize[3] = par.numSamples;
     ImageType4DS::RegionType sampleRegion;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(sampleIdx);
     ImageType4DS::Pointer samplePtr = ImageType4DS::New();
     samplePtr->SetRegions( sampleRegion );
     samplePtr->Allocate();
     samplePtr->FillBuffer( 0 );

     // sampling.
     for (unsigned bid = 0; bid < par.burnin; bid ++) {
     	  SamplingSub(grpPtr, subPtr, maskPtr, par);
     }
     
     for (unsigned mcid = 0; mcid < par.numSamples; mcid ++) {
     	  // skip interval samples. (at least sampling once)
     	  for (unsigned skipid = 0; skipid < interval; skipid ++) {
     	       SamplingSub(grpPtr, subPtr, maskPtr, par);
     	  }

     	  SaveRunningSample(subPtr, samplePtr, maskPtr, mcid);
     	  printf("sample %i saved.\n", mcid);
     }

     // save all samples in a 4d file
     WriterType4DS::Pointer writer = WriterType4DS::New();
     writer->SetInput(samplePtr);
     writer->SetFileName(allSampleFile);

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

     std::cout << "gentestlabel(): File " << allSampleFile << " saved.\n";
}

// save the running sample (one MC sample) to the MC sample pool.
int SaveRunningSample(ImageType3DChar::Pointer subPtr,
		      ImageType4DS::Pointer samplePtr,
		      ImageType3DChar::Pointer maskPtr,
		      unsigned mcid)
{
     // IteratorType3DChar subIt(subPtr, subPtr->GetLargestPossibleRegion());
     IteratorType3DCharIdx subIt(subPtr, subPtr->GetLargestPossibleRegion());
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());
     ImageType4DS::IndexType sampleIdx;
     ImageType3DChar::IndexType subIdx;

     for (subIt.GoToBegin(), maskIt.GoToBegin(); !subIt.IsAtEnd(); ++ subIt, ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       subIdx = subIt.GetIndex();
	       sampleIdx[0] = subIdx[0];
	       sampleIdx[1] = subIdx[1];
	       sampleIdx[2] = subIdx[2];
	       sampleIdx[3] = mcid;

	       samplePtr->SetPixel(sampleIdx, subIt.Get() + 1 );
	  }
     }
     return 0;
}


