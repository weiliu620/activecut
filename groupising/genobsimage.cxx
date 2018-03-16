#include "common_mcem.h"
using std::cout;
using std::string;

// Init random generator as a global object.
twister_base_gen_type mygenerator(42u);
int main (int argc, char* argv[])
{
     unsigned verbose = 0;
     float variance = 0.1;
     std::string obsimage, trueimage;
     unsigned seed = 0;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "create observed noise image.")
	  ("labelmap,l", po::value<std::string>(&trueimage)->default_value("trueimage.nii"), "label map with labels from 1 to K.")
	  ("variance,r", po::value<float>(&variance)->default_value(0.1), "variance of the Gaussian component for each label.")
	  ("obs,o", po::value<std::string>(&obsimage)->default_value("obsimage.nii"), "output observed noised image.")
	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("verbose,v", po::value<unsigned int>(&verbose)->default_value(0), 
	   "Verbose level.");

     
     po::options_description cmdline_options;
     cmdline_options.add(generic);

     po::positional_options_description p;
     p.add("in", -1);

     po::variables_map vm;        
     po::store(po::command_line_parser(argc, argv).
	       options(cmdline_options).positional(p).run(), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") || argc == 1) {
	       std::cout << "Create group label map.\n";
	       cout << cmdline_options << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  cout << e.what() << "\n";
	  return 1;
     }    


     // read in true image
     ReaderType3DChar::Pointer trueReader = ReaderType3DChar::New();
     trueReader->SetFileName(trueimage);
     trueReader->Update();

     // we convert 1-based label to 0-based label here.
     AddConstantToImageFilterType::Pointer minusoneFilter = AddConstantToImageFilterType::New();
     minusoneFilter->SetInput(trueReader->GetOutput() );
     minusoneFilter->SetConstant(-1);
     minusoneFilter->Update();
     ImageType3DChar::Pointer truePtr = minusoneFilter->GetOutput();

     ImageType3DChar::RegionType trueRegion;
     trueRegion = trueReader->GetOutput()->GetLargestPossibleRegion();
     ImageType3DChar::SizeType trueSize = trueRegion.GetSize();
     IteratorType3DChar trueIt(truePtr, trueRegion);
     ImageType3DChar::IndexType trueIdx;     

     // Allocate memory for observed images.
     ImageType3DFloat::Pointer obsPtr = ImageType3DFloat::New();
     ImageType3DFloat::IndexType start;
     start[0] = 0;
     start[1] = 0;
     start[2] = 0;

     ImageType3DFloat::SizeType obsSize; 
     obsSize = trueSize;

     ImageType3DFloat::RegionType obsRegion;
     obsRegion.SetSize(obsSize);
     obsRegion.SetIndex(start);
     obsPtr->SetRegions(obsRegion);
     obsPtr->Allocate();
     obsPtr->FillBuffer( -1 );
     
     IteratorType3DFloat obsIt(obsPtr, obsPtr->GetLargestPossibleRegion() );

     mygenerator.seed(static_cast<unsigned int>(seed));     

     // normal distribution generator.
     boost::normal_distribution<> normal_dist(0, 1); // Normal distribution.
     boost::variate_generator<twister_base_gen_type&, boost::normal_distribution<> > nor(mygenerator, normal_dist);

     // Add noise.
     for (obsIt.GoToBegin(), trueIt.GoToBegin(); !obsIt.IsAtEnd(); ++obsIt, ++ trueIt) {
	  if ( trueIt.Get() >= 0) {
	       obsIt.Set( trueIt.Get() + nor() * variance);
	  }
     }


     // write back. to 1-based image.
     typedef itk::AddConstantToImageFilter <ImageType3DFloat, float, ImageType3DFloat> plusoneFilterType;

     plusoneFilterType::Pointer myfilter = plusoneFilterType::New();
     myfilter->SetInput(obsPtr);
     myfilter->SetConstant(1);
     myfilter->Update();
     
     WriterType3DFloat::Pointer writer = WriterType3DFloat::New();
     writer->SetInput(myfilter->GetOutput() );
     writer->SetFileName(obsimage);

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

     std::cout << "save3dchar(): File " << obsimage << " saved.\n";

     return 0;
}
