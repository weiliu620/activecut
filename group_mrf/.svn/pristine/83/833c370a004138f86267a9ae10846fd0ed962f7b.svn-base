#include "commonalt.h"
using std::cout;
using std::cin;
using std::string;
twister_base_gen_type mygenerator(42u);

int main (int argc, char* argv[])
{
     unsigned verbose = 0;
     float tsvar = 0;
     unsigned seed = 0;
     unsigned timeSeriesLength = 0;
     std::string observedimage, trueimage;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Generate 4D noise time series file.")
	  ("tsvar,s", po::value<float>(&tsvar)->default_value(0.5), "The  variance of noise added on mean time series")	
	  ("true,t", po::value<std::string>(&trueimage)->default_value("trueimage.nii"), "input true image file. Only size information is used to create noise nifti file.")  
	  ("obs", po::value<std::string>(&observedimage)->default_value("observedimage.nii"), "output noise image file")
	  ("timepoints,p", po::value<unsigned>(&timeSeriesLength)->default_value(300), "Time series length.")	
	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")  
	  ("verbose,v", po::value<unsigned int>(&verbose)->default_value(0), 
	   "Verbose level.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: generateimage [options]\n";
	       cout << mydesc << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  cout << e.what() << "\n";
	  return 1;
     }    

     // init random generator.
     mygenerator.seed(static_cast<unsigned int>(seed));

     // read in true label image
     ReaderType3DChar::Pointer trueImageReader = ReaderType3DChar::New();
     trueImageReader->SetFileName(trueimage);
     trueImageReader->Update();
     ImageType3DChar::Pointer trueImagePtr = trueImageReader->GetOutput();

     ImageType3DChar::RegionType trueImageRegion;
     trueImageRegion = trueImageReader->GetOutput()->GetLargestPossibleRegion();
     ImageType3DChar::SizeType trueImageSize = trueImageRegion.GetSize();
     IteratorType3DChar trueImageIt(trueImagePtr, trueImageRegion);
     ImageType3DChar::IndexType trueImageIdx;

     // Allocate memory for fmri image.
     ImageType4DFloat::Pointer fmriPtr = ImageType4DFloat::New();
     ImageType4DFloat::SizeType fmriSize;
     ImageType4DFloat::IndexType fmriIdx;
     fmriSize[0] = trueImageSize[0];
     fmriSize[1] = trueImageSize[1];
     fmriSize[2] = trueImageSize[2];
     fmriSize[3] = timeSeriesLength;

     fmriIdx.Fill(0);
     ImageType4DFloat::RegionType fmriRegion;
     fmriRegion.SetSize(fmriSize);
     fmriRegion.SetIndex(fmriIdx);
     fmriPtr->SetRegions(fmriRegion);
     fmriPtr->Allocate();
     fmriPtr->FillBuffer( 0 );


     // Uniform real random generator.
     boost::normal_distribution<double> normal_dist(0, sqrt(tsvar)); // normal distribution
     boost::variate_generator<twister_base_gen_type&, boost::normal_distribution<> > nor(mygenerator, normal_dist);

     unsigned tsIdx = 0;
     short clsIdx = 0;
     for (trueImageIt.GoToBegin(); !trueImageIt.IsAtEnd(); ++ trueImageIt) {
     	  clsIdx = trueImageIt.Get() - 1;
     	  if(clsIdx >= 0) {

	       trueImageIdx = trueImageIt.GetIndex();
     	       fmriIdx[0] = trueImageIdx[0];
     	       fmriIdx[1] = trueImageIdx[1];
     	       fmriIdx[2] = trueImageIdx[2];
     	       for (fmriIdx[3] = 0; fmriIdx[3] < fmriSize[3]; fmriIdx[3] ++) {
		    tsIdx = fmriIdx[3];
		    fmriPtr->SetPixel(fmriIdx,  nor());
     	       } 
	  } // clsIdx >= 0
     }

     
     WriterType4DFloat::Pointer writer = WriterType4DFloat::New();
     writer->SetInput(fmriPtr);
     writer->SetFileName(observedimage);
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

     std::cout << "File " << observedimage << " saved.\n";
     return 0;
}
