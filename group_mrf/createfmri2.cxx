#include "commonalt.h"
#include <utilalt.h>
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
     unsigned numClusters = 0;
     std::string observedimage, trueimage, meantsfile;
     std::string noiseSrc, noiseFile;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "given label map and mean time series file, generate time series. The observed time series given cluster label is then by adding a gaussian noise (tsvar) on the mean time series.")
	  ("noisesource,n", po::value<std::string>(&noiseSrc)->default_value("internal"), "source of additive noise time series added on the mean time series. Can be internal or external. External means the noise can be spatially dependent (Gaussian blurred noise")
	  ("noisefile,f", po::value<std::string>(&noiseFile)->default_value("noise.nii.gz"), "Additive noise file.")

	  ("tsvar,s", po::value<float>(&tsvar)->default_value(0.5), "The  variance of noise added on mean time series")	  
	  ("obs", po::value<std::string>(&observedimage)->default_value("observedimage.nii"), "output noise image file")
	  ("true,t", po::value<std::string>(&trueimage)->default_value("trueimage.nii"), "input true image file")
	  ("timepoints,p", po::value<unsigned>(&timeSeriesLength)->default_value(300), "Time series length.")	
	  ("meants,m", po::value<std::string>(&meantsfile)->default_value("meants.txt"), "Mean time series file. Each row is a time series.")
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

     // if reading noise from file
     ImageType4DFloat::Pointer noisePtr;
     if (noiseSrc.compare("external") == 0) {
	  // read in noise data.
	  ReaderType4DFloat::Pointer noiseReader = ReaderType4DFloat::New();
	  noiseReader->SetFileName(noiseFile);
	  noiseReader->Update();
	  noisePtr = noiseReader->GetOutput();
	  ImageType4DFloat::SizeType noiseSize = noisePtr->GetLargestPossibleRegion().GetSize();
	  ImageType4DFloat::IndexType noiseIdx;    
     }	  

     // compute number of components.
     typedef itk::MinimumMaximumImageCalculator <ImageType3DChar>
	  ImageCalculatorFilterType;

     ImageCalculatorFilterType::Pointer imageCalculatorFilter
	  = ImageCalculatorFilterType::New ();

     imageCalculatorFilter->SetImage(trueImagePtr);
     imageCalculatorFilter->ComputeMaximum();
     numClusters = imageCalculatorFilter->GetMaximum();

     if (verbose >= 1) {
	  printf("numClusters = %i\n", numClusters);
     }


     // Allocate memory for mean time series.
     std::vector< vnl_vector<float> > allMu(numClusters);
     for (unsigned clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  allMu.at(clsIdx).set_size(timeSeriesLength);
	  allMu.at(clsIdx).fill(0);
     }

     // read mean time series file.
     std::string tsValue;
     std::ifstream meantsStream(meantsfile.c_str() );
     if (meantsStream.is_open() ) {
	  for (unsigned clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	       for (unsigned timeIdx = 0; timeIdx < timeSeriesLength; timeIdx ++) {
		    meantsStream >> tsValue;
		    allMu[clsIdx][timeIdx] = boost::lexical_cast<float> (tsValue);

	       }
	       // printVnlVector(allMu[clsIdx], 5);
	  }
     }
     else {
	  std::cout << "can not open file "  << meantsfile << std::endl;
     }




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


     // normal distribution.
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

		    if (noiseSrc.compare("internal") == 0) {
			 fmriPtr->SetPixel(fmriIdx, allMu[clsIdx][tsIdx] + nor());
		    }
		    else if (noiseSrc.compare("external") == 0) {
			 fmriPtr->SetPixel(fmriIdx, allMu[clsIdx][tsIdx] + noisePtr->GetPixel(fmriIdx) );
		    }
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

