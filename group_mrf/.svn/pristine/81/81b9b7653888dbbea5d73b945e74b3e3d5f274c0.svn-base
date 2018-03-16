#include "commonalt.h"
using std::cout;
using std::cin;
using std::string;

int main (int argc, char* argv[])
{
     unsigned verbose = 0;
     float kappa;
     unsigned timeSeriesLength = 0;
     unsigned numClusters = 0;
     std::string observedimage, trueimage, maskimage, meantsfile;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "given label map and mean time series, generate time series in von Mises-Fisher distribution.")
	  ("kappa,k", po::value<float>(&kappa)->default_value(50), "concentration parameter of von Mise. i.e. kappa.")	  
	  ("obs", po::value<std::string>(&observedimage)->default_value("observedimage.nii"), "output noise image file")
	  ("true,t", po::value<std::string>(&trueimage)->default_value("trueimage.nii"), "input true image file")
	  ("meants,m", po::value<std::string>(&meantsfile)->default_value("meants.txt"), "Mean time series file. Each row is a time series.")
	  ("timepoints,p", po::value<unsigned>(&timeSeriesLength)->default_value(300), "Time series length.")	  
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

	  
     vnl_vector<float> mu(timeSeriesLength), 
	  northPole(timeSeriesLength), 
	  vmfSample(timeSeriesLength);
     northPole.fill(0);
     northPole[timeSeriesLength - 1] = 1;
     short clsIdx = 0;

     // read mean time series file.
     std::string tsValue;
     std::ifstream meantsStream(meantsfile.c_str() );
     if (meantsStream.is_open() ) {
	  for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	       for (unsigned timeIdx = 0; timeIdx < timeSeriesLength; timeIdx ++) {
		    meantsStream >> tsValue;
		    allMu[clsIdx][timeIdx] = boost::lexical_cast<float> (tsValue);

	       }
	  }
     }
     else {
	  std::cout << "can not open file "  << meantsfile << std::endl;
     }

     // Get rotation matrix that transform northPole to current
     // mu.
     std::vector< vnl_matrix<float> > allRotationMat(numClusters);
     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  allRotationMat.at(clsIdx).set_size(timeSeriesLength, timeSeriesLength);
	  getRotationMat(allRotationMat.at(clsIdx), northPole, allMu.at(clsIdx));
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


     // sampling from von Mises-Fisher     
     for (trueImageIt.GoToBegin(); !trueImageIt.IsAtEnd(); ++ trueImageIt)
     {
     	  clsIdx = trueImageIt.Get() - 1;
     	  if(clsIdx >= 0) {
     	       generateVonMiseWood(vmfSample, kappa);    
     	       vmfSample = allRotationMat.at(clsIdx) * vmfSample;	       

     	       // save samples in the fmri data structure.
     	       trueImageIdx = trueImageIt.GetIndex();
     	       fmriIdx[0] = trueImageIdx[0];
     	       fmriIdx[1] = trueImageIdx[1];
     	       fmriIdx[2] = trueImageIdx[2];
     	       for (fmriIdx[3] = 0; fmriIdx[3] < fmriSize[3]; fmriIdx[3] ++) {
     		    fmriPtr->SetPixel(fmriIdx, vmfSample[fmriIdx[3]]);
     	       } 
     	  } // end if
     } // end for

     
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
