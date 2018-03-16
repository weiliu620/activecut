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

     std::string observedimage, trueimage, maskimage;
     std::string getMean;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "produce help message")
	  ("timepoints,p", po::value<unsigned>(&timeSeriesLength)->default_value(300), "Time series length.")	  
	  ("kappa", po::value<float>(&kappa)->default_value(50), "concentration parameter of von Mise. i.e. kappa.")	  

	  ("numclusters,k", po::value<unsigned>(&numClusters)->default_value(4), "Number of clusters. ")	  
	  ("obs", po::value<std::string>(&observedimage)->default_value("observedimage.nii"), "output noise image file")
	  ("true,t", po::value<std::string>(&trueimage)->default_value("trueimage.nii"), "input true image file")
	  ("getmean,g", po::value<std::string>(&getMean)->default_value("auto"), "how to construct mean vector: manual, auto.")
//	  ("mask,m", po::value<std::string>(&maskimage)->default_value("maskimage.nii"), "input maskimage file")
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
	  if (vm.count("help")) {
	       std::cout << "Usage: generateimage [options]\n";
	       cout << cmdline_options << "\n";
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

     // Allocate memory for fmri image.
     ImageType4DFloat::Pointer fmriPtr = ImageType4DFloat::New();
     ImageType4DFloat::SizeType fmriSize;
     ImageType4DFloat::IndexType fmriIdx;
     fmriSize[0] = trueImageSize[0];
     fmriSize[1] = trueImageSize[1];
     fmriSize[2] = trueImageSize[2];
     fmriSize[3] = timeSeriesLength;

     fmriIdx[0] = 0;
     fmriIdx[1] = 0;
     fmriIdx[2] = 0;
     fmriIdx[3] = 0;

     ImageType4DFloat::RegionType fmriRegion;
     fmriRegion.SetSize(fmriSize);
     fmriRegion.SetIndex(fmriIdx);
     fmriPtr->SetRegions(fmriRegion);
     fmriPtr->Allocate();


     vnl_vector<float> mu(timeSeriesLength), 
     	  northPole(timeSeriesLength), 
     	  vmfSample(timeSeriesLength);

     std::vector< vnl_matrix<float> > allRotationMat(numClusters);
     std::vector< vnl_vector<float> > allMu(numClusters);
     northPole.fill(0);
     northPole[timeSeriesLength - 1] = 1;


     short clsIdx = 0;

     // construct mean direction (mu) vector.
     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  allMu.at(clsIdx).set_size(timeSeriesLength);
	  allMu.at(clsIdx).fill(0);
	  if (getMean.compare("auto") == 0 ) {
	       // construct mean direction vectors. For k'th cluster center
	       // vector, let its k'th component be one, and others be
	       // zero.
	       if (clsIdx < timeSeriesLength) {
		    allMu.at(clsIdx)[clsIdx] = 1;
	       }
	       else {
		    allMu.at(clsIdx)[std::fmod(float(clsIdx), float(timeSeriesLength))] = -1;
	       }
	  }
	  else if (getMean.compare("manual") == 0 && timeSeriesLength == 3) {
	       printf("Give x and y value of cluster [%i]'s mean (separate by space):\n", clsIdx);
	       cin >> allMu.at(clsIdx)[0] >> allMu.at(clsIdx)[1];
//	       scanf ("%f %f\n", &allMu.at(clsIdx)[0], &allMu.at(clsIdx)[1]);
	       allMu.at(clsIdx)[2] = sqrt(1 - allMu.at(clsIdx)[0] * allMu.at(clsIdx)[0] - allMu.at(clsIdx)[1] * allMu.at(clsIdx)[1]);
	  }
	  else {
	       printf("To construct mean, choose auto for > 3 length time series, and choose manual for 3D space visualizaiton. exit. ");
	       exit (1);
	  }

	  printVnlVector(allMu.at(clsIdx), 100);

	  // Get rotation matrix that transform northPole to current
	  // mu.
	  allRotationMat.at(clsIdx).set_size(timeSeriesLength, timeSeriesLength);
	  getRotationMat(allRotationMat.at(clsIdx), northPole, allMu.at(clsIdx));
     }

     // sampling from von Mises-Fisher     
     for (trueImageIt.GoToBegin(); !trueImageIt.IsAtEnd(); ++ trueImageIt)
     {
     	  clsIdx = trueImageIt.Get();
     	  if(clsIdx >= 0) {
     	       generateVonMiseWood(vmfSample, kappa);    
     	       vmfSample = allRotationMat.at(clsIdx) * vmfSample;	       
	       if (getMean.compare("manual") == 0 && timeSeriesLength == 3 && verbose >= 1) {
		    printf("%u ", clsIdx);
		    printVnlVector(vmfSample, 100);
	       }
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
