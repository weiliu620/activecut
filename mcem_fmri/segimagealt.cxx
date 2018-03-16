#include "commonalt.h"
#include "MCModel.h"
int ReadInitLabel(std::string initLabelFile, ImageType4DChar::Pointer samplePtr);

using std::cout;
using std::endl;

twister_base_gen_type mygenerator(42u);

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     float beta = 0.1;
     int burnin = 1;
     unsigned seed = 0;

     float ASE = 1;
     float ascentAlpha = 0;
     float ascentGamma = 0;
     float minllinc = 0;
     float addSamplePercent = 0;

     unsigned int numSamples = 1;
     unsigned int maxNumSamples = 1;
     unsigned int maxEmIter = 1;

     int numClusters = 5;
     bool acceptParEst = 1;
     bool convergeIndicator = 0;
     bool estBeta = 0;

     double meanKappa = 0;
     double stdKappa = 0;

     std::string observedFile, recoveredFile, maskFile, trueFile, initLabelFile;
     std::string parfile;
     std::string config_file;
     std::string samplingMethod;
     std::string estMethod;
     std::string cheat;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "produce help message")
	  ("samplingmethod,s", po::value<std::string>(&samplingMethod)->default_value("mc"),
	   "Sampling method in E step. Can be choosen from \"mc\", \"icm\", and \"gc\".")
	  ("estimattionmethod,e", po::value<std::string>(&estMethod)->default_value("icm"),
	   "Estimation method in E step. Can be choosen from \"ann\", \"icm\", and \"gc\".")
	  ("burnin,b", po::value<int>(&burnin)->default_value(20),
	   "number of scans for burn-in period. ")
	  ("numSamples,n", po::value<unsigned int>(&numSamples)->default_value(10),
	   "initial number of Monte Carlo samples. ")
	  ("maxnumscan", po::value<unsigned int>(&maxNumSamples)->default_value(50),
	   "max number of Monte Carlo samples. ")
	  ("maxemiter", po::value<unsigned int>(&maxEmIter)->default_value(100),
	   "max number of EM iteration. ")
	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("estbeta", po::value<bool>(&estBeta)->default_value(1),
	   "if estimate beta in MRF prior. choose 1 (for estimate) or 0 (no estimate).")
	  ("numClusters,k", po::value<int>(&numClusters)->default_value(4),
	       "Number of labels. Default is 5.")
	  ("meankappa", po::value<double>(&meanKappa)->default_value(300),
	   "Hyper-parameter for kappa: mean")
	  ("stdkappa", po::value<double>(&stdKappa)->default_value(0.1),
	   "Hyper-parameter for kappa: standard deviation")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")
	  ("input,i", po::value<std::string>(&observedFile)->default_value("observedimage.nii"),
	   "noised image file")
	  ("init", po::value<std::string>(&initLabelFile)->default_value("initlabel.nii"),
	   "Initial label map image.")
	  ("out,o", po::value<std::string>(&recoveredFile)->default_value("finalLabeled.nii"),
	   "output labeled image file")
	  ("trueimage,t", po::value<std::string>(&trueFile)->default_value("trueimage.nii"), "trueimage. Only for cheating...")
	  ("cheat", po::value<std::string>(&cheat)->default_value("no"),
	   "Init with true parameter values. \"yes\" or \"no\". ")
	  ("parfile", po::value<std::string>(&parfile)->default_value("parfile.txt"),
	   "file that contains true parameters. Only used when cheat = yes.");

     po::options_description cmdline_options;
     cmdline_options.add(generic);

     po::positional_options_description p;
     p.add("in", -1);

     po::variables_map vm;        
     po::store(po::command_line_parser(argc, argv).
	       options(cmdline_options).positional(p).run(), vm);
     po::notify(vm);    

     try {

	  if (vm.count("help") | argc == 1) {
	       std::cout << "Usage: segimagealt [options]\n";
	       std::cout << cmdline_options << "\n";
	       return 0;
	  }
	  
	  if (vm.count("maxNumSamples")) {
	       maxNumSamples = vm["maxNumSamples"].as<unsigned>();
	       if (maxNumSamples < numSamples) {
		    printf("maxNumSamples must be greater than numSamples.\n");
		    exit(1);
	       }
	  }
     }

     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // init random generator.
//     mygenerator.seed(static_cast<unsigned int>(std::time(0)));
     mygenerator.seed(static_cast<unsigned int>(seed));

     // read in image data.
     ReaderType4DFloat::Pointer fmriReader = ReaderType4DFloat::New();
     fmriReader->SetFileName(observedFile);
     fmriReader->Update();
     ImageType4DFloat::Pointer imagePtr = fmriReader->GetOutput();
     ImageType4DFloat::SizeType fmriSize = imagePtr->GetLargestPossibleRegion().GetSize();

     // // read in mask image
     // ReaderType3DFloat::Pointer maskReader = ReaderType3DFloat::New();
     // maskReader->SetFileName(maskFile);
     // maskReader->Update();
     // ImageType3DFloat::Pointer maskPtr = maskReader->GetOutput();

     // ImageType3DFloat::RegionType maskRegion;
     // maskRegion = maskReader->GetOutput()->GetLargestPossibleRegion();
     // ImageType3DFloat::SizeType maskSize = maskRegion.GetSize();
     // IteratorType3DFloat maskIt(maskPtr, maskRegion);
     // ImageType3DFloat::IndexType maskIdx;     

     // if (maskSize[0] != fmriSize[0] 
     // 	 || maskSize[1] != fmriSize[1] 
     // 	 || maskSize[2] != fmriSize[2]) {
     // 	  cout << "mask image and true label image have different size. Need double check before masking. Exit. " << endl;
     // }

     // Allocate memory for saving samples of MC.
     ImageType4DChar::Pointer samplePtr = ImageType4DChar::New();
     ImageType4DChar::IndexType start;
     start.Fill(0);

     ImageType4DChar::SizeType sampleSize;
     sampleSize[0] = fmriSize[0];
     sampleSize[1] = fmriSize[1];
     sampleSize[2] = fmriSize[2];

     if (samplingMethod.compare("mc") == 0 || samplingMethod.compare("nchain") == 0 ) {
	  sampleSize[3] = maxNumSamples; // number of MC samples.
     }
     else if (samplingMethod.compare("icm") == 0 || samplingMethod.compare("gc") == 0) {
	  sampleSize[3] = 1;
	  numSamples = 1;
     }
     else {
	  printf("sampling method (parameter \"sampling\" must be either \"mc\" or \"icm\" or \"gc\". \n");
	  exit (1);
     }
     ImageType4DChar::RegionType sampleRegion;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(start);
     samplePtr->SetRegions(sampleRegion);
     samplePtr->Allocate();
     samplePtr->FillBuffer(0);
     samplePtr->SetOrigin(imagePtr->GetOrigin());
     samplePtr->SetOrigin(imagePtr->GetOrigin());
     samplePtr->SetSpacing(imagePtr->GetSpacing());
     samplePtr->SetDirection(imagePtr->GetDirection());
     ImageType4DChar::PointType myorigin = samplePtr->GetOrigin();
     printf("samplePtr origion: [%f][%f][%f]\n", myorigin[0], myorigin[1], myorigin[2]);

     myorigin = imagePtr->GetOrigin();
     printf("samplePtr origion: [%f][%f][%f]\n", myorigin[0], myorigin[1], myorigin[2]);


     // read initial label map, and also naturally set -1 to outside of the mask.
     ReadInitLabel(initLabelFile, samplePtr);
     saveExtractVolume(samplePtr, sampleSize[3] - 1, recoveredFile);

     MCModel mcmodel(imagePtr,
		     samplePtr,
		     numClusters,
		     numSamples,
		     meanKappa,
		     stdKappa,
		     samplingMethod,
		     verbose);


     printf("mcmodel object init done.\n");
     mcmodel.printSelf("normal");

     // init with true parameters value?
     if (cheat.compare("yes") == 0) {
     	  mcmodel.initWithTruth(parfile, trueFile, samplePtr, imagePtr);
     	  printf("Force parameters init'd to true value.\n");
     	  mcmodel.printSelf("normal");
	       
     }
     else if (cheat.compare("no") == 0) {
     	  // doning nothing, since mcmodel is already initialized by constructor.
     }
     else {
     	  printf(" \"cheat\" parameter must be init'd with either yes or no. exit.");
     	  exit(1);
     }
	  
     // EM iteration begins.
     int emIter = 0;

     do {
     	  do {
     	       printf("\nEM iteration %d begin.\n", emIter);

     	       // Expectation step.
     	       if (samplingMethod.compare("mc") == 0) {
     		    printf("Expectation step. Monte Carlo Method, one long chain:\n");
		    mcmodel.estep(samplePtr, imagePtr, burnin);
     	       }
     	       else if(samplingMethod.compare("nchain") == 0) {
     		    printf("Expectation step. Monte Carlo, N indepedent chain:\n");
		    if (emIter == 0) {
			 mcmodel.nChainInit(samplePtr, "ml");
		    }
		    else {
			 mcmodel.nChainInit(samplePtr, "largestll");
		    }


     		    mcmodel.nChainSampling(samplePtr, imagePtr, burnin);	       

     	       }
     	       else if (samplingMethod.compare("icm") == 0) {
     		    printf("Expectation step. ICM Method:\n");
     		    mcmodel.icmEstep(samplePtr, imagePtr, "lastone", burnin);
     	       }
     	       else if (samplingMethod.compare("gc") == 0) {
     		    // mcmodel.estimateLabelsGraphcuts(imagePtr, samplePtr);
     	       }
     	       else {
     		    // do nothing. No need for safety check.
     	       }

     	       printf("After E step sampling:\n");
     	       mcmodel.printSelf("normal");


     	       // Maximization Step.
     	       printf("\nEM iteration %d Maximization step:\n", emIter);

     	       // Estimate likelihood parameters: mu and sigma.
     	       mcmodel.estimateMuSigmaWrapper(samplePtr, imagePtr, meanKappa, stdKappa);
     
	       // Estimate prior distribution parameters: beta.
	       if (estBeta) {
		    if (emIter == 0) {
			 mcmodel.estimatePriorParameter(0.01, 5, samplePtr, imagePtr);
		    }
		    else {
			 mcmodel.estPriorPramNewton(samplePtr, imagePtr,  5);
		    }
	       }
	       // mcmodel.estPriorPramNewton(samplePtr, imagePtr,  5);

     	       printf("After parameter estimation:\n");
     	       mcmodel.printSelf("normal");

     	       if (samplingMethod.compare("nchain") == 0) {
     		    // Compute asymptotic std error of delta Q.
     		    // mcmodel.updateSingleSampleLL(imagePtr, samplePtr);
     		    // ASE = mcmodel.computeASE(imagePtr, samplePtr);
     		    // printf("ASE = %f\n", ASE);
     	       }

     	       // print parameters for analsys.
     	       if (verbose >= 1) {
     		    printf("TABLELOG %i ", emIter);
     		    mcmodel.printSelf("table");	  
     	       }
	       
     	       emIter++;

     	       // if previous parameter estimatioin is good by ascent based method.
     	       if (emIter > 1 && samplingMethod.compare("nchain") == 0) {
     		    // acceptParEst = mcmodel.ascentAcceptEstimation(samplePtr, ascentAlpha, ASE, addSamplePercent);
     	       }
     	       else {
     		    acceptParEst = 1;
     	       }
     	  }
     	  while(!acceptParEst);

     	  // With known parameters, estimate labels and save in last
     	  // slice of smaplePtr.


		    // DBG
		    saveimage4dchar(samplePtr, "samples.nii");
		    saveExtractVolume(samplePtr, 0, "onesample.nii");
		    // end of DBG

     	  if (estMethod.compare("gc") == 0) {
     	       // mcmodel.estimateLabelsGraphcuts(imagePtr, samplePtr);
     	  }
     	  else if (estMethod.compare("icm") == 0) {
     	       if (samplingMethod.compare("nchain") == 0) {
     		    mcmodel.icmEstep(samplePtr, imagePtr, "all", burnin);
     	       }
     	       else {
     		    mcmodel.icmEstep(samplePtr, imagePtr, "lastone", burnin);		     

     	       }
     	  }
     	  else if (estMethod.compare("ann") == 0) {
     	       // mcmodel.annealing(imagePtr, samplePtr, 20);	       
     	  }
     	  else {
     	       // Something must be wrong. 
     	  }


     	  if (emIter == 1) {
     	       convergeIndicator = 0;
     	  }
     	  else if (emIter > maxEmIter) {
     	       convergeIndicator = 1;
     	  }
     	  else if (samplingMethod.compare("nchain") == 0) {
     	       // convergeIndicator = mcmodel.ascentConvergeTest(samplePtr, ascentGamma, ASE, minllinc);
     	  }
     	  else {
     	       convergeIndicator = 0;
     	  }

	  // save final results
	  saveExtractVolume(samplePtr, sampleSize[3] - 1, recoveredFile);
     }
     while (!convergeIndicator);

     if (verbose >= 1) {
	  mcmodel.printMeanKappa();
     }
		 
     // // printparameters for Monte Carlo analsys.
     // printf("FINALTABLE %i ", emIter);
     // mcmodel.print("table");	  


}

// Init last sample with initial label map.
int ReadInitLabel(std::string initLabelFile, ImageType4DChar::Pointer samplePtr)
{

     // init label file
     ReaderType3DChar::Pointer initReader = ReaderType3DChar::New();
     initReader->SetFileName(initLabelFile);
     ImageType3DChar::Pointer initPtr = initReader->GetOutput();
     initPtr->Update();
     IteratorType3DChar initIt(initPtr, initPtr->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType initIdx;

     // samplePtr.
     // ImageType4DChar::RegionType sampleRegion;
     ImageType4DChar::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize() ;
     ImageType4DChar::IndexType sampleIdx;

     for (initIt.GoToBegin(); !initIt.IsAtEnd(); ++ initIt) {
	  initIdx = initIt.GetIndex();
	  sampleIdx[0] = initIdx[0];
	  sampleIdx[1] = initIdx[1];
	  sampleIdx[2] = initIdx[2];
	  
	  for (sampleIdx[3] = 0; sampleIdx[3] < sampleSize[3]; sampleIdx[3] ++) {
	       samplePtr->SetPixel(sampleIdx, initIt.Get() - 1);
	  }
     }
     
}
