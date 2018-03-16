#include "commonalt.h"
#include "MCModel.h"

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     float beta = 0.1;
     int burnin = 1;

     float scale = 0;
     float offset = 0;
     float sigma = 2;
     float adjustRate = 1;
     float ASE = 1;
     float ascentAlpha = 0;
     float ascentGamma = 0;
     float minllinc = 0;
     float addSamplePercent = 0;

     unsigned int numScan = 1;
     unsigned int maxNumScan = 1;
     unsigned int maxEmIter = 1;
     int nlabels = 5;
     bool acceptParEst = 1;
     bool convergeIndicator = 0;


     std::string observedimage, recoveredimage;
     std::string config_file;
     std::string samplingMethod;
     std::string estMethod;
     std::string cheat;
     

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "produce help message")
	  ("beta", po::value<float>(&beta)->default_value(0.1), 
	   "MRF smoothness prior. When cheat = no, used as initnial value of EM. When cheat = yes, used as true value.")
	  ("sigma", po::value<float>(&sigma)->default_value(2),
	   "std deviation of noise signal. Used only for true value when cheat = yes.")
	  ("samplingmethod,m", po::value<std::string>(&samplingMethod)->default_value("mc"),
	   "Sampling method in E step. Can be choosen from \"mc\", \"icm\", and \"gc\".")
	  ("estimattionmethod,e", po::value<std::string>(&estMethod)->default_value("icm"),
	   "Sampling method in E step. Can be choosen from \"ann\", \"icm\", and \"gc\".")
	  ("burnin,b", po::value<int>(&burnin)->default_value(20),
	   "number of scans for burn-in period. ")
	  ("numScan,s", po::value<unsigned int>(&numScan)->default_value(20),
	   "initial number of Monte Carlo samples. ")
	  ("adjustrate", po::value<float>(&adjustRate)->default_value(1.0),
	   "When joint log-likelihood begins to decrease, increase number of Monte Carlo samples by multipling it by a constant adjustrate. By default adjustrate = 1, i.e. no adjustment.")
	  ("maxnumscan", po::value<unsigned int>(&maxNumScan)->default_value(50),
	   "max number of Monte Carlo samples. ")
	  ("maxemiter", po::value<unsigned int>(&maxEmIter)->default_value(100),
	   "max number of EM iteration. ")

	  ("nlabels,l", po::value<int>(&nlabels)->default_value(5),
	       "Number of labels. Default is 5.")
	  ("verbose", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")

	  ("in", po::value<std::string>(&observedimage)->default_value("observedimage.nii"),
	   "noised image file")
	  ("out", po::value<std::string>(&recoveredimage)->default_value("finalLabeled.nii"),
	   "output labeled image file")
	  ("scale", po::value<float>(&scale)->default_value(10),
	   "likelihood function mu = scale * label + offset. Only used when cheat = yes.")
	  ("offset", po::value<float>(&offset)->default_value(100),
	   "likelihood function mu = scale * label + offset. Only used when cheat = yes.")
	  ("cheat", po::value<std::string>(&cheat)->default_value("no"),
	   "Init with true parameter values. \"yes\" or \"no\". ")
	  ("ascentalpha", po::value<float>(&ascentAlpha)->default_value(0.05),
	   "Significant level in ascent-based method, to add more MC samples.")
	  ("ascentgamma", po::value<float>(&ascentGamma)->default_value(0.05),
	   "Significant level in ascent-based method, to stop EM.")
	  ("minllinc", po::value<float>(&minllinc)->default_value(5),
	   "Minimal increasement of joint log-likelihood for EM to continue.")
	  ("addsamplepercent", po::value<float>(&addSamplePercent)->default_value(0.5),
	   "Percent of samples to add if ascent based method does not accept parameters.");

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
	       std::cout << cmdline_options << "\n";
	       return 0;
	  }
	  
	  if (vm.count("maxNumScan")) {
	       maxNumScan = vm["maxNumScan"].as<unsigned>();
	       if (maxNumScan < numScan) {
		    printf("maxNumScan must be greater than numScan.\n");
		    exit(1);
	       }
	  }
     }

     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // read in image data.
     ReaderType::Pointer reader = ReaderType::New();
     reader->SetFileName(observedimage);
     reader->Update();
     ImageType2D::Pointer imagePtr = reader->GetOutput();
     ImageType2D::SizeType size =
          reader->GetOutput()->GetLargestPossibleRegion().GetSize();

     // Allocate memory for saving samples of MC.
     ImageType3D::Pointer samplePtr = ImageType3D::New();
     ImageType3D::IndexType start;
     start[0] = 0;
     start[1] = 0;
     start[2] = 0;
     ImageType3D::SizeType sampleSize;
     sampleSize[0] = size[0];
     sampleSize[1] = size[1];

     if (samplingMethod.compare("mc") == 0 || samplingMethod.compare("nchain") == 0 ) {
	  sampleSize[2] = maxNumScan; // number of MC samples.
     }
     else if (samplingMethod.compare("icm") == 0 || samplingMethod.compare("gc") == 0) {
	  sampleSize[2] = 1;
	  numScan = 1;
     }
     else {
	  printf("sampling method (parameter \"sampling\" must be either \"mc\" or \"icm\" or \"gc\". \n");
	  exit (1);
     }
     ImageType3D::RegionType region;
     region.SetSize(sampleSize);
     region.SetIndex(start);
     samplePtr->SetRegions(region);
     samplePtr->Allocate();

     MCModel mcmodel(imagePtr,
		     samplePtr,
		     nlabels,
		     numScan,
		     beta,
		     samplingMethod,
		     verbose);


     printf("mcmodel object init done.\n");
     mcmodel.print("normal");

     // init with true parameters value?
     if (cheat.compare("yes") == 0) {
	  mcmodel.initWithTruth(sigma, offset, scale, beta);
	  printf("Force parameters init'd to true value.\n");
	  mcmodel.print("normal");
	       
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
		    printf("Expectation step. Monte Carlo Method:\n");
		    mcmodel.estep(imagePtr, samplePtr, burnin);
	       }
	       else if(samplingMethod.compare("nchain") == 0) {
		    printf("Expectation step. Monte Carlo, N indepedent chain:\n");
		    mcmodel.nChainInit(samplePtr, "largestll");
		    mcmodel.nChainSampling(imagePtr, samplePtr, burnin);	       
	       }
	       else if (samplingMethod.compare("icm") == 0) {
		    printf("Expectation step. ICM Method:\n");
		    mcmodel.icmEstep(imagePtr, samplePtr, "lastone", burnin);
	       }
	       else if (samplingMethod.compare("gc") == 0) {
		    mcmodel.estimateLabelsGraphcuts(imagePtr, samplePtr);
	       }
	       else {
		    // do nothing. No need for safety check.
	       }

	       printf("After E step sampling:\n");
	       mcmodel.print("normal");
	       
	       // Maximization Step.
	       printf("\nEM iteration %d Maximization step:\n", emIter);

	       // Estimate likelihood parameters: mu and sigma.
	       mcmodel.estimateMu(samplePtr, imagePtr);
	       mcmodel.estimateSigma(samplePtr, imagePtr);
     
	       // Estimate prior distribution parameters: beta.
	       if (emIter == 0) {
		    mcmodel.estimatePriorParameter(0.1, 5, samplePtr);
	       }
	       else {
		    mcmodel.estPriorPramNewton(samplePtr, 5);
	       }
	  

	       saveExtractSlice(samplePtr, sampleSize[2]-1, recoveredimage);

	       printf("After parameter estimation:\n");
	       mcmodel.print("normal");

	       if (samplingMethod.compare("nchain") == 0) {
		    // Compute asymptotic std error of delta Q.
		    mcmodel.updateSingleSampleLL(imagePtr, samplePtr);
		    ASE = mcmodel.computeASE(imagePtr, samplePtr);
		    printf("ASE = %f\n", ASE);
	       }

	       // print parameters for analsys.
	       if (verbose >= 1) {
		    printf("TABLELOG %i ", emIter);
		    mcmodel.print("table");	  
	       }
	       
	       emIter++;

	       // if previous parameter estimatioin is good by ascent based method.
	       if (emIter > 1 
		   && (samplingMethod.compare("nchain") == 0) || (samplingMethod.compare("mc") == 0)) {
		    acceptParEst = mcmodel.ascentAcceptEstimation(samplePtr, ascentAlpha, ASE, addSamplePercent);
	       }
	       else {
		    acceptParEst = 1;
	       }
	  }
	  while(!acceptParEst);

	  // With known parameters, estimate labels and save in last
	  // slice of smaplePtr.

	  if (estMethod.compare("gc") == 0) {
	       mcmodel.estimateLabelsGraphcuts(imagePtr, samplePtr);
	  }
	  else if (estMethod.compare("icm") == 0) {
	       if (samplingMethod.compare("nchain") == 0) {
		    mcmodel.icmEstep(imagePtr, samplePtr, "all", burnin);
	       }
	       else {
		    mcmodel.icmEstep(imagePtr, samplePtr, "lastone", burnin);		     
	       }
	  }
	  else if (estMethod.compare("ann") == 0) {
	       mcmodel.annealing(imagePtr, samplePtr, 20);	       
	  }
	  else if (estMethod.compare("no") == 0) {
	       printf("No estimating labels after parameter estimation.\n");
	  }
	  else {
	       printf("estimation emthod incorrect. exit.\n");
	       exit (1);
	  }

	  saveimage3d(samplePtr, "sample.nii");

	  if (emIter == 1) {
	       convergeIndicator = 0;
	  }
	  else if (emIter > maxEmIter) {
	       convergeIndicator = 1;
	  }
	  else if (samplingMethod.compare("nchain") == 0) {
	       convergeIndicator = mcmodel.ascentConvergeTest(samplePtr, ascentGamma, ASE, minllinc);
	  }
	  else {
	       convergeIndicator = 0;
	  }
     }
     while (!convergeIndicator);

		 
     // print parameters for Monte Carlo analsys.
     printf("FINALTABLE %i ", emIter);
     mcmodel.print("table");	  
     saveExtractSlice(samplePtr, sampleSize[2]-1, recoveredimage);

}


