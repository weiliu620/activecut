#include "common.h"

twister_base_gen_type mygenerator(42u);

int main(int argc, char* argv[])
{
     int _DBG = 0;
     float alpha = 0; 
     float beta = 0.7, beta_old = 0;
     int burnin = 5;
     int label = 0;
     float scale = 0;
     float offset = 0;

     int numScan = 10;
     float sigma = 2;
     int nlabels = 5, k = 0;

     string observedimage, recoveredimage;
     string config_file;
     
     double jll = 0; // joint log-likelihood log p(f, d) = log (f) + log (d|f);
     double priorll = 0, cll = 0;
     double relerr = 0;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help", "produce help message")
	  ("version,v", "print version string")
	  ("config,c", po::value<string>(&config_file)->default_value("mcem.cfg"),
	   "name of a file of a configuration.")
	  ("burnin,b", po::value<int>(&burnin)->default_value(10),
	   "number of scans for burn-in period. Default 10.")
	  ("observed", po::value<string>(), "output noise image file")
	  ("recovered", po::value<string>(), "output recovered image file");

     // Declare a group of options that will be allowed only in configuration file.
     po::options_description config("Configuration options (can be given at only config file");
     config.add_options()
	  ("alpha", po::value<float>(&alpha)->default_value(0),
	   "Set MRF parameter alpha, used as initial value of EM.  Default 0.")
	  ("beta", po::value<float>(&beta)->default_value(0.7), 
	   "Set MRF parameter beta, as initnial value of EM.  Default 0.7.")
	  ("scale", po::value<float>(&scale)->default_value(10),
	   "likelihood function mu = scale * label + offset. Default 10.")
	  ("offset", po::value<float>(&offset)->default_value(100),
	   "likelihood function mu = scale * label + offset. Default 100.")
	  ("nlabels", po::value<int>(&nlabels)->default_value(5),
	       "Number of labels. Default is 5.")
	  ("sigma", po::value<float>()->default_value(2), 
	   "Set noise stand deviation sigma. Used as initial value of EM. Default is 2.")
	  ("debuglevel", po::value<int>(&_DBG)->default_value(0), 
	   "set debug level. 0 for no debug output (default). 3 for most debug output.");

     po::options_description cmdline_options;
     cmdline_options.add(generic);

     po::options_description config_file_options;
     config_file_options.add(config);

     po::positional_options_description p;
     p.add("in", -1);

     po::variables_map vm;        
     po::store(po::command_line_parser(argc, argv).
	       options(cmdline_options).positional(p).run(), vm);
     po::notify(vm);    

     try {
	  ifstream ifs(config_file.c_str());
	  if (!ifs) {
	       cout << "can not open config file: " << config_file << "\n";
	       return 0;
	  } else
	  {
	       store(parse_config_file(ifs, config_file_options), vm);
	       notify(vm);
	  }

	  if (vm.count("help")) {
	       cout << "Usage: generateimage [options]\n";
	       cout << cmdline_options << "\n";
	       return 0;
	  }
     
	  if (vm.count("version")) {
	       cout << " version 1.0.\n";
	       return 0;
	  }

	  if (vm.count("observed")) {
	       observedimage = vm["observed"].as<string>();
	  } else {
	       observedimage = "observedimage.nii";

	  }
	  
	  if (vm.count("recovered")) {
	       recoveredimage = vm["recovered"].as<string>();
	  } else {
	       recoveredimage = "recoveredimage.nii";
	  }
     }

     catch(exception& e)
     {
	  cout << e.what() << "\n";
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
     sampleSize[2] = numScan; // number of MC samples.
     ImageType3D::RegionType region;
     region.SetSize(sampleSize);
     region.SetIndex(start);
     samplePtr->SetRegions(region);
     samplePtr->Allocate();

     // Init the cluster component parameters.
     CompType* cls  = new CompType[nlabels];
     CompType* cls_old  = new CompType[nlabels];
     for (label = 0; label < nlabels; label ++) {
	  
	  cls[label].mu = label * scale + offset;
	  cls[label].sigma = sigma;
     }

     // EM iteration begins.
     int emiter = 0;
     double maxStepLength = 0;
     do {
	  emiter++;
	  printf("EM iteration %d begin.\n", emiter);

	  // Expectation step.
	  printf("Expectation step:\n");
	  estep(imagePtr, samplePtr, cls, nlabels, burnin, numScan, beta);

	  // Maximization Step.
	  printf("Maximization step:\n");
	  // Estimate likelihood parameters: mu and sigma.
	  memcpy(cls_old, cls, sizeof(CompType) * nlabels);

	  est_ll(samplePtr, imagePtr, cls, nlabels);
     
	  // Estimate prior distribution parameters: beta.
	  //     ll_over_beta(samplePtr, nlabels);
	  maxStepLength = abs(beta_old - beta);
	  beta_old = beta;
	  beta =  golden_section_search(0.01, 5, samplePtr, nlabels);
//	  beta = descent(samplePtr, nlabels, beta, maxStepLength);

	  // Get new value.
	  priorll = eval_ll(samplePtr, nlabels, beta);
	  cll = eval_cll(samplePtr, imagePtr, cls);
	  jll = priorll + cll;
	  
	  // print parameters and log-likelihood value.
	  printcls(cls, nlabels);
	  printf("beta = %f, prior ll = %f, cll = %f, jll = %f\n", beta, priorll, cll, jll);
	  printf("EM iteration %d done.\n", emiter);

	  // Stop criteria.
	  relerr = 0;
	  for (k = 0; k < nlabels; k++) {
	       relerr = relerr + abs((cls[k].mu - cls_old[k].mu)/cls[k].mu);
	       relerr = relerr + abs((cls[k].sigma - cls_old[k].sigma)/cls[k].sigma);
	       }
	  relerr = relerr + abs((beta - beta_old)/beta);
	  saveimage3d(samplePtr, "mcSamples.nii");     
	  
	  
     }
     while (relerr > 0.001);

     saveimage3d(samplePtr, "mcSamples.nii");     
     delete[] cls, cls_old;

}

int est_ll(ImageType3D::Pointer samplePtr,
	   ImageType2D::Pointer imagePtr,
	   CompType* cls,
	   int numClusters)
{
     int k = 0; // cluster number.

     ImageType3D::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;     
     ImageType2D::IndexType idx;     
     
     // Compute number of points in each cluster.
     for (k = 0; k < numClusters; k ++) {
	  cls[k].numPoints = 0;
	  cls[k].mu = 0;
	  cls[k].sigma = 0;
	  
     }
     
     for (sampleIdx[0] = 0;sampleIdx[0] < sampleSize[0]; sampleIdx[0] ++) {
	  for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1] ++) {
	       for (sampleIdx[2] = 0; sampleIdx[2] < sampleSize[2]; sampleIdx[2] ++) {
		    idx[0] = sampleIdx[0];
		    idx[1] = sampleIdx[1];
		    k = samplePtr->GetPixel(sampleIdx);
		    cls[k].mu += imagePtr->GetPixel(idx);
		    cls[k].numPoints ++;
	       } 
	  } 
     }
     // compute mu_k
     for (k = 0; k < numClusters; k ++) {
	  cls[k].mu = cls[k].mu/cls[k].numPoints;
     }
     
     // Compute sigma_k.
     for (sampleIdx[0] = 0;sampleIdx[0] < sampleSize[0]; sampleIdx[0] ++) {
	  for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1] ++) {
	       for (sampleIdx[2] = 0; sampleIdx[2] < sampleSize[2]; sampleIdx[2] ++) {
		    idx[0] = sampleIdx[0];
		    idx[1] = sampleIdx[1];
		    k = samplePtr->GetPixel(sampleIdx);
		    cls[k].sigma += (imagePtr->GetPixel(idx) - cls[k].mu) * (imagePtr->GetPixel(idx) - cls[k].mu);
	       } 
	  } 
     }
     
     for (k = 0; k < numClusters; k ++) {
	  cls[k].sigma = sqrt(cls[k].sigma/cls[k].numPoints);
     }
     return 0;
}
	   
int ll_over_beta(ImageType3D::Pointer samplePtr,
	       int numClusters)
{
     // open a file to write Q function value over all beta.
     FILE* pFile;
     pFile = fopen("prior_LL_overbeta.txt", "w");

     for (double beta = 0.001; beta < 5; beta = beta + 0.2) {
	  
	  fprintf(pFile, "%f %f\n", beta, eval_ll(samplePtr, numClusters, beta));
     }

     fclose(pFile);
     return 0;
}

// Given a beta, return the prior log-likelihood.
double eval_ll(ImageType3D::Pointer samplePtr,
		  int numClusters,
		  double beta)
{
     int k = 0; // cluster number.
     int S = 0, curS = 0; //  neighbors sum.
     double sumexp = 0;
     double Q = 0; // log likelihood (prior).
     int thisLabel = 0;
     int curLabel = 0;

     ImageType3D::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;     
     ImageType3D::IndexType nidx;     

     Q = 0;

     for (sampleIdx[2] = 0; sampleIdx[2] < sampleSize[2]; sampleIdx[2] ++) {
	  nidx[2] = sampleIdx[2];
	  for (sampleIdx[0] = 0; sampleIdx[0] < sampleSize[0]; sampleIdx[0] ++) {
	       for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1] ++) {
		    curLabel = samplePtr->GetPixel(sampleIdx);
		    sumexp = 0;
		    for (thisLabel = 0; thisLabel < numClusters; thisLabel++) {
			 // Compute S_i
			 S = 0; 
			 if (sampleIdx[0] > 0) {
			      nidx[0] = sampleIdx[0] - 1;
			      nidx[1] = sampleIdx[1];
			      S = S + int(thisLabel != samplePtr->GetPixel(nidx));
			 }
			 
			 if  (sampleIdx[0] < sampleSize[0] - 1) {
			      nidx[0] = sampleIdx[0] +1;
			      nidx[1] = sampleIdx[1];
			      S = S + int(thisLabel != samplePtr->GetPixel(nidx));
			 }
			 
			 if  (sampleIdx[1] > 0) {
			      nidx[0] = sampleIdx[0];
			      nidx[1] = sampleIdx[1] - 1;
			      S = S + int(thisLabel != samplePtr->GetPixel(nidx));
			 }

			 if  (sampleIdx[1] < sampleSize[1] - 1) {
			      nidx[0] = sampleIdx[0];
			      nidx[1] = sampleIdx[1] + 1;
			      S = S + int(thisLabel != samplePtr->GetPixel(nidx));
			 }
			 sumexp = sumexp + exp(-beta * S);
			 if (thisLabel == curLabel) {
			      curS = S;
			 }
		    }
			 
		    Q += (-beta * curS - log(sumexp))/double(sampleSize[2]);
	       }
	  }
     }

     
     return Q;
}

// Maximization prior log-likelihood over beta. return beta.
double golden_section_search(double a, 
			     double b,
			     ImageType3D::Pointer samplePtr,
			     int numClusters)
{
     if (a >= b) {
	  printf("golden_section_search: a must be smaller than b.\n");
	  exit(1);
     }
     double tau = (sqrt(5) - 1)/2;
     double x1 = 0, x2 = 0;
     double f1 = 0, f2 = 0;
     x1 = a + (1 - tau) * (b - a);
     f1 =  - eval_ll(samplePtr, numClusters, x1);
     x2 = a + tau * (b - a);

     // compute the minimal value of negagive log-likelihood.
     f2 =  - eval_ll(samplePtr, numClusters, x2);
     while (abs(b - a) > 0.0001) {
	  if (f1 > f2) {
	       a = x1;
	       x1 = x2;
	       f1 = f2;
	       x2 = a + tau * (b - a);
	       f2 =  - eval_ll(samplePtr, numClusters, x2);
	  }
	  else {
	       b = x2;
	       x2 = x1;
	       f2 = f1;
	       x1 = a + (1 - tau) * (b - a);
	       f1 =  - eval_ll(samplePtr, numClusters, x1);
	  }
	  printf("a = %f, x1 = %f, x2 = %f, b = %f, f1 = %f, f2 = %f\n", a, x1, x2, b, f1, f2);
     }
     return (b+a)/2;
}

int estep(ImageType2D::Pointer imagePtr,
	  ImageType3D::Pointer samplePtr,
	  CompType* cls,
	  int numClusters,
	  int burnin,
	  int numScan,
	  double beta)

{
     int scan = 0;
     double p_acpt = 0;
     double denergy = 0;
     int cand = 1;
     int label = 0;
     int nLabel = 0; // neighobr label
     float intensity = 0; // current pixel intensity.

     ImageType2D::IndexType idx;     
     ImageType2D::IndexType nidx;    

     ImageType3D::IndexType idx3;      
//     ImageType3D::IndexType sampleIdx;     
     ImageType3D::SizeType sampleSize = samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType2D::SizeType size =
          imagePtr->GetLargestPossibleRegion().GetSize();

     // Allocate memory for a single sample of MC.
     ImageType2D::Pointer labelPtr = ImageType2D::New();
     ImageType2D::IndexType start2d;
     start2d[0] = 0;
     start2d[1] = 0;
     ImageType2D::RegionType region2d;
     region2d.SetSize(size);
     region2d.SetIndex(start2d);
     labelPtr->SetRegions(region2d);
     labelPtr->Allocate();

     //////////////////////////////////////////////////////////////
     ///////        Define random number generator.        ////////
     //////////////////////////////////////////////////////////////

     mygenerator.seed(static_cast<unsigned int>(std::time(0)));

     // Uniform integer generator.
     boost::uniform_int<> uni_int(1, numClusters); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     // Init the label. We save the init'd labels into the last slices of samplePtr.
     idx3[2] = sampleSize[2] - 1;
     for (idx3[0] = 0; idx3[0] < sampleSize[0]; idx3[0] ++) {
     	  for (idx3[1] = 0; idx3[1] < sampleSize[1]; idx3[1]++) {
	       label = roll_die() - 1;	       
	       samplePtr->SetPixel(idx3, label);
     	  }
     }

     // copy the init'd labels to labelPtr.
     idx3[2] = sampleSize[2] - 1;
     for (idx3[0] = 0; idx3[0] < sampleSize[0]; idx3[0] ++) {
     	  for (idx3[1] = 0; idx3[1] < sampleSize[1]; idx3[1]++) {
	       idx[0] = idx3[0];
	       idx[1] = idx3[1];
	       labelPtr->SetPixel(idx, samplePtr->GetPixel(idx3));
     	  }
     }

     // sampling Markov Random Field. 
     scan = 0;
     while (scan < burnin + numScan) {
	  

     	  for (idx[0] = 0; idx[0] < size[0]; idx[0]++) {
     	       for (idx[1] = 0; idx[1] < size[1]; idx[1]++) {
		    cand = roll_die() - 1;
		    label = labelPtr->GetPixel(idx);
     		    denergy = 0;
     		    if (idx[0] > 0) {
			 nidx[0] = idx[0] - 1;
			 nidx[1] = idx[1];
			 nLabel = labelPtr->GetPixel(nidx); // neighbor label
			 denergy = denergy + int(cand != nLabel) - int(label != nLabel);
		    }

     		    if (idx[0] < size[0]-1) {
			 nidx[0] = idx[0] + 1;
			 nidx[1] = idx[1];
			 nLabel = labelPtr->GetPixel(nidx); // neighbor label
			 denergy = denergy + int(cand != nLabel) - int(label != nLabel);
		    }

     		    if (idx[1] > 0) {
			 nidx[0] = idx[0];
			 nidx[1] = idx[1] - 1;
			 nLabel = labelPtr->GetPixel(nidx); // neighbor label
			 denergy = denergy + int(cand != nLabel) - int(label != nLabel);
		    }

     		    if (idx[1] < size[1]-1) {
			 nidx[0] = idx[0];
			 nidx[1] = idx[1] + 1;
			 nLabel = labelPtr->GetPixel(nidx); // neighbor label
			 denergy = denergy + int(cand != nLabel) - int(label != nLabel);
		    }

     		    denergy = beta * denergy;

		    // likelihood energy.
		    intensity = imagePtr->GetPixel(idx);
		    denergy = denergy
			 + (intensity - cls[cand].mu)*(intensity - cls[cand].mu)/(2*cls[cand].sigma * cls[cand].sigma) + log(cls[cand].sigma) // Candidate's likelihood energy.
			 - (intensity - cls[label].mu)*(intensity - cls[label].mu)/(2*cls[label].sigma * cls[label].sigma) - log(cls[label].sigma); // Current label's likelihood energy.
		    
     		    // if energy change less than zero, just accept
     		    // candidate. otherwise accept with exp(- energy
     		    // change).

     		    if (denergy <= 0) {
			 
			 labelPtr->SetPixel(idx, cand);
     		    }
     		    else {
     			 p_acpt = exp(-denergy);
     			 if (uni() < p_acpt) {
			      labelPtr->SetPixel(idx, cand);
     			 }
     		    }
     	       }
     	  }
	  
	  if (scan >= burnin) {
	       // save current labels into samplePtr.
	       idx3[2] = scan - burnin;
	       for (idx3[0] = 0; idx3[0] < sampleSize[0]; idx3[0] ++) {
		    for (idx3[1] = 0; idx3[1] < sampleSize[1]; idx3[1]++) {
			 idx[0] = idx3[0];
			 idx[1] = idx3[1];
			 samplePtr->SetPixel(idx3, labelPtr->GetPixel(idx));
		    }
	       }
	  }

	  printf("Scan %i...", ++scan);
     }
     printf("\n");
     return 0;
}

// Given mu and sigma, evaluate the conditional log-likelihood P(d | f).
double eval_cll(ImageType3D::Pointer samplePtr,
		ImageType2D::Pointer imagePtr,
		CompType* cls)

{
     ImageType3D::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3D::IndexType sampleIdx;     
     ImageType2D::IndexType idx;     
     double cll = 0;
     int k = 0;
     
     for (sampleIdx[0] = 0;sampleIdx[0] < sampleSize[0]; sampleIdx[0] ++) {
	  for (sampleIdx[1] = 0; sampleIdx[1] < sampleSize[1]; sampleIdx[1] ++) {
	       for (sampleIdx[2] = 0; sampleIdx[2] < sampleSize[2]; sampleIdx[2] ++) {
		    idx[0] = sampleIdx[0];
		    idx[1] = sampleIdx[1];
		    k = samplePtr->GetPixel(sampleIdx);
		    cll = cll + (1/double(sampleSize[2])) * (
			 pow((imagePtr->GetPixel(idx) - cls[k].mu), 2)/(-2*pow(cls[k].sigma, 2))
			 - log(cls[k].sigma)
			 );
	       } 
	  } 
     }
     return cll;
}
