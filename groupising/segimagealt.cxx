#include "commonalt.h"
#include "MCModel.h"
int SplitMRI(ImageType5DFloat::Pointer inPtr,      
	     std::vector<VnlVectorImageType::Pointer> & fmriVec);


using std::cout;
using std::endl;

twister_base_gen_type mygenerator(42u);


int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     float beta = 0.1;
     int burnin = 1;
     unsigned seed = 0;

     unsigned int numSamples = 1;
     unsigned int EMIter = 1;

     int numClusters = 4;
     float initTemp = 0.5;
     float finalTemp = 0.01;

     float alpha = 0.5;
     float beta_g = 0.5;
     float beta_z = 0.5;

     double meanKappa = 0;
     double stdKappa = 0;

     std::string fmriPath, initGrplabel, outGrpLabel, outGrpProb, outSubBase;
     std::string samplePrefix;
     std::string samplingMethod;
     std::string estMethod;
     std::string cheat;
     std::string parfile;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "Inference on group mrf.")
	  ("burnin,b", po::value<int>(&burnin)->default_value(20),
	   "number of scans for burn-in period. ")
	  ("inittemp", po::value<float>(&initTemp)->default_value(0.5),
	   "Initial temperature for annealing.")
	  ("finaltemp", po::value<float>(&finalTemp)->default_value(0.01),
	   "Final temperature for annealing to stop")
	  ("numSamples,n", po::value<unsigned int>(&numSamples)->default_value(20),
	   "initial number of Monte Carlo samples. ")
	  ("emiter", po::value<unsigned int>(&EMIter)->default_value(30),
	   "initial number of Monte Carlo samples. ")

	  ("betag", po::value<float>(&beta_g)->default_value(0.5),
	   "group level pairwise interation term")

	  ("betaz", po::value<float>(&beta_z)->default_value(0.5),
	   "subject level pairwise interation term")
	  ("alpha", po::value<float>(&alpha)->default_value(0.5),
	   "connection between group lebel and individual subject level.")

	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("numClusters,k", po::value<int>(&numClusters)->default_value(6),
	       "Number of clusters. Default is 6.")
	  ("meankappa", po::value<double>(&meanKappa)->default_value(100),
	   "Hyper-parameter for kappa: mean")
	  ("stdkappa", po::value<double>(&stdKappa)->default_value(20),
	   "Hyper-parameter for kappa: standard deviation")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")

	  ("fmripath,f", po::value<std::string>(&fmriPath)->default_value("."),
	   "noised image file")

	  ("grouplabel,g", po::value<std::string>(&outGrpLabel)->default_value("outgrplabel.nii"), "output group label file. nii format.")

	  ("sampleprefix", po::value<std::string>(&samplePrefix)->default_value("subject"), "Monte Carlo samples file name prefix")

	   ("subbasename", po::value<std::string>(&outSubBase)->default_value("subject"),
	    "Output individual subject's base name (with path).")

	   ("initgrouplabel,i", po::value<std::string>(&initGrplabel)->default_value("grouplabel.nii"), "Initial group level label map. Also used as mask file.");

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
     }

     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // init random generator.
     mygenerator.seed(static_cast<unsigned int>(seed));

	   
     boost::filesystem::path fmriPathVar(fmriPath);
     boost::filesystem::directory_iterator fmriPathEnd;


     SeriesReaderType5DFloat::Pointer fmriReader = SeriesReaderType5DFloat::New();

     for (boost::filesystem::directory_iterator fmriPathIt(fmriPathVar); fmriPathIt != fmriPathEnd; ++fmriPathIt) {

     	  fmriReader->AddFileName( (*fmriPathIt).path().string());	  
     	  cout <<  "add " << (*fmriPathIt).path().string() << "\n";
     }

     fmriReader->Update();
     ImageType5DFloat::SizeType fmriSize = fmriReader->GetOutput()->GetLargestPossibleRegion().GetSize();
     unsigned numSubs = fmriSize[4];

     // Allocate memory for saving vectorized fmri image.
     std::vector<VnlVectorImageType::Pointer>  fmriVec(numSubs);
     SplitMRI(fmriReader->GetOutput(), fmriVec);

     // read in initial grp label. Also used as mask. 1-based label converted to
     // 0-based label.
     ReaderType3DChar::Pointer grpReader = ReaderType3DChar::New();
     grpReader->SetFileName(initGrplabel);
     AddConstantToImageFilterType::Pointer myfilter = AddConstantToImageFilterType::New();
     myfilter->InPlaceOff();
     myfilter->SetInput( grpReader->GetOutput());
     myfilter->SetConstant(-1);
     ImageType3DChar::Pointer grpPtr = myfilter->GetOutput();
     grpPtr->Update();

     // original 1-based group map only used for mask.
     ImageType3DChar::Pointer maskPtr = grpReader->GetOutput();


     // Allocate memory for saving samples of MC.
     std::vector< std::vector<ImageType3DChar::Pointer> >  
	  sampleVec(numSubs, std::vector<ImageType3DChar::Pointer>(numSamples) );
     
     // Allocate memory for sampleVec, and apply mask to sampleVec and init
     // samples same with initial group label.
     InitSamples(sampleVec, grpPtr);

     MCModel mcmodel(fmriVec,
     		     sampleVec,
		     grpPtr,
		     maskPtr,
		     numClusters,
		     alpha,
		     beta_g,
		     beta_z,
     		     verbose);

     if (verbose >= 1) {
	  SaveSamples(sampleVec, grpPtr, samplePrefix);
	  save3dchar(grpPtr, outGrpLabel);
     }

     mcmodel.estimateMu(grpPtr, fmriVec, maskPtr);
     mcmodel.estimateKappaWithPrior(meanKappa, stdKappa, fmriSize[3]);
     mcmodel.printSelf("normal");
     
     for (unsigned short emIterIdx = 0; emIterIdx < EMIter; emIterIdx ++) {

	  if (verbose >= 1) {
	       printf("EM iteration %i begin:\n", emIterIdx);
	  }

	  mcmodel.SetTemperature(initTemp * pow( (finalTemp / initTemp), float(emIterIdx) / float(EMIter) ) );

	  // sampling both group and individual subject label map.
	  mcmodel.mcsampling(grpPtr, fmriVec, maskPtr, sampleVec, burnin);

	  if (verbose >= 1) {
	       printf("EM iteration %i, sampling done. \n", emIterIdx);
	       SaveSamples(sampleVec, maskPtr, samplePrefix);
	       save3dchar(grpPtr, outGrpLabel);
	  }

	  // estimate vMF parameters mu, kappa.
	  mcmodel.EstimateMuFromSample(sampleVec, fmriVec, maskPtr);
	  mcmodel.estimateKappaWithPrior(meanKappa, stdKappa, fmriSize[3]);

	  printf("EM iteration %i, parameter estimation done. \n", emIterIdx);
	  mcmodel.printSelf("normal");
	  fflush(stdout);
	  
	  mcmodel.llhood(grpPtr, sampleVec, fmriVec, maskPtr);

	  if (verbose >= 1) {
	       for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
		    std::string subFileNum = boost::lexical_cast<std::string> (subIdx+1);
		    std::string thisSubFilename(outSubBase);
		    // add sub number.
		    thisSubFilename.append(subFileNum);
		    thisSubFilename.append(".nii.gz");
		    save3dchar(sampleVec[subIdx][numSamples - 1], thisSubFilename);
	       }

	       // Group label map.
	       printf("begin print group label map.\n");
	       save3dchar(grpPtr, outGrpLabel);
	       printf("end print group label map.\n");
	  }
     }
     
     // EM done. Save group map and sub map.
     
     // Individual subjects label map.
     for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	  std::string subFileNum = boost::lexical_cast<std::string> (subIdx+1);
	  std::string thisSubFilename(outSubBase);
	  // add sub number.
	  thisSubFilename.append(subFileNum);
	  thisSubFilename.append(".nii.gz");
	  save3dchar(sampleVec[subIdx][numSamples - 1], thisSubFilename);
     }

     // Group label map.
     save3dchar(grpPtr, outGrpLabel);
}

int SplitMRI(ImageType5DFloat::Pointer inPtr,      
	     std::vector<VnlVectorImageType::Pointer> & fmriVec)
{
     VnlVectorImageType::IndexType fmriIdx;
     fmriIdx.Fill(0);
     VnlVectorImageType::SizeType fmriSize;


     ImageType5DFloat::IndexType inIdx;
     inIdx.Fill(0);
     ImageType5DFloat::SizeType inSize = inPtr->GetLargestPossibleRegion().GetSize();
     fmriSize[0] = inSize[0];
     fmriSize[1] = inSize[1];
     fmriSize[2] = inSize[2];
     VnlVectorImageType::RegionType fmriRegion;
     fmriRegion.SetSize(fmriSize);
     fmriRegion.SetIndex(fmriIdx);

     unsigned numSubs = inSize[4];
     vnl_vector<float> timeSeries( inSize[3], 0 );
     for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	  // first allocate memory for fmriPtr.
	  fmriVec[subIdx] = VnlVectorImageType::New();
	  fmriVec[subIdx]->SetRegions( fmriRegion );
	  fmriVec[subIdx]->Allocate();
	  fmriVec[subIdx]->FillBuffer( timeSeries );

	  IteratorTypeVnlVector fmriIt(fmriVec[subIdx], fmriVec[subIdx]->GetLargestPossibleRegion() );

	  // copy time series from inPtr to fmriPtr.
	  inIdx[4] = subIdx;
	  for (fmriIt.GoToBegin(); !fmriIt.IsAtEnd(); ++ fmriIt) {
	       fmriIdx = fmriIt.GetIndex();
	       inIdx[0] = fmriIdx[0];
	       inIdx[1] = fmriIdx[1];
	       inIdx[2] = fmriIdx[2];
	       
	       for (inIdx[3] = 0; inIdx[3] < inSize[3]; inIdx[3] ++) {
		    timeSeries[inIdx[3]] = inPtr->GetPixel( inIdx );
	       }

	       fmriIt.Set ( timeSeries );
	  } // fmriIt
     } // subIdx

     return 0;
}
