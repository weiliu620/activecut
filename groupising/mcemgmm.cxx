#include "common_mcem.h"
#include "gmmodel_mcem.h"

using std::cout;
using std::endl;

twister_base_gen_type mygenerator(42u);

unsigned ReadObsToVec(std::string obspath,
		      std::vector<ImageType3DFloat::Pointer> & obsVec,
		      ImageType3DChar::Pointer maskPtr);

int InitGrp(std::vector< ImageType3DChar::Pointer > & grpVec,
	    ImageType3DChar::Pointer initGrpPtr);

int SaveGrpSamples(std::vector<ImageType3DChar::Pointer>  &  grpVec,
		   ImageType3DChar::Pointer maskPtr,
		   std::string outFilename);


int main(int argc, char* argv[])
{
     ParStruct par;
     unsigned int EMIter = 1;
     unsigned seed = 0;
     float Q12 = 0; // -log P(G) - log(P(Z|G)
     float Q3 = 0; // -log P(X|Z)

     std::string obspath, initGrplabel, outGrpLabel, outGrpProb, outSubBase;
     std::string samplePrefix, grpSampleName, outSubPrefix;
     std::string opmethod;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Inference on group mrf.")
	  ("emiter", po::value<unsigned int>(&EMIter)->default_value(30),
	   "initial number of Monte Carlo samples. ")
	  ("numClusters,k", po::value<unsigned>(&par.numClusters)->default_value(6),
	       "Number of clusters. Default is 6.")
	  ("burnin,b", po::value<unsigned>(&par.burnin)->default_value(20),
	   "number of scans for burn-in period. ")
	  ("numSamples,n", po::value<unsigned >(&par.numSamples)->default_value(20),
	   "initial number of Monte Carlo samples. ")
	  ("inittemp", po::value<float>(&par.initTemp)->default_value(1),
	   "Initial temperature for annealing.")
	  ("finaltemp", po::value<float>(&par.finalTemp)->default_value(1),
	   "Final temperature for annealing to stop")
	  ("verbose,v", po::value<unsigned short>(&par.verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")
	  ("betag", po::value<float>(&par.betag)->default_value(0.5),
	   "group level pairwise interation term")

	  ("betaz", po::value<float>(&par.betaz)->default_value(0.5),
	   "subject level pairwise interation term")
	  ("alpha", po::value<float>(&par.alpha)->default_value(0.5),
	   "connection between group lebel and individual subject level.")

	  ("gamma", po::value<float>(&par.gamma)->default_value(1),
	   " weights of lilelihood term.")

	  ("opmethod", po::value<std::string>(&opmethod)->default_value("mc"),
	   "optimization method. 'mc' or 'gc'. ")

	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")

	  ("obspath,f", po::value<std::string>(&obspath)->default_value("."),
	   "noised image file")

	  ("sampleprefix", po::value<std::string>(&samplePrefix)->default_value("subject"), "Monte Carlo samples file name prefix")

	  ("grpsamplename", po::value<std::string>(&grpSampleName)->default_value("grpSample.nii.gz"), "Monte Carlo group label samples file name.")

	  ("grouplabel,g", po::value<std::string>(&outGrpLabel)->default_value("outgrplabel.nii"), "output group label file. nii format.")

	  ("outsubprefix", po::value<std::string>(&outSubPrefix)->default_value("outSub"), "output subject label map file prefix.")

	   ("initgrouplabel,i", po::value<std::string>(&initGrplabel)->default_value("grouplabel.nii"), "Initial group level label map. Also used as mask file.");

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

     // define empty group sample vector.
     std::vector<ImageType3DChar::Pointer> grpVec(par.numSamples);
     
     ReaderType3DChar::Pointer initGrpReader = ReaderType3DChar::New();
     initGrpReader->SetFileName(initGrplabel);

     // original 1-based group map only used for mask.
     ImageType3DChar::Pointer maskPtr = initGrpReader->GetOutput();
     maskPtr->Update();

     AddConstantToImageFilterType::Pointer myfilter = AddConstantToImageFilterType::New();

     // keep original input image unchanged.
     myfilter->InPlaceOff();
     myfilter->SetInput( initGrpReader->GetOutput());
     myfilter->SetConstant(-1);
     myfilter->GetOutput()->Update();

     InitGrp(grpVec, myfilter->GetOutput() );

     // observed images.
     std::vector<ImageType3DFloat::Pointer>  obsVec;
     par.numSubs = ReadObsToVec(obspath, obsVec, maskPtr);
     VnlVectorImageType::SizeType obsSize = obsVec[0]->GetLargestPossibleRegion().GetSize();

     // Allocate memory for saving samples of MC.
     std::vector< std::vector<ImageType3DChar::Pointer> >  
	  sampleVec(par.numSubs, std::vector<ImageType3DChar::Pointer>(par.numSamples) );
     // Allocate memory for sampleVec, and apply mask to sampleVec and init
     // samples same with initial group label.
     InitSamples(sampleVec, myfilter->GetOutput());

     GMModel gmmodel(obsVec,
		     grpVec,
		     sampleVec,
		     maskPtr,
		     par);

     if(opmethod.compare("gc") == 0) {
	  gmmodel.InitGraph(obsVec, grpVec, sampleVec, maskPtr);
     }
     save3dchar(grpVec[par.numSamples-1], outGrpLabel);
     if (par.verbose >= 1) {
     }

     // Since group and sub all init with same label map, just estimate
     // parameters mu and sigma with same function in M step.
     gmmodel.estimatePar(sampleVec, obsVec, maskPtr);
     // gmmodel.EstimatePriorPar(sampleVec, grpVec, maskPtr, "init");

     gmmodel.printSelf("colname");
     gmmodel.printSelf("table");

     for (unsigned short emIterIdx = 0; emIterIdx < EMIter; emIterIdx ++) {

	  if (par.verbose >= 1) {
	       printf("EM iteration %i begin:\n", emIterIdx + 1);
	  }

	  gmmodel.SetTemperature(par.initTemp * pow( (par.finalTemp / par.initTemp), float(emIterIdx) / float(EMIter) ) );

	  if (opmethod.compare("mc") == 0) {
	       // sampling both group and individual subject label map.
	       gmmodel.Sampling(grpVec, obsVec, maskPtr, sampleVec, par.burnin);
	  }
	  else if(opmethod.compare("gc") == 0) {
	       gmmodel.graphcuts(obsVec, grpVec, sampleVec, maskPtr);
	  }
	  else {
	       fprintf(stderr, "opmethod have to be either 'mc' or 'gc'.\n");
	  }


	  if (par.verbose >=1 ) {
	       if (opmethod.compare("mc") == 0 ) {
		    SaveSamples(sampleVec, maskPtr, samplePrefix);
		    SaveGrpSamples(grpVec, maskPtr, grpSampleName);

	       }
	       else if (opmethod.compare("gc") == 0) {
		    
	       }
	  }

	  if (par.verbose >= 1) {
	       printf("EM iteration %i, sampling done. \n", emIterIdx + 1);
	  }

	  // estimate parameters mu, kappa, and prior parameters (betag, betaz,
	  // alpha,e tc).
	  if (opmethod.compare("mc") == 0) {
	       gmmodel.estimatePar(sampleVec, obsVec, maskPtr);
	       Q12 = gmmodel.EstimatePriorPar(sampleVec, grpVec, maskPtr, "normal");
	  }
	  else if(opmethod.compare("gc") == 0) {
	       
	  }

	  printf("EM iteration %i, parameter estimation done. \n", emIterIdx + 1 );
	  gmmodel.printSelf("table");
	  save3dchar(grpVec[par.numSamples-1], outGrpLabel);


	  // compute log likelihood of Q2 -- log P(Z | G)
	  if (par.verbose >= 1) {
	       float Q3 = gmmodel.Q3(sampleVec, obsVec, maskPtr);
	       printf("After parameter estimation, (Q1+Q2) = %.3f, Q3 = %.3f, (Q1+Q2+Q3) = %.3f\n", Q12, Q3, (Q12+Q3));
	  }


	  // estimate labels by ICM. This is used as final label output, and
	  // also as initial label map for sampling.
	  gmmodel.ICM(grpVec, obsVec, maskPtr, sampleVec, 0.01); 
	  if (opmethod.compare("mc") == 0 ) {
	       SaveSamples(sampleVec, maskPtr, samplePrefix);
	       SaveGrpSamples(grpVec, maskPtr, grpSampleName);
	       save3dchar(grpVec[par.numSamples-1], outGrpLabel);
	  }
	  else if (opmethod.compare("gc") == 0) {
	       // TBA.
	  }


     }
     
     // Save final resutls for both group label map and subjects label
     // map. 
     save3dchar(grpVec[par.numSamples-1], outGrpLabel);
     for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  std::string subFileNum = boost::lexical_cast<std::string> (subIdx+1);
	  std::string thisSubFilename(outSubPrefix);
	  // add sub number.
	  thisSubFilename.append(subFileNum);
	  thisSubFilename.append(".nii.gz");
	  save3dchar(sampleVec[subIdx][par.numSamples-1], thisSubFilename);
     }
     
}


// Init all samples with grpPtr, i.. initGrp
int InitSamples(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec, 
		ImageType3DChar::Pointer grpPtr)
{
     // group.
     ImageType3DChar::SizeType grpSize = grpPtr->GetLargestPossibleRegion().GetSize();
     ConstIteratorType3DChar grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );

     // samples.
     ImageType3DChar::IndexType start;
     start.Fill(0);
     ImageType3DChar::SizeType sampleSize;
     sampleSize = grpSize;

     ImageType3DChar::RegionType sampleRegion;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(start);

     for (unsigned subIdx = 0; subIdx < sampleVec.size(); subIdx ++) {
	  for (unsigned mcIdx = 0; mcIdx < sampleVec[0].size(); mcIdx ++) {
	       sampleVec[subIdx][mcIdx] = ImageType3DChar::New();
	       sampleVec[subIdx][mcIdx]->SetRegions(sampleRegion);
	       sampleVec[subIdx][mcIdx]->Allocate();
	       sampleVec[subIdx][mcIdx]->FillBuffer(0);

	       sampleVec[subIdx][mcIdx]->SetOrigin(grpPtr->GetOrigin());
	       sampleVec[subIdx][mcIdx]->SetSpacing(grpPtr->GetSpacing());
	       sampleVec[subIdx][mcIdx]->SetDirection(grpPtr->GetDirection());

	       
	       IteratorType3DChar sampleIt(sampleVec[subIdx][mcIdx],
					   sampleVec[subIdx][mcIdx]->GetLargestPossibleRegion() );
	       
	       for (grpIt.GoToBegin(), sampleIt.GoToBegin(); !grpIt.IsAtEnd(); ++ grpIt, ++ sampleIt) {
		    sampleIt.Set( grpIt.Get() );
	       }
	  } // mcIdx
     } // subIdx.
     return 0;
}


// Init all images in grpVec with initGrpPtr. 
int InitGrp(std::vector< ImageType3DChar::Pointer > & grpVec,
	    ImageType3DChar::Pointer initGrpPtr)
{
     // source.
     ConstIteratorType3DChar initGrpIt(initGrpPtr, initGrpPtr->GetLargestPossibleRegion() );

     // destination.
     ImageType3DChar::IndexType start;
     start.Fill(0);
     ImageType3DChar::SizeType grpSize = initGrpPtr->GetLargestPossibleRegion().GetSize();

     ImageType3DChar::RegionType grpRegion;
     grpRegion.SetSize(grpSize);
     grpRegion.SetIndex(start);

     for (unsigned mcIdx = 0; mcIdx < grpVec.size(); mcIdx ++) {
	  grpVec[mcIdx] = ImageType3DChar::New();
	  grpVec[mcIdx]->SetRegions(grpRegion);
	  grpVec[mcIdx]->Allocate();
	  grpVec[mcIdx]->FillBuffer(-1);
	       
	  IteratorType3DChar grpIt(grpVec[mcIdx],
				   grpVec[mcIdx]->GetLargestPossibleRegion() );
	       
	  for (initGrpIt.GoToBegin(), grpIt.GoToBegin(); !grpIt.IsAtEnd(); ++ grpIt, ++ initGrpIt) {
	       grpIt.Set( initGrpIt.Get() );
	  }
     } // mcIdx
     return 0;
}


unsigned ReadObsToVec(std::string obspath,
		      std::vector<ImageType3DFloat::Pointer> & obsVec,
		      ImageType3DChar::Pointer maskPtr)
{
     // some of the code below refers to boost filesystem tutorial tut04.cpp.
     
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     // default constructed iterator acts as a end iterator.
     boost::filesystem::path obspathVar(obspath);
     boost::filesystem::directory_iterator obspathEnd;

     ReaderType3DFloat::Pointer obsReader = ReaderType3DFloat::New();

     
     // prepare for sorting directory entries.
     typedef std::vector<boost::filesystem::path> PathVec;
     PathVec  sortedEntries;
     copy(boost::filesystem::directory_iterator(obspath), boost::filesystem::directory_iterator(), std::back_inserter(sortedEntries) );
     sort(sortedEntries.begin(), sortedEntries.end() );

     // Get number of subjects.
     unsigned numSubs = 0;
     // for (boost::filesystem::directory_iterator obspathIt(obspathVar); obspathIt != obspathEnd; ++obspathIt) {
     // 	  numSubs ++;
     // }

     numSubs = sortedEntries.size();
     obsVec.resize(numSubs);

     ImageType3DFloat::Pointer singleObsPtr = ImageType3DFloat::New();

     // allocate memory for obsPtr.
     ImageType3DFloat::IndexType obsIdx;
     obsIdx.Fill ( 0 );
     ImageType3DFloat::SizeType obsSize;

     // get image size.
     boost::filesystem::directory_iterator obspathIt(obspathVar);
     obsReader->SetFileName( (*obspathIt).path().string() );
     obsReader->Update();
     singleObsPtr = obsReader->GetOutput();
     obsSize = singleObsPtr->GetLargestPossibleRegion().GetSize();

     // region
     ImageType3DFloat::RegionType obsRegion;
     obsRegion.SetSize (obsSize);
     obsRegion.SetIndex( obsIdx );



     unsigned subIdx = 0;

     // for (boost::filesystem::directory_iterator obspathIt(obspathVar); obspathIt != obspathEnd; ++obspathIt) {
     for (PathVec::const_iterator obspathIt(sortedEntries.begin() ); obspathIt != sortedEntries.end(); ++ obspathIt) {


	  // Destination.
	  obsVec[subIdx] = ImageType3DFloat::New();
	  obsVec[subIdx]->SetRegions( obsRegion );
	  obsVec[subIdx]->Allocate();
	  obsVec[subIdx]->FillBuffer ( 0 );
	  IteratorType3DFloat obsIt(obsVec[subIdx], obsVec[subIdx]->GetLargestPossibleRegion() );

	  // source.
	  // obsReader->SetFileName( (*obspathIt).path().string() );
	  obsReader->SetFileName( (*obspathIt).string() );
	  obsReader->Update();
	  singleObsPtr = obsReader->GetOutput();
	  IteratorType3DFloat singleObsIt(singleObsPtr, singleObsPtr->GetLargestPossibleRegion() );
	  std::cout <<  "add " << (*obspathIt).string() << "\n";

	  // read in data for this subject.

	  for (maskIt.GoToBegin(), obsIt.GoToBegin(), singleObsIt.GoToBegin(); !obsIt.IsAtEnd(); ++ obsIt, ++maskIt, ++singleObsIt){
	       if (maskIt.Get() > 0) {
		    obsIt.Set( singleObsIt.Get() );
	       } // maskIt > 0
	  } // maskIt
	  
	  // done with read data for this sub. 
	  subIdx ++;
     }
     return numSubs;
}


