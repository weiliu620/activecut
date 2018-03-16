#include "commonalt.h"
#include "MCModel.h"
#include <utilalt.h>

unsigned ReadObsToVec(std::string fmripath,
		      std::vector<VnlVectorImageType::Pointer> & fmriVec,
		      ImageType3DChar::Pointer maskPtr);

int InitGrp(std::vector< ImageType3DChar::Pointer > & grpVec,
	    ImageType3DChar::Pointer initGrpPtr);

int LLEnergy(ParStruct & par, MCModel & mcmodel,
	     std::vector<VnlVectorImageType::Pointer> & fmriVec,
	     std::vector< std::vector<ImageType3DChar::Pointer> >  & sampleVec,
	     ImageType3DChar::Pointer maskPtr);

int InitSubSamples(std::string srcpath,
		   std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec, 
		   ImageType3DChar::Pointer maskPtr);

int InitSamples(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec, 
		ImageType3DChar::Pointer grpPtr,
     		ImageType3DChar::Pointer maskPtr);

twister_base_gen_type mygenerator(42u);

int main(int argc, char* argv[])
{
     ParStruct par;
     unsigned seed = 0;
     unsigned int EMIter = 1;

     std::string fmriPath, initGrplabel, outGrpLabel, outGrpProb, outSubBase, 
	  grpSampleName, samplePrefix, initsubpath, estprior, initsame, groupprob, subprobbase;

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Inference on group mrf.")
	  ("burnin,b", po::value<unsigned>(&par.burnin)->default_value(20),
	   "number of scans for burn-in period. ")
	  ("inittemp", po::value<float>(&par.initTemp)->default_value(1),
	   "Initial temperature for annealing.")
	  ("finaltemp", po::value<float>(&par.finalTemp)->default_value(1),
	   "Final temperature for annealing to stop")
	  ("numSamples,n", po::value<unsigned int>(&par.numSamples)->default_value(20),
	   "initial number of Monte Carlo samples. ")
	  ("emiter", po::value<unsigned int>(&EMIter)->default_value(30),
	   "initial number of Monte Carlo samples. ")

	  ("betag", po::value<float>(&par.betag)->default_value(0.5),
	   "group level pairwise interation term")
	  ("betaz", po::value<float>(&par.betaz)->default_value(0.8),
	   "subject level pairwise interation term")
	  ("alpha", po::value<float>(&par.alpha)->default_value(0.7),
	   "connection between group lebel and individual subject level.")
	  ("gamma", po::value<float>(&par.gamma)->default_value(1),
	   " weights of lilelihood term.")

	  ("estprior", po::value<std::string>(&estprior)->default_value("yes"), "whether to estimate alpha and beta, the prior parameter.")
	  ("initsame", po::value<std::string>(&initsame)->default_value("yes"), "whether to init subject same to group label. choose yes to init same, otherwise init with subject label map.")

	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("numClusters,k", po::value<unsigned>(&par.numClusters)->default_value(6),
	   "Number of clusters. Default is 6.")

	   ("initgrouplabel,i", po::value<std::string>(&initGrplabel)->default_value("grouplabel.nii"), "Initial group level label map. Also used as mask file.")

	   ("initsubpath,t", po::value<std::string>(&initsubpath)->default_value("subpath"), "Initial subject label maps path")

	  ("fmripath,f", po::value<std::string>(&fmriPath)->default_value("."), "noised image file")

	  ("sampleprefix", po::value<std::string>(&samplePrefix)->default_value("subject"), "Monte Carlo samples file name prefix")

	   ("subbasename", po::value<std::string>(&outSubBase)->default_value("subject"), "Output individual subject's base name (with path).")

	   ("subprobbase", po::value<std::string>(&subprobbase)->default_value("subject"), "Output individual subject's probability map name (with path).")

	  ("grouplabel,g", po::value<std::string>(&outGrpLabel)->default_value("outgrplabel.nii"), "output group label file. nii format.")

	  ("groupprob,g", po::value<std::string>(&groupprob)->default_value("outgrplabel.nii"), "output group probability map.")

	  ("grpsamplename", po::value<std::string>(&grpSampleName)->default_value("grpSample.nii.gz"), "Monte Carlo group label samples file name.")


	  ("verbose,v", po::value<unsigned short>(&par.verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: segimagealt [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // init random generator.
     mygenerator.seed(static_cast<unsigned int>(seed));

     // read in initial grp label. Also used as mask. 1-based label converted to
     // 0-based label.
     ReaderType3DChar::Pointer grpReader = ReaderType3DChar::New();
     grpReader->SetFileName(initGrplabel);
     AddConstantToImageFilterType::Pointer myfilter = AddConstantToImageFilterType::New();
     myfilter->InPlaceOff();
     myfilter->SetInput( grpReader->GetOutput());
     myfilter->SetConstant(-1);
     myfilter->Update();

     // define empty group sample vector.
     std::vector<ImageType3DChar::Pointer> grpVec(par.numSamples);
     InitGrp(grpVec, myfilter->GetOutput() );

     // original 1-based group map only used for mask.
     ImageType3DChar::Pointer maskPtr = grpReader->GetOutput();

     // Allocate memory for saving vectorized fmri image.
     std::vector<VnlVectorImageType::Pointer>  fmriVec;


     // read observed fmri data into vector image.
     par.numSubs = ReadObsToVec(fmriPath, fmriVec, maskPtr);
     VnlVectorImageType::SizeType fmriSize = fmriVec[0]->GetLargestPossibleRegion().GetSize();

     // Allocate memory for saving subjects MC samples.
     std::vector< std::vector<ImageType3DChar::Pointer> >  
	  sampleVec(par.numSubs, std::vector<ImageType3DChar::Pointer>(par.numSamples) );


     
     // Init all subject samples.
     if (initsame.compare("yes") == 0) {
	  InitSamples(sampleVec, myfilter->GetOutput(), maskPtr);
     }
     else if (initsame.compare("no") == 0) {
	       InitSubSamples(initsubpath, sampleVec, maskPtr);
	  }
     else {
	  fprintf(stderr, "main(): initsame either set to yes or no.\n");
	  exit(1);
     }

     if (par.verbose >= 1) {
     	  SaveSamples(sampleVec, maskPtr, samplePrefix);
	  SaveGrpSamples(grpVec, maskPtr, grpSampleName);
     }

     MCModel mcmodel(fmriVec, sampleVec, grpVec, maskPtr, par);

     mcmodel.EstimateMuFromSample(sampleVec, fmriVec, maskPtr);
     mcmodel.estimateKappa();
     if (par.verbose >= 1)  mcmodel.printSelf("normal");

     mcmodel.ComputeEnergy(sampleVec, grpVec, maskPtr, fmriVec);

     for (unsigned short emIterIdx = 0; emIterIdx < EMIter; emIterIdx ++) {

     	  printf("EM iteration %i begin:\n", emIterIdx + 1);

	  mcmodel.SetTemperature(par.initTemp * pow( (par.finalTemp / par.initTemp), float(emIterIdx) / float(EMIter) ) );

     	  // sampling both group and individual subject label map.
     	  mcmodel.mcsampling(grpVec, fmriVec, maskPtr, sampleVec, par.burnin);

     	  if (par.verbose >= 0) {
     	       printf("EM iteration %i, sampling done. \n", emIterIdx + 1);
	       mcmodel.ComputeEnergy(sampleVec, grpVec, maskPtr, fmriVec);
     	  }

     	  if (par.verbose >= 1) {
     	       SaveSamples(sampleVec, maskPtr, samplePrefix);
	       SaveGrpSamples(grpVec, maskPtr, grpSampleName);
	  }

     	  // estimate vMF parameters mu, kappa.
     	  printf("EM iteration %i, parameter estimation begin. \n", emIterIdx + 1);
     	  mcmodel.EstimateMuFromSample(sampleVec, fmriVec, maskPtr);
     	  mcmodel.estimateKappa();

	  // estimate alpha, beta.
	  if (estprior.compare("yes") == 0) {
	       mcmodel.EstimatePriorPar(sampleVec, grpVec, maskPtr);
	  }

     	  printf("EM iteration %i, parameter estimation done. \n", emIterIdx + 1);
	  mcmodel.ComputeEnergy(sampleVec, grpVec, maskPtr, fmriVec);
     	  if (par.verbose >= 1)  mcmodel.printSelf("normal");
     }

     
     // EM done. Save group map and sub map.
     // Individual subjects label map.
     for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
     	  std::string subFileNum = boost::lexical_cast<std::string> (subIdx+1);
	  if (subIdx < 9) {
	       subFileNum.insert(0,"0");
	  }
     	  std::string thisSubFilename(outSubBase);
     	  // add sub number.
     	  thisSubFilename.append(subFileNum);

     	  thisSubFilename.append(".nii.gz");
     	  save3dcharInc(sampleVec[subIdx][par.numSamples - 1], thisSubFilename);
     }

     // Group label map.
     save3dcharInc(grpVec[par.numSamples-1], outGrpLabel);
     return 0;
}

int LLEnergy(ParStruct & par, MCModel & mcmodel,
	     std::vector<VnlVectorImageType::Pointer> & fmriVec,
	     std::vector< std::vector<ImageType3DChar::Pointer> >  & sampleVec,
	     ImageType3DChar::Pointer maskPtr)
{
     vnl_matrix<float> Q3(par.numSamples, par.numSubs, 0);
     
     for (unsigned mcIdx = 0; mcIdx < par.numSamples; mcIdx ++) {
	  for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	       Q3[mcIdx][subIdx] = mcmodel.ComputeQ3(fmriVec[subIdx], sampleVec[subIdx][mcIdx], maskPtr, subIdx);
	  }
     }
     printf("E(X|Z):\n");
     printVnlMatrixSci(Q3, 100);
     return 0;
}

unsigned ReadObsToVec(std::string srcpath,
		      std::vector<VnlVectorImageType::Pointer> & fmriVec,
		      ImageType3DChar::Pointer maskPtr)
{
     // some of the code below refers to boost filesystem tutorial tut04.cpp.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     // prepare for sorting directory entries.
     typedef std::vector<boost::filesystem::path> PathVec;
     PathVec  sortedEntries;
     copy(boost::filesystem::directory_iterator(srcpath), boost::filesystem::directory_iterator(), std::back_inserter(sortedEntries) );
     sort(sortedEntries.begin(), sortedEntries.end() );

     // Get number of subjects.
     unsigned numSubs  = sortedEntries.size();
     fmriVec.resize(numSubs);

     // get time series length.
     boost::filesystem::directory_iterator obspathIt(srcpath);
     ReaderType4DFloat::Pointer srcReader = ReaderType4DFloat::New();
     srcReader->SetFileName( (*obspathIt).path().string() );
     srcReader->Update();
     ImageType4DF::Pointer srcPtr = srcReader->GetOutput();
     ImageType4DF::SizeType srcSize = srcReader->GetOutput()->GetLargestPossibleRegion().GetSize();
     ImageType4DF::IndexType srcIdx;

     VnlVectorImageType::SizeType destSize = maskPtr->GetLargestPossibleRegion().GetSize();
     VnlVectorImageType::IndexType destIdx;
     destIdx.Fill(0);
     VnlVectorImageType::RegionType destRegion;
     destRegion.SetSize(destSize);
     destRegion.SetIndex(destIdx);


     unsigned subIdx = 0;
     vnl_vector<float> timeSeries( srcSize[3], 0 );
     vnl_vector<float> zeroVec( srcSize[3], 0 );
     for (PathVec::const_iterator srcpathIt(sortedEntries.begin() ); srcpathIt != sortedEntries.end(); ++ srcpathIt) {

	  // Source.
	  srcReader->SetFileName( (*srcpathIt).string() );
	  srcReader->Update();
	  srcPtr = srcReader->GetOutput();

	  // Destination.
	  fmriVec[subIdx] = VnlVectorImageType::New();
	  fmriVec[subIdx]->SetRegions( destRegion );
	  fmriVec[subIdx]->Allocate();
	  fmriVec[subIdx]->FillBuffer ( zeroVec );

	  // transfer data.
	  IteratorTypeVnlVector destIt(fmriVec[subIdx], fmriVec[subIdx]->GetLargestPossibleRegion() );

	  for (maskIt.GoToBegin(), destIt.GoToBegin(); !destIt.IsAtEnd(); ++ maskIt, ++ destIt) {
	       if (maskIt.Get() > 0) {
		    destIdx = destIt.GetIndex();
		    srcIdx[0] = destIdx[0];
		    srcIdx[1] = destIdx[1];
		    srcIdx[2] = destIdx[2];
	       
		    for (srcIdx[3] = 0; srcIdx[3] < srcSize[3]; srcIdx[3] ++) {
			 timeSeries[srcIdx[3]] = srcPtr->GetPixel( srcIdx );
		    }

		    destIt.Set ( timeSeries );
	       } // mask > 0
	  } // maskIt
	  subIdx ++;
     }
     return numSubs;
}


/* Read initial subject label map from files. */
int InitSubSamples(std::string srcpath,
		   std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec, 
		   ImageType3DChar::Pointer maskPtr)
{
     // some of the code below refers to boost filesystem tutorial tut04.cpp.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     // prepare for sorting directory entries.
     typedef std::vector<boost::filesystem::path> PathVec;
     PathVec  sortedEntries;
     copy(boost::filesystem::directory_iterator(srcpath), boost::filesystem::directory_iterator(), std::back_inserter(sortedEntries) );
     sort(sortedEntries.begin(), sortedEntries.end() );
     unsigned numSubs  = sortedEntries.size();

     printf("InitSubSamples(): numSubs = %i\n", numSubs);

     // Allocate samples memory space.
     ImageType3DChar::IndexType start;
     start.Fill(0);
     ImageType3DChar::SizeType sampleSize = maskPtr->GetLargestPossibleRegion().GetSize();
     ImageType3DChar::RegionType sampleRegion;
     sampleRegion.SetSize(sampleSize);
     sampleRegion.SetIndex(start);

     // source.
     ReaderType3DChar::Pointer srcReader = ReaderType3DChar::New();
     ImageType3DChar::Pointer srcPtr;

     for (unsigned subIdx = 0; subIdx < sampleVec.size(); subIdx ++) {
	  // source.
	  std::cout << "Init sub " << subIdx + 1 << " all samples with " << sortedEntries[subIdx].string() << std::endl;
	  srcReader->SetFileName( sortedEntries[subIdx].string() );
	  srcReader->Update();
	  srcPtr = srcReader->GetOutput();
	  IteratorType3DChar srcIt(srcPtr, srcPtr->GetLargestPossibleRegion() );

	  for (unsigned mcIdx = 0; mcIdx < sampleVec[0].size(); mcIdx ++) {
	       sampleVec[subIdx][mcIdx] = ImageType3DChar::New();
	       sampleVec[subIdx][mcIdx]->SetRegions(sampleRegion);
	       sampleVec[subIdx][mcIdx]->Allocate();
	       sampleVec[subIdx][mcIdx]->FillBuffer(-1);

	       // destination.
	       IteratorType3DChar sampleIt(sampleVec[subIdx][mcIdx],
					   sampleVec[subIdx][mcIdx]->GetLargestPossibleRegion() );
	       
	       for (maskIt.GoToBegin(), srcIt.GoToBegin(), sampleIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt, ++ srcIt, ++ sampleIt) {
		    if (maskIt.Get() > 0) {
			 sampleIt.Set( srcIt.Get() -1 );
		    }
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

     // For now do not use mask. Since initGrpPtr have -1 outside of mask, and
     // this will be copied to groupVec unchanged. (May use mask in the future.
     for (unsigned mcIdx = 0; mcIdx < grpVec.size(); mcIdx ++) {
	  grpVec[mcIdx] = ImageType3DChar::New();
	  grpVec[mcIdx]->SetRegions(grpRegion);
	  grpVec[mcIdx]->Allocate();
	  grpVec[mcIdx]->FillBuffer(-1);

	  grpVec[mcIdx]->SetOrigin(initGrpPtr->GetOrigin());
	  grpVec[mcIdx]->SetSpacing(initGrpPtr->GetSpacing());
	  grpVec[mcIdx]->SetDirection(initGrpPtr->GetDirection());
	       
	  IteratorType3DChar grpIt(grpVec[mcIdx],
				   grpVec[mcIdx]->GetLargestPossibleRegion() );
	       
	  for (initGrpIt.GoToBegin(), grpIt.GoToBegin(); !grpIt.IsAtEnd(); ++ grpIt, ++ initGrpIt) {
	       grpIt.Set( initGrpIt.Get() );
	  }
     } // mcIdx
     return 0;
}

// used by segimagealt. for allocating memory, initialize, and also used as mask.
int InitSamples(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec, 
		ImageType3DChar::Pointer grpPtr,
		ImageType3DChar::Pointer maskPtr)
{
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

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
	       sampleVec[subIdx][mcIdx]->FillBuffer(-1);

	       sampleVec[subIdx][mcIdx]->SetOrigin(grpPtr->GetOrigin());
	       sampleVec[subIdx][mcIdx]->SetSpacing(grpPtr->GetSpacing());
	       sampleVec[subIdx][mcIdx]->SetDirection(grpPtr->GetDirection());

	       
	       IteratorType3DChar sampleIt(sampleVec[subIdx][mcIdx],
					   sampleVec[subIdx][mcIdx]->GetLargestPossibleRegion() );
	       
	       for (maskIt.GoToBegin(), grpIt.GoToBegin(), sampleIt.GoToBegin(); !grpIt.IsAtEnd(); ++maskIt, ++ grpIt, ++ sampleIt) {
		    if (maskIt.Get() > 0) {
			 sampleIt.Set( grpIt.Get() );
		    }
	       }
	  } // mcIdx
     } // subIdx.
     return 0;
}
