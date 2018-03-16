#include "commonalt.h"
#include "vmmodel.h"
#include <utilalt.h>

int InitGrp(std::vector< ImageType3DChar::Pointer > & grpVec,
	    ImageType3DChar::Pointer initGrpPtr);

int SplitMRI(ImageType5DFloat::Pointer inPtr,      
	     std::vector<VnlVectorImageType::Pointer> & fmriVec);

using std::cout;
using std::endl;

twister_base_gen_type mygenerator(42u);


int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     int burnin = 1;
     unsigned seed = 0;

     unsigned int numSamples = 1;
     unsigned int EMIter = 1;
     int numClusters = 4;
     float initTemp = 0.5;
     float finalTemp = 0.01;
     float beta_g = 0.5, gamma = 100;

     std::string fmriPath, initGrplabel, outGrpLabel, outGrpProb, outSubBase;
     std::string grpSampleFilename;
     std::string samplingMethod;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Inference on group mrf.")
	  ("burnin,b", po::value<int>(&burnin)->default_value(20),
	   "number of scans for burn-in period. ")
	  ("inittemp", po::value<float>(&initTemp)->default_value(1),
	   "Initial temperature for annealing.")
	  ("finaltemp", po::value<float>(&finalTemp)->default_value(1),
	   "Final temperature for annealing to stop")
	  ("numSamples,n", po::value<unsigned int>(&numSamples)->default_value(20),
	   "initial number of Monte Carlo samples. ")
	  ("emiter", po::value<unsigned int>(&EMIter)->default_value(30),
	   "initial number of Monte Carlo samples. ")
	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("numClusters,k", po::value<int>(&numClusters)->default_value(6),
	       "Number of clusters. Default is 6.")
	  ("betag", po::value<float>(&beta_g)->default_value(0.5),
	   "group level pairwise interation term")
	  ("gamma", po::value<float>(&gamma)->default_value(100),
	   "likelihood normalization term")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")

	  ("fmripath,f", po::value<std::string>(&fmriPath)->default_value("."),
	   "noised image file")

	  ("grouplabel,g", po::value<std::string>(&outGrpLabel)->default_value("outgrplabel.nii"), "output group label file. nii format.")

	  ("grpsamples", po::value<std::string>(&grpSampleFilename)->default_value("grpSample.nii.gz"), "group label Monte Carlo samples file")


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

	   
     boost::filesystem::path fmriPathVar(fmriPath);
     boost::filesystem::directory_iterator fmriPathEnd;

     ImageType5DFloat::IndexType imageIdx;          

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


     // read in initial grp label. Also used as mask.
     std::vector<ImageType3DChar::Pointer> grpVec(numSamples);
     
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

     // for (unsigned mcIdx = 0; mcIdx < numSamples; mcIdx ++) {
     // 	  std::string mcString = boost::lexical_cast<std::string> (mcIdx+1);
     // 	  std::string thisFilename = samplePrefix;
     // 	  thisFilename.append(mcString);
     // 	  thisFilename.append(".nii.gz");
     // 	  save3dcharInc(grpVec[mcIdx], thisFilename);
     // }

     SaveGrpSamples(grpVec, maskPtr, grpSampleFilename);
     

     VMModel vmmodel(fmriVec,
		     grpVec,
		     maskPtr,
		     numClusters,
		     beta_g,
		     gamma,
     		     verbose);


     vmmodel.EstimateMu(grpVec, fmriVec, maskPtr);
     vmmodel.estimateKappa(fmriSize[3]);
     vmmodel.printSelf("normal");
     
     for (unsigned short emIterIdx = 0; emIterIdx < EMIter; emIterIdx ++) {

	  if (verbose >= 1) {
	       printf("EM iteration %i begin:\n", emIterIdx);
	  }

	  vmmodel.SetTemperature(initTemp * pow( (finalTemp / initTemp), float(emIterIdx) / float(EMIter) ) );

	  vmmodel.mcsampling(grpVec, fmriVec, maskPtr, burnin);


	  if (verbose >= 1) {
	       save3dcharInc(grpVec[numSamples -1], outGrpLabel);
	       SaveGrpSamples(grpVec, maskPtr, grpSampleFilename);
	       printf("EM iteration %i, sampling done. \n", emIterIdx);
	  }

	  // estimate vMF parameters mu, kappa.
	  vmmodel.EstimateMu(grpVec, fmriVec, maskPtr);
	  vmmodel.estimateKappa(fmriSize[3]);

	  vmmodel.ICM(grpVec, fmriVec, maskPtr, burnin);
	  save3dcharInc(grpVec[numSamples -1], outGrpLabel);

	  printf("EM iteration %i, parameter estimation done. \n", emIterIdx);
	  vmmodel.printSelf("normal");
     }

     vmmodel.ICM(grpVec, fmriVec, maskPtr, 10);

     // Save Group label map.
     save3dcharInc(grpVec[numSamples -1], outGrpLabel);
     
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
	  grpVec[mcIdx]->FillBuffer(0);
	       
	  IteratorType3DChar grpIt(grpVec[mcIdx],
				   grpVec[mcIdx]->GetLargestPossibleRegion() );
	       
	  for (initGrpIt.GoToBegin(), grpIt.GoToBegin(); !grpIt.IsAtEnd(); ++ grpIt, ++ initGrpIt) {
	       grpIt.Set( initGrpIt.Get() );
	  }
     } // mcIdx
     return 0;
}

// rearrange 5d fmri as a vector image.
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
