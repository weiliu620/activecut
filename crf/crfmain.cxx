#include <common.h>
#include <utility.h>

twister_base_gen_type mygenerator(42u);

int main(int argc, char* argv[])
{
     ParStruct par;
     unsigned seed = 0;
     unsigned int emIter = 1;

     std::string dataFile, iLabelFile, oLabelFile, sampleFile;
     std::string initsubpath;
     std::string fmriPath, initGrplabel, cumusampledir, rsampledir;

     bool estprior = false;
     bool initsame = false;

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Inference by conditional random field.")
	  // ("nthreads", po::value<unsigned>(&par.nthreads)->default_value(100),
	  //  "number of threads to run for sampling")
	  ("burnin,b", po::value<unsigned>(&par.burnin)->default_value(20),
	   "number of scans for burn-in period. ")
	  ("sweepperthread", po::value<unsigned>(&par.sweepPerThread)->default_value(1),
	   "number of scans for burn-in period. ")
	  ("inittemp", po::value<float>(&par.initTemp)->default_value(1),
	   "Initial temperature for annealing.")
	  ("finaltemp", po::value<float>(&par.finalTemp)->default_value(1),
	   "Final temperature for annealing to stop")
	  ("numSamples,n", po::value<unsigned int>(&par.numSamples)->default_value(20),
	   "Number of Monte Carlo samples. ")
	  ("emiter", po::value<unsigned int>(&emIter)->default_value(30),
	   "Number of EM iteration. ")

	  ("alpha", po::value<float>(&par.alpha)->default_value(0.7),
	   "connection between group lebel and individual subject level.")

	  ("beta", po::value<float>(&par.beta)->default_value(0.5),
	   "pairwise interation term")

	  ("gamma", po::value<float>(&par.gamma)->default_value(1),
	   " weights of lilelihood term.")

	  ("alpha0", po::value<float>(&par.alpha0)->default_value(0.7),
	   "Mean of prior on alpha.")

	  ("sigma2", po::value<float>(&par.sigma2)->default_value(1),
	   "variance of prior on alpha")

	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")

	  ("numClusters,k", po::value<unsigned>(&par.numClusters)->default_value(6),
	   "Number of clusters. Default is 6.")

	  ("estprior", po::bool_switch(&estprior), 
	   "whether to estimate alpha and beta, the prior parameter. Default no.")

	  ("initsame", po::bool_switch(&initsame), 
	   "whether to init subject same to group label. Default is no.")

	  ("weightbetadata", po::bool_switch(&par.weightbetadata), 
	   "whether to weight beta by voxel correlation. Default is no.")

	  ("weightbetadist", po::bool_switch(&par.weightbetadist), 
	   "whether to weight beta by voxel's spatial distance. Default is no.")

	   ("initgrouplabel,i", po::value<std::string>(&initGrplabel)->default_value("grouplabel.nii"), 
	    "Initial group level label map (Intensity value 1-K. Also used as mask file.")

	   ("initsubpath,t", po::value<std::string>(&initsubpath)->default_value("./subpath"), 
	    "Initial subject label maps path. The files in this path must has the same file names as the fMRI files.")

	  ("fmripath,f", po::value<std::string>(&fmriPath)->default_value("."), 
	   "noised image file")

	  ("cumusampledir", po::value<std::string>(&cumusampledir)->default_value("./cumusamples"), 
	   "Monte Carlo samples file dir.")

	  ("rsampledir", po::value<std::string>(&rsampledir)->default_value("./rsamples"), 
	   "running samples file dir. Also can be the ultimate output dir.")

	  ("verbose,v", po::value<unsigned short>(&par.verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: gropumrf [options]\n";
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

     // output label init'd.
     ImageType3DChar::Pointer labelPtr = myfilter->GetOutput();

     // original 1-based group map only used for mask.
     ImageType3DChar::Pointer maskPtr = grpReader->GetOutput();
     par.maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     par.cumusampledir = cumusampledir;
     par.rsampledir = rsampledir;

     // par.numSubs = ComputeNumSubs(fmriPath);
     printf("number of subjects: %i.\n", par.numSubs);
     par.temperature = par.initTemp;

     // init parameters so BuildDataMap knows what to do.
     par.vmm.sub.resize(par.numSubs);
}

unsigned ComputeNumSubs(std::string srcpath)
{
     // prepare for sorting directory entries.
     typedef std::vector<boost::filesystem::path> PathVec;
     PathVec  sortedEntries;
     copy(boost::filesystem::directory_iterator(srcpath), boost::filesystem::directory_iterator(), std::back_inserter(sortedEntries) );
     sort(sortedEntries.begin(), sortedEntries.end() );

     // Get number of subjects.
     unsigned numSubs  = sortedEntries.size();
     return numSubs;
}


vnl_sparse_matrix <double> & con_map,

int BuildGraph(vnl_sparse_matrix <double> & G,
	       lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
	       unsigned numSubs, 
	       ImageType3DChar::Pointer maskPtr)

{
     //mask
     IteratorType3DCharIdx maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     ImageType3DChar::IndexType neiIdx;
     
     // map xyz to graph node id.
     ImageType4Int::IndexType nodeMapIdx;
     nodeMapIdx.Fill(0);
     ImageType3DChar::SizeType maskSize =  maskPtr->GetLargestPossibleRegion().GetSize();
     ImageType4Int::SizeType nodeMapSize;
     nodeMapSize[0] = maskSize[0];
     nodeMapSize[1] = maskSize[1];
     nodeMapSize[2] = maskSize[2];
     nodeMapSize[3] = numSubs + 1;

     ImageType4Int::RegionType nodeMapRegion;
     nodeMapRegion.SetSize(nodeMapSize);
     nodeMapRegion.SetIndex(nodeMapIdx);
     ImageType4Int::Pointer nodeMapPtr = ImageType4Int::New();
     nodeMapPtr->SetRegions(nodeMapRegion);
     nodeMapPtr->Allocate();
     nodeMapPtr->FillBuffer(0);


     // comptue number of points in one subjects. 
     unsigned sub_pts = 0;
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       sub_pts ++;
	  }	  
     }

     // we need to set the sparse matrix size here. Otherwise the assigment will
     // fail. Not sure why this happens. (ideally it should just add an entry).

     unsigned n_pts = sub_pts * (numSubs + 1); // including group level.
     G.set_size(n_pts, n_pts);

     // fill in subjects nodes.
     for (unsigned short subIdx = 0; subIdx < numSubs; subIdx ++) {
	  for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	       if (maskIt.Get() > 0) {
		    curNode = theGraph.addNode();
		    coordMap[curNode].idx = maskIt.GetIndex();
		    coordMap[curNode].subid = subIdx;
		    nodeMapIdx[0] = coordMap[curNode].idx[0];
		    nodeMapIdx[1] = coordMap[curNode].idx[1];
		    nodeMapIdx[2] = coordMap[curNode].idx[2];
		    nodeMapIdx[3] = subIdx;
		    nodeMapPtr->SetPixel(nodeMapIdx, theGraph.id(curNode) );
	       }
	  }
     }
