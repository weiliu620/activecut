#include <common.h>
#include <utility.h>

twister_base_gen_type mygenerator(42u);
using namespace lemon;

int BuildGraph(lemon::SmartGraph & G, 
	       lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coord,
	       ImageType3DChar::Pointer maskPtr);

int BuildDataMap(lemon::SmartGraph & G, 
		 lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coord,
		 lemon::SmartGraph::NodeMap<vnl_vector<float> > & tsMap,
		 std::string fmriImageFile,
		 ParStruct & par);

int EstimateMu(lemon::SmartGraph & G, 
	       lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coord,
	       lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
	       lemon::SmartGraph::NodeMap<vnl_vector<float> > & tsMap,
	       ParStruct & par);

int EstimateKappa(ParStruct & par);

int Sampling(lemon::SmartGraph & G, 
	     lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coord,
	     lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
	     lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
	     lemon::SmartGraph::NodeMap<vnl_vector<float> > & tsMap,
	     ParStruct & par);

int CompSampleEnergy(lemon::SmartGraph & G, 
		     lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coord,
		     lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
		     lemon::SmartGraph::NodeMap<vnl_vector<float> > & tsMap,
		     ParStruct & par);

int main(int argc, char* argv[])
{
     ParStruct par;
     unsigned seed = 0;
     unsigned int emIter = 1;

     std::string initLabelFile, fmriImageFile, cumuSampleFile, rSampleFile;
     bool estprior = false;
     bool initsame = false;

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
	   "Number of Monte Carlo samples. ")
	  ("emiter", po::value<unsigned int>(&emIter)->default_value(30),
	   "Number of EM iteration. ")
	  ("beta", po::value<float>(&par.beta)->default_value(0.5),
	   "pairwise interation term")
	  ("gamma", po::value<float>(&par.gamma)->default_value(1),
	   " weights of lilelihood term.")
	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")

	  ("numClusters,k", po::value<unsigned>(&par.numClusters)->default_value(6),
	   "Number of clusters. Default is 6.")
	  ("estprior", po::bool_switch(&estprior), 
	   "whether to estimate alpha and beta, the prior parameter. Default no.")
	   ("init,i", po::value<std::string>(&initLabelFile)->default_value("initlabelfile.nii.gz"), 
	    "Initial label map (Intensity value 1-K. Also used as mask file.")
	  ("fmri,f", po::value<std::string>(&fmriImageFile)->default_value("fmriimage.nii.gz"), 
	   "noised image file")

	  ("cumusample,c", po::value<std::string>(&cumuSampleFile)->default_value("cumusamplefile.nii.gz"), 
	   "Cumulative sample file.")
	  ("rsample", po::value<std::string>(&rSampleFile)->default_value("rsamplefile.nii.gz"), 
	   "running samples file.")

	  ("verbose,v", po::value<unsigned short>(&par.verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: mcem_lemon [options]\n";
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

     // read in initial label map. Also used as mask. 1-based label converted to
     // 0-based label.
     ReaderType3DChar::Pointer initReader = ReaderType3DChar::New();
     initReader->SetFileName(initLabelFile);
     AddConstantToImageFilterType::Pointer myfilter = AddConstantToImageFilterType::New();
     myfilter->InPlaceOff();
     myfilter->SetInput( initReader->GetOutput());
     myfilter->SetConstant(-1);
     myfilter->Update();

     // output label init'd.
     ImageType3DChar::Pointer labelPtr = myfilter->GetOutput();
     // original 1-based group map only used for mask.
     ImageType3DChar::Pointer maskPtr = initReader->GetOutput();

     par.maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     
     // Build the graph.
     lemon::SmartGraph G;
     lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> coord(G);

     BuildGraph(G, coord, maskPtr);

     // Build Data map. 
     lemon::SmartGraph::NodeMap<vnl_vector<float> > tsMap(G);
     BuildDataMap(G, coord, tsMap, fmriImageFile, par);
     
     // At this point, we know time series length so do the remaining
     // initializaion of par.vmm.
     par.vmm.comp.resize(par.numClusters);
     for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	  par.vmm.comp[clsIdx].mu.set_size(par.tsLength);
	  par.vmm.comp[clsIdx].mu = 0;
	  par.vmm.comp[clsIdx].meanNorm = 0;
	  par.vmm.comp[clsIdx].kappa = 0;
	  par.vmm.comp[clsIdx].numPts = 0;
	  par.vmm.comp[clsIdx].prop = 0;
     }

     // define cumulative sample map. Each node is assigned a vector. The k'th
     // elements of the vector is the number of samples that falls into cluster
     // k.
     lemon::SmartGraph::NodeMap<std::vector<unsigned short> > cumuSampleMap(G);
     for (SmartGraph::NodeIt nodeIt(G); nodeIt !=INVALID; ++ nodeIt) {
     	  cumuSampleMap[nodeIt].resize(par.numClusters, 0);
     }

     // init the cumuSampleMap just to estimate parameters. The k'th bit is set
     // to the number of samples M, since I suppose the M samples are init'd
     // same to the initial group or subject label map.
     for (SmartGraph::NodeIt nodeIt(G); nodeIt !=INVALID; ++ nodeIt) {
	  cumuSampleMap[nodeIt].at(labelPtr->GetPixel(coord[nodeIt])) = par.numSamples;
     }

     SaveCumuSample(G, coord, cumuSampleMap, par, cumuSampleFile);

     // define a running sample map and init it with the input labels. This is
     // used for sampling only. And when the whoel scan is done, this map is
     // saved into the cumuSampleMap.
     lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > rSampleMap(G);
     for (SmartGraph::NodeIt nodeIt(G); nodeIt !=INVALID; ++ nodeIt) {
	  rSampleMap[nodeIt].resize(par.numClusters, 0);
	  rSampleMap[nodeIt].set(labelPtr->GetPixel(coord[nodeIt]), true);
     }

     EstimateMu(G, coord, cumuSampleMap, tsMap, par);
     EstimateKappa(par);	  

     // EM starts here.
     for (unsigned short emIterIdx = 0; emIterIdx < emIter; emIterIdx ++) {     
     	  printf("EM iteration %i begin:\n", emIterIdx + 1);
	  par.temperature = par.initTemp * pow( (par.finalTemp / par.initTemp), float(emIterIdx) / float(emIter) );
	  Sampling(G, coord, cumuSampleMap, rSampleMap, tsMap, par);
	  SaveCumuSample(G, coord, cumuSampleMap, par, cumuSampleFile);

	  PrintPar(1, par);
	  CompSampleEnergy(G, coord, rSampleMap, tsMap, par);

     	  // estimate vMF parameters mu, kappa.
     	  printf("EM iteration %i, parameter estimation begin. \n", emIterIdx + 1);
	  EstimateMu(G, coord, cumuSampleMap, tsMap, par);
	  EstimateKappa(par);
	  PrintPar(1, par);
     }

     return 0;
}

int BuildGraph(lemon::SmartGraph & G, 
	       lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coord,
	       ImageType3DChar::Pointer maskPtr)
{
     //mask
     IteratorType3DCharIdx maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     
     // a data structure to map xyz coordinate to node id.
     ImageType3DU::IndexType nodeMapIdx;
     nodeMapIdx.Fill(0);
     ImageType3DU::SizeType nodeMapSize = maskPtr->GetLargestPossibleRegion().GetSize();

     ImageType3DU::RegionType nodeMapRegion;
     nodeMapRegion.SetSize(nodeMapSize);
     nodeMapRegion.SetIndex(nodeMapIdx);
     ImageType3DU::Pointer nodeMapPtr = ImageType3DU::New();
     nodeMapPtr->SetRegions(nodeMapRegion);
     nodeMapPtr->Allocate();
     nodeMapPtr->FillBuffer(0);

     // Add nodes.

     lemon::SmartGraph::Node curNode, neiNode;
     
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
     	  if (maskIt.Get() > 0) {
	       curNode = G.addNode();
	       coord[curNode] = maskIt.GetIndex();
	  
	       nodeMapIdx = maskIt.GetIndex();
	       nodeMapPtr->SetPixel(nodeMapIdx, G.id(curNode));
	  }
     }
     printf("BuildGraph(): number of Nodes: %i\n", G.maxNodeId()+1 );

     // Add edges.
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType maskNeiIt(radius, maskPtr, maskPtr->GetLargestPossibleRegion() );

     std::array<unsigned int, 6 > neiIdxSet = {{14, 12, 16, 10, 22, 4}}; 
     // std::array<unsigned int, 26> neiIdxSet = {{ 0,1,2,3,4,5,6,7,8,9,10,11,12,//no 13
     // 						 14,15,16,17,18,19,20,21,22,23,24,25,26 }};
     // std::array<unsigned int, 18 > neiIdxSet = {{1,3,4,5,7,9,10,11,12,
     // 						14,15,16,17,19,21,22,23,25}};
     
     int curNodeId = 0, neiNodeId = 0;
     auto neiIdxIt = neiIdxSet.begin();

     for (maskNeiIt.GoToBegin(); !maskNeiIt.IsAtEnd(); ++ maskNeiIt) {
	  if (maskNeiIt.GetCenterPixel() > 0) {
	       curNodeId = nodeMapPtr->GetPixel( maskNeiIt.GetIndex() );
	       curNode = G.nodeFromId( curNodeId );
	       for (neiIdxIt=neiIdxSet.begin(); neiIdxIt < neiIdxSet.end(); neiIdxIt++) {
		    if (maskNeiIt.GetPixel(*neiIdxIt) > 0) {
			 neiNodeId = nodeMapPtr->GetPixel( maskNeiIt.GetIndex(*neiIdxIt) );
			 // make sure each edge is added only once.
			 if (neiNodeId > curNodeId) {
			      neiNode = G.nodeFromId( neiNodeId );
			      G.addEdge(curNode, neiNode);
			 } // if
		    } // maskIt
	       } // neiIt
	  } // maskNeiIt > 0
     } // maskNeiIt
     printf("BuildGraph(): Number of edges: %i\n", G.maxEdgeId() + 1);
     return 0;
}


int BuildDataMap(lemon::SmartGraph & G, 
		 lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coord,
		 lemon::SmartGraph::NodeMap<vnl_vector<float> > & tsMap,
		 std::string fmriImageFile,
		 ParStruct & par)

{
     ReaderType4DFloat::Pointer fmriReader = ReaderType4DFloat::New();
     fmriReader->SetFileName( fmriImageFile );
     fmriReader->Update();
     ImageType4DF::Pointer fmriPtr = fmriReader->GetOutput();
     ImageType4DF::SizeType fmriSize = fmriPtr->GetLargestPossibleRegion().GetSize();
     ImageType4DF::IndexType fmriIdx;
     par.tsLength = fmriSize[3];

     for (SmartGraph::NodeIt nodeIt(G); nodeIt !=INVALID; ++ nodeIt) {
	  tsMap[nodeIt].set_size(par.tsLength);
	  fmriIdx[0] = coord[nodeIt][0];
	  fmriIdx[1] = coord[nodeIt][1];
	  fmriIdx[2] = coord[nodeIt][2];
	  
	  for (fmriIdx[3] = 0; fmriIdx[3] < par.tsLength; fmriIdx[3] ++) {
	       tsMap[nodeIt][fmriIdx[3]] = fmriPtr->GetPixel(fmriIdx);
	  }
     }
     return 0;
}


int EstimateMu(lemon::SmartGraph & G, 
	       lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coord,
	       lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
	       lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
	       ParStruct & par)
{
     // reset all mu and numPts to zero.
     for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	  par.vmm.comp[clsIdx].mu = 0;
	  par.vmm.comp[clsIdx].numPts = 0;
     }

     for (SmartGraph::NodeIt nodeIt(G); nodeIt !=INVALID; ++ nodeIt) {
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       par.vmm.comp[clsIdx].mu += tsMap[nodeIt] * cumuSampleMap[nodeIt][clsIdx];
	       par.vmm.comp[clsIdx].numPts += cumuSampleMap[nodeIt][clsIdx];
	  }
     }

     for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {     
	  if(par.vmm.comp[clsIdx].numPts > 0) {
	       par.vmm.comp[clsIdx].meanNorm = par.vmm.comp[clsIdx].mu.two_norm() / par.vmm.comp[clsIdx].numPts;
	       par.vmm.comp[clsIdx].mu.normalize();
	  }
	  else {
	       par.vmm.comp[clsIdx].meanNorm = 0;
	  }
     }
	       
     printf("EstimateMu(). Done.\n");
     return 0;
}


int EstimateKappa(ParStruct & par)
{
     float kappa = 0, kappa_new = 0;
     double Ap = 0;
     float RBar = 0;
     float Dim = par.tsLength;

     for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	  if (par.vmm.comp[clsIdx].numPts > 0){

	       RBar = par.vmm.comp[clsIdx].meanNorm;
	       kappa_new = RBar * (Dim - RBar * RBar) / (1 - RBar * RBar);
	       unsigned iter = 0;
	       do {
		    iter ++;
		    kappa = kappa_new;
		    Ap = exp(logBesselI(Dim/2, kappa) - logBesselI(Dim/2 - 1, kappa));
		    kappa_new = kappa - ( (Ap - RBar) / (1 - Ap * Ap - (Dim - 1) * Ap / kappa)  );
		    if (par.verbose >= 3) {
			 printf("    cls[%i] kappa: %3.1f -> %3.1f\n", clsIdx + 1, kappa, kappa_new);
		    }
	       }
	       while(vnl_math_abs(kappa_new - kappa) > 0.01 * kappa && iter < 5);
	       par.vmm.comp[clsIdx].kappa = kappa_new;
	  } // numPts > 0
	  else {
	       par.vmm.comp[clsIdx].kappa = 0;
	  }

     } // clsIdx

     printf("EstimateKappa(). Done.\n");
     return 0;
}


int Sampling(lemon::SmartGraph & G, 
	     lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coord,
	     lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
	     lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
	     lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
	     ParStruct & par)

{
     // reset cumuSampleMap to zero.
     for (SmartGraph::NodeIt nodeIt(G); nodeIt !=INVALID; ++ nodeIt) {
	  cumuSampleMap[nodeIt].assign(par.numClusters, 0);
     }

     // define random generator.
     boost::random::mt19937 gen(std::time(0));
     boost::random::uniform_int_distribution<> roll_die(0, par.numClusters - 1);

     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(gen, uni_dist);

     std::vector<double> vmfLogConst(par.numClusters);
     float myD = par.tsLength;
     double const Pi = 4 * atan(1);
     unsigned clsIdx = 0;
     for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	  if (par.vmm.comp[clsIdx].kappa > 1) {
	       vmfLogConst[clsIdx] = (myD/2 - 1) * log (par.vmm.comp[clsIdx].kappa) - myD/2 * log(2*Pi) -  logBesselI(myD/2 - 1, par.vmm.comp[clsIdx].kappa);
	  }
	  else {
	       vmfLogConst[clsIdx] = 25.5; // the limit when kappa -> 0
	  }
     } // for clsIdx

     float oldPriorEngy = 0, newPriorEngy = 0, oldDataEngy = 0, newDataEngy = 0, denergy;
     boost::dynamic_bitset<> cand;
     cand.resize(par.numClusters);
     float p_acpt = 0;
     // cl is current label, and nl is new label.
     unsigned short cl = 0, nl = 0;
     lemon::SmartGraph::Node curNode;
     unsigned flipEngy = 0, numUpdLabels = 0;;

     // Sampling starts here.
     printf("Sampling() starts scan...\n");
     for (unsigned scanIdx = 0; scanIdx < par.burnin + par.numSamples; scanIdx ++) {
	  printf("%i, ", scanIdx+1);
	  fflush(stdout);
	  flipEngy = 0, numUpdLabels = 0;;
	  for (SmartGraph::NodeIt nodeIt(G); nodeIt !=INVALID; ++ nodeIt) {
	       // compute old and new prior energy.
	       curNode = nodeIt;
	       oldPriorEngy = 0;
	       newPriorEngy = 0;
	       cand.reset();
	       cand[roll_die(gen)] = 1;

	       if (par.verbose >= 3) {
	       printf("curNode:[%ld %ld %ld]\n", coord[curNode][0], coord[curNode][1], coord[curNode][2]);
	       }

	       for (SmartGraph::IncEdgeIt edgeIt(G, curNode); edgeIt != INVALID; ++ edgeIt) {
		    if (par.verbose >= 3) {
			 printf("  edgeId: %d, neiNodeId: %d, [%ld %ld %ld]\n", G.id(edgeIt), G.id(G.runningNode(edgeIt)), coord[G.runningNode(edgeIt)][0], coord[G.runningNode(edgeIt)][1], coord[G.runningNode(edgeIt)][2]);
		    }

	       	    oldPriorEngy += par.beta * (double)(rSampleMap[G.baseNode(edgeIt)] != rSampleMap[G.runningNode(edgeIt)]);
	       	    newPriorEngy += par.beta * (double)(cand != rSampleMap[G.runningNode(edgeIt)]);
	       }
	       denergy = newPriorEngy - oldPriorEngy;
	       
	       cl = rSampleMap[curNode].find_first(); // current label.
	       nl = cand.find_first(); // new label.
	       oldDataEngy = - vmfLogConst[cl] - par.vmm.comp[cl].kappa * inner_product(tsMap[curNode], par.vmm.comp[cl].mu);
	       newDataEngy = - vmfLogConst[nl] - par.vmm.comp[nl].kappa * inner_product(tsMap[curNode], par.vmm.comp[nl].mu);
	       
	       denergy = denergy + par.gamma * (newDataEngy - oldDataEngy);

	       if (denergy * par.gamma * (newDataEngy - oldDataEngy) < 0) {
		    // that means after introducing prior, the posterior engy changed sign compared to likelihood.
		    flipEngy ++;
	       }
	       denergy = denergy / par.temperature;


	       if (denergy < 0) {
		    rSampleMap[curNode] = cand;
		    numUpdLabels ++;
	       }
	       else if (denergy > 0) {
		    p_acpt = exp(-denergy);
		    if (uni() < p_acpt) {
			 rSampleMap[curNode] = cand;
			 numUpdLabels ++;
		    }
	       } // else if
	       else {
		    // candidate = current label. No change.
	       }
	  } // nodeIt

	  if (par.verbose >=3) {
	       printf("scan: %d, flipEngy: %d, numUpdLabels: %d\n", scanIdx+1, flipEngy, numUpdLabels);
	  }

	  // after burnin peroid, save it to correct place.
	  if (scanIdx >= par.burnin) {
	       for (SmartGraph::NodeIt nodeIt(G); nodeIt !=INVALID; ++ nodeIt) {
		    cumuSampleMap[nodeIt][rSampleMap[nodeIt].find_first()] ++;
	       }
	  }

     } // scanIdx
     printf("Done.\n");
     return 0;
}

int CompSampleEnergy(lemon::SmartGraph & G, 
		     lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coord,
		     lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
		     lemon::SmartGraph::NodeMap<vnl_vector<float>> & tsMap,
		     ParStruct & par)
{
     double samplePriorEng = 0, sampleDataEng = 0;
     boost::dynamic_bitset<> curBit(par.numClusters), runningBit(par.numClusters);
     double c = 0, ccur = 0, M0 = 0;
     unsigned runningLabel = 0;
     for (SmartGraph::NodeIt nodeIt(G); nodeIt !=INVALID; ++ nodeIt) {
	  curBit = rSampleMap[nodeIt];
	  M0 = 0;
	  for (runningLabel = 0; runningLabel < par.numClusters; runningLabel ++) {
	       runningBit.reset(), runningBit[runningLabel] = 1;
	       c = 0;
	       for (SmartGraph::IncEdgeIt edgeIt(G, nodeIt); edgeIt != INVALID; ++ edgeIt) {
		    c = c - par.beta * (double)(runningBit != rSampleMap[G.runningNode(edgeIt)]) ;
	       } // incEdgeIt
	       M0 += exp (c);

	       if (runningBit == curBit) {
		    ccur = c;
	       }
	  } // runningBit.
	  samplePriorEng += (ccur - log(M0));
     } // nodeIt

     // convert log likelihood to energy.
     samplePriorEng = - samplePriorEng;

     std::vector<double> vmfLogConst(par.numClusters);
     float myD = par.tsLength;
     double const Pi = 4 * atan(1);
     unsigned clsIdx = 0;
     for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	  if (par.vmm.comp[clsIdx].kappa > 1) {
	       vmfLogConst[clsIdx] = (myD/2 - 1) * log (par.vmm.comp[clsIdx].kappa) - myD/2 * log(2*Pi) -  logBesselI(myD/2 - 1, par.vmm.comp[clsIdx].kappa);
	  }
	  else {
	       vmfLogConst[clsIdx] = 25.5; // the limit when kappa -> 0
	  }
     } // for clsIdx

     for (SmartGraph::NodeIt nodeIt(G); nodeIt !=INVALID; ++ nodeIt) {
	  clsIdx = rSampleMap[nodeIt].find_first();
	  sampleDataEng += vmfLogConst[clsIdx] + par.vmm.comp[clsIdx].kappa * inner_product(tsMap[nodeIt], par.vmm.comp[clsIdx].mu);
     } // nodeIt

     // convert log likelihood to energy.
     sampleDataEng = - sampleDataEng;

     printf("samplePriorEng: %E, sampleDataEng: %E, single sample energy (- log P(X,Y) ): %E.\n", samplePriorEng, sampleDataEng, samplePriorEng + sampleDataEng);
     return 0;
}
