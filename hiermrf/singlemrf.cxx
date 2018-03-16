#include <common.h>
#include <utility.h>

twister_base_gen_type mygenerator(42u);

int BuildGraph(lemon::SmartGraph & theGraph, 
	       lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coordMap,
	       ImageType3DChar::Pointer maskPtr);

int BuildDataMap(lemon::SmartGraph & theGraph, 
		 lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coordMap,
		 lemon::SmartGraph::NodeMap<float > & dataMap,
		 std::string dataFile);


int EstimateMuSigma(lemon::SmartGraph & theGraph, 
		    lemon::SmartGraph::NodeMap<float > & dataMap,
		    lemon::SmartGraph::NodeMap< std::vector< unsigned char > > & labelMap,
		    SWParType & par);

int Sampling(lemon::SmartGraph & theGraph, 
	     lemon::SmartGraph::NodeMap<float > & dataMap,
	     lemon::SmartGraph::NodeMap< std::vector< unsigned char > > & labelMap,
	     SWParType & par);

int SWSampling(lemon::SmartGraph & theGraph, 
	       lemon::SmartGraph::NodeMap<float > & dataMap,
	       lemon::SmartGraph::NodeMap< std::vector< unsigned char > > & labelMap,
	       lemon::SmartGraph::EdgeMap<bool > & edgeMap,
	       SWParType & par);

int ComputeExtField(std::vector< std::vector<double> > & labelWeightSet,
		    lemon::SmartGraph::NodeMap<unsigned> & compMap,
		    lemon::SmartGraph & theGraph, 
		    lemon::SmartGraph::NodeMap< std::vector< unsigned char > > & labelMap,
		    lemon::SmartGraph::NodeMap<float > & dataMap,
		    SWParType & par);

unsigned short Delta(unsigned xi, unsigned xj);
int PrintParam(SWParType & par, unsigned plevel);


using namespace lemon;
int main(int argc, char* argv[])
{
     SWParType par;
     unsigned seed = 0;
     unsigned int emIter = 1;

     // std::string fmriPath, initGrplabel, outGrpLabel, outGrpProb, outSubBase, 
     // 	  grpSampleName, samplePrefix, initsubpath, groupprob, subbaseprob;
     std::string dataFile, iLabelFile, oLabelFile, sampleFile;

     bool usesw = false;

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Inference on group mrf.")
	  ("burnin,b", po::value<unsigned>(&par.burnin)->default_value(20),
	   "number of scans for burn-in period. ")
	  ("inittemp", po::value<double>(&par.initTemp)->default_value(1),
	   "Initial temperature for annealing.")
	  ("finaltemp", po::value<double>(&par.finalTemp)->default_value(1),
	   "Final temperature for annealing to stop")
	  ("numSamples,n", po::value<unsigned int>(&par.numSamples)->default_value(20),
	   "Number of Monte Carlo samples. ")
	  ("emiter", po::value<unsigned int>(&emIter)->default_value(30),
	   "Number of EM iteration. ")

	  ("beta", po::value<double>(&par.beta)->default_value(0.5),
	   "pairwise interation term")

	  ("gamma", po::value<double>(&par.gamma)->default_value(1),
	   " weights of lilelihood term.")

	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("numClusters,k", po::value<unsigned>(&par.numClusters)->default_value(6),
	   "Number of clusters. Default is 6.")

	  ("usesw,s", po::bool_switch(&usesw), "whether to use Swendsen-Wang for sampling.")

	   ("ilabelfile,i", po::value<std::string>(&iLabelFile)->default_value("iLabelFilenii.gz"), "Initiallabel map. Also used as mask file.")

	  ("datafile,f", po::value<std::string>(&dataFile)->default_value("datafile"), "data File. Either time series or scalar image file.")

	  ("samplefile", po::value<std::string>(&sampleFile)->default_value("subject"), "Monte Carlo samples file name prefix")

	  ("outlabel,g", po::value<std::string>(&oLabelFile)->default_value("oLabelFile.nii.gz"), "output label file. must be nii or nii.gz.")

	  ("verbose,v", po::value<unsigned short>(&par.verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: singlemrf [options]\n";
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

     // read in initial label. Also used as mask. 1-based label converted to
     // 0-based label.
     ReaderType3DChar::Pointer labelReader = ReaderType3DChar::New();
     labelReader->SetFileName(iLabelFile);
     AddConstantToImageFilterType::Pointer myfilter = AddConstantToImageFilterType::New();
     myfilter->InPlaceOff();
     myfilter->SetInput( labelReader->GetOutput());
     myfilter->SetConstant(-1);
     myfilter->Update();

     // output label init'd.
     ImageType3DChar::Pointer labelPtr = myfilter->GetOutput();

     // original 1-based group map only used for mask.
     ImageType3DChar::Pointer maskPtr = labelReader->GetOutput();

     // Construct graph.
     lemon::SmartGraph theGraph;
     lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> coordMap(theGraph); // node --> coordinates.
     BuildGraph(theGraph, coordMap, maskPtr);

     // Build Data map.
     lemon::SmartGraph::NodeMap< float > dataMap(theGraph); // node --> data.
     BuildDataMap(theGraph, coordMap, dataMap, dataFile);

     // define and init label map.
     lemon::SmartGraph::NodeMap<std::vector<unsigned char> > labelMap(theGraph);
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  labelMap[nodeIt].assign(par.numSamples, labelPtr->GetPixel(coordMap[nodeIt]));
     }

     // Define edge map.
     lemon::SmartGraph::EdgeMap<bool > edgeMap(theGraph);

     // init parameters.
     par.comp.resize(par.numClusters);
     par.totalPts = theGraph.maxNodeId()+1;
     par.temperature = 1;

     // estimate mu and sigma of Gaussian Mixture Model.
     EstimateMuSigma(theGraph, dataMap, labelMap, par);

     for (unsigned short emIdx = 0; emIdx < emIter; emIdx ++) {
	  if (par.verbose >= 1) {
	       printf("EM iteration %i begin:\n", emIdx +1 );
	  }

	  if (usesw) {
	       SWSampling(theGraph, dataMap, labelMap, edgeMap, par);
	  }
	  else {
	       Sampling(theGraph, dataMap, labelMap, par);
	  }

	  // estimate mu and sigma of Gaussian Mixture Model.
	  EstimateMuSigma(theGraph, dataMap, labelMap, par);

	  PrintParam(par, 1);
     }

     // save graph labels to image.
     SaveGraphToImage(theGraph, coordMap, labelMap, maskPtr, oLabelFile);

     // save samples to image.
     SaveSamplesToImage(theGraph, coordMap, labelMap, maskPtr, sampleFile);

     return 0;
	  
}

int BuildGraph(lemon::SmartGraph & theGraph, 
	       lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coordMap,
	       ImageType3DChar::Pointer maskPtr)
{
     //mask
     IteratorType3DCharIdx maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType maskIdx, neiIdx;
     
     // map xyz to graph node id.
     ImageType3DU::IndexType nodeMapIdx;
     nodeMapIdx.Fill(0);
     ImageType3DChar::SizeType maskSize =  maskPtr->GetLargestPossibleRegion().GetSize();
     ImageType3DU::SizeType nodeMapSize;
     nodeMapSize[0] = maskSize[0];
     nodeMapSize[1] = maskSize[1];
     nodeMapSize[2] = maskSize[2];

     ImageType3DU::RegionType nodeMapRegion;
     nodeMapRegion.SetSize(nodeMapSize);
     nodeMapRegion.SetIndex(nodeMapIdx);
     ImageType3DU::Pointer nodeMapPtr = ImageType3DU::New();
     nodeMapPtr->SetRegions(nodeMapRegion);
     nodeMapPtr->Allocate();
     nodeMapPtr->FillBuffer(0);
     
     // Add nodes.
     lemon::SmartGraph::Node curNode, neiNode;

     // Fill in group level nodes.
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
     	  if (maskIt.Get() > 0) {
     	       curNode = theGraph.addNode();
	       maskIdx = maskIt.GetIndex();
	       coordMap[curNode] = maskIdx;
	       
	       nodeMapIdx = coordMap[curNode];
	       nodeMapPtr->SetPixel(nodeMapIdx, theGraph.id(curNode) );
	  }
     }


     printf("BuildGraph(): number of Nodes: %i\n", theGraph.maxNodeId()+1 );
     printf("BuildGraph(): number of Nodes: %i\n", lemon::countNodes(theGraph) );

     // Add edges. 

     // Define neighborhood iterator
     MyBoundCondType constCondition;
     constCondition.SetConstant(-1);     
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType maskNeiIt(radius, maskPtr, maskPtr->GetLargestPossibleRegion() );
     maskNeiIt.OverrideBoundaryCondition(&constCondition);
     
     // xplus, xminus, yplus, yminus, zplus, zminus
     std::array<unsigned int, 6> neiIdxSet = {14, 12, 16, 10, 22, 4}; 

     int curNodeId = 0, neiNodeId = 0;
     auto neiIdxIt = neiIdxSet.begin();

     // edge within grp.
     for (maskNeiIt.GoToBegin(); !maskNeiIt.IsAtEnd(); ++ maskNeiIt) {
	  if (maskNeiIt.GetCenterPixel() > 0) {
	       maskIdx = maskNeiIt.GetIndex();
	       nodeMapIdx = maskIdx;
	       curNodeId = nodeMapPtr->GetPixel(nodeMapIdx);
	       // curNode is group node.
	       curNode = theGraph.nodeFromId( curNodeId );

	       for (neiIdxIt=neiIdxSet.begin(); neiIdxIt < neiIdxSet.end(); neiIdxIt++) {
		    if (maskNeiIt.GetPixel(*neiIdxIt) > 0) {
			 neiIdx = maskNeiIt.GetIndex(*neiIdxIt);
			 nodeMapIdx = neiIdx;
			 neiNodeId = nodeMapPtr->GetPixel(nodeMapIdx);
			 // make sure each edge is added only once.
			 if (neiNodeId > curNodeId) {
			      neiNode = theGraph.nodeFromId( neiNodeId );
			      theGraph.addEdge(curNode, neiNode);
			 } // if
		    } // maskIt
	       } // neiIt
	  } // maskNeiIt > 0
     } // maskNeiIt

     printf("BuildGraph(): number of edges: %i\n", theGraph.maxEdgeId() + 1);

     // test code. Test nodeMap and coordMap.
     
     for (int nodeId = 0; nodeId < 10; nodeId ++) {
	  maskIdx = coordMap[theGraph.nodeFromId( nodeId )];
	  printf("nodeId = %i, Idx = [%i, %i, %i]\n", nodeId, maskIdx[0], maskIdx[1], maskIdx[2]);
     }

     for (int nodeId = theGraph.maxNodeId() - 10; nodeId <= theGraph.maxNodeId(); nodeId ++) {
	  maskIdx = coordMap[theGraph.nodeFromId( nodeId )];
	  printf("nodeId = %i, Idx = [%i, %i, %i]\n", nodeId, maskIdx[0], maskIdx[1], maskIdx[2]);
     }

     lemon::SmartGraph::Edge edge;     
     for (int edgeId = 0; edgeId < 10; edgeId ++) {
	  edge = theGraph.edgeFromId(edgeId);
	  printf("edgeId = %i, btw node %i and %i\n", edgeId, theGraph.id(theGraph.u(edge)), theGraph.id(theGraph.v(edge)));
     }
     for (int edgeId = theGraph.maxEdgeId() - 10; edgeId <= theGraph.maxEdgeId(); edgeId ++) {
	  edge = theGraph.edgeFromId(edgeId);
	  printf("edgeId = %i, btw node %i and %i\n", edgeId, theGraph.id(theGraph.u(edge)), theGraph.id(theGraph.v(edge)));
     }


     return 0;
}


int BuildDataMap(lemon::SmartGraph & theGraph, 
		 lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coordMap,
		 lemon::SmartGraph::NodeMap<float > & dataMap,
		 std::string dataFile)
{
     ReaderType3DFloat::Pointer dataReader = ReaderType3DFloat::New();
     dataReader->SetFileName(dataFile);
     ImageType3DFloat::Pointer dataPtr = dataReader->GetOutput();
     dataPtr->Update();
     ImageType3DFloat::IndexType dataIdx;

     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
     	  dataIdx = coordMap[nodeIt];
	  dataMap[nodeIt] = dataPtr->GetPixel(coordMap[nodeIt]);
     }

     return 0;
}


int EstimateMuSigma(lemon::SmartGraph & theGraph, 
		    lemon::SmartGraph::NodeMap<float > & dataMap,
		    lemon::SmartGraph::NodeMap< std::vector< unsigned char > > & labelMap,
		    SWParType & par)
{

     unsigned mcIdx = 0; // Monte Carlo sample index.
     unsigned short clsIdx = 0;

     // reset all mu and sigma to zero.
     for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
     	  par.comp[clsIdx].mu = 0;
     	  par.comp[clsIdx].sigma = 0;
     	  par.comp[clsIdx].numPts = 0;
     }

     // estimate Pi_k
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
     	  for (mcIdx = 0; mcIdx < par.numSamples; mcIdx ++) {
     	       clsIdx = labelMap[nodeIt][mcIdx];
     	       par.comp[clsIdx].numPts ++;
     	  }
     }
     
     for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
     	  par.comp[clsIdx].prop = par.comp[clsIdx].numPts / (par.totalPts * par.numSamples);
     }

     // Estimate mu.
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
     	  for (mcIdx = 0; mcIdx < par.numSamples; mcIdx ++) {
     	       clsIdx = labelMap[nodeIt][mcIdx];
     	       par.comp[clsIdx].mu += dataMap[nodeIt];
     	  }
     }
     
     // normalize mu by dividing by K.
     for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
     	  par.comp[clsIdx].mu = par.comp[clsIdx].mu / par.comp[clsIdx].numPts;
     }

     // Estimate sigma
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
     	  for (mcIdx = 0; mcIdx < par.numSamples; mcIdx ++) {
     	       clsIdx = labelMap[nodeIt][mcIdx];
     	       par.comp[clsIdx].sigma += pow( (dataMap[nodeIt] - par.comp[clsIdx].mu), 2);
     	  }
     }

     // normalize sigma
     for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
     	  par.comp[clsIdx].sigma = sqrt(par.comp[clsIdx].sigma / par.comp[clsIdx].numPts);
     }

     return 0;
}

int Sampling(lemon::SmartGraph & theGraph, 
	     lemon::SmartGraph::NodeMap<float > & dataMap,
	     lemon::SmartGraph::NodeMap< std::vector< unsigned char > > & labelMap,
	     SWParType & par)
{
     // Uniform integer generator.
     boost::random::uniform_int_distribution<> uni_int(0, par.numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::random::uniform_int_distribution<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     unsigned scanIdx = 0;
     float oldPriorEngy = 0, newPriorEngy = 0, oldDataEngy = 0, newDataEngy = 0, denergy;
     unsigned short curLabel = 0, cand = 0;
     float p_acpt = 0;
     unsigned M = par.numSamples;

     for (scanIdx = 0 ;scanIdx < par.burnin + par.numSamples; scanIdx ++) {
	  for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {

	       // compute old prior energy.
	       oldPriorEngy = 0;

	       for (SmartGraph::IncEdgeIt edgeIt(theGraph, nodeIt); edgeIt != INVALID; ++ edgeIt) {
		    oldPriorEngy += Delta( labelMap[theGraph.u(edgeIt)][M-1], labelMap[theGraph.v(edgeIt)][M-1] );
	  }
	       oldPriorEngy = - par.beta * oldPriorEngy;

	       // compute candidate's prior energy.
	       cand = roll_die();
	       newPriorEngy = 0;
	       for (SmartGraph::IncEdgeIt edgeIt(theGraph, nodeIt); edgeIt != INVALID; ++ edgeIt) {
		    newPriorEngy += Delta( cand, labelMap[theGraph.runningNode(edgeIt)][M-1] );

	       }
	       newPriorEngy = - par.beta * newPriorEngy;
	       
	       // compute old data energy.
	       curLabel = labelMap[nodeIt][M-1];
	       oldDataEngy = log ( par.comp[curLabel].sigma) + pow( dataMap[nodeIt] - par.comp[curLabel].mu, 2) / (2 * pow(par.comp[curLabel].sigma, 2) );

	       // compute new data energy.
	       newDataEngy = log ( par.comp[cand].sigma) + pow( dataMap[nodeIt] - par.comp[cand].mu, 2) / (2 * pow(par.comp[cand].sigma, 2) );
	       
	       // change of energy.
	       denergy = (newPriorEngy + newDataEngy) - (oldPriorEngy + oldDataEngy);
	       denergy = denergy / par.temperature;
	       
	       if (denergy <= 0) {
		    labelMap[nodeIt][M-1] = cand;
	       }
	       else {
		    p_acpt = exp(-denergy);
		    if (uni() < p_acpt) {
			 labelMap[nodeIt][M-1] = cand;
		    }
	       }
	  } // NodeIt

	  // after burnin pariod, save it to correct place.
	  if (scanIdx >= par.burnin) {
	       for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
		    labelMap[nodeIt][scanIdx - par.burnin] = labelMap[nodeIt][M-1];
	       }	       
	  }
     }

     return 0;
}


unsigned short Delta(unsigned xi, unsigned xj)
{
     return ((xi == xj) ? 1:0);
}


int SWSampling(lemon::SmartGraph & theGraph, 
	       lemon::SmartGraph::NodeMap<float > & dataMap,
	       lemon::SmartGraph::NodeMap< std::vector< unsigned char > > & labelMap,
	       lemon::SmartGraph::EdgeMap<bool > & edgeMap,
	       SWParType & par) 
{
     unsigned numComp = 0;
     unsigned M = par.numSamples;

     // Define Bernoulli distribution.
     boost::bernoulli_distribution<> Ber_Dist(1 - exp(-par.beta));
     boost::variate_generator<twister_base_gen_type&, boost::bernoulli_distribution<> > Ber(mygenerator, Ber_Dist);

     // Uniform integer generator.
     boost::random::uniform_int_distribution<> uni_int(0, par.numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::random::uniform_int_distribution<> > roll_die(mygenerator, uni_int);

     // unfair die.
     boost::random::discrete_distribution<> disc_dist; 
     // boost::variate_generator<twister_base_gen_type&, boost::random::discreate_distribution<> > roll_die(mygenerator, uni_int);

     // Define component Map. The values is in [0, C-1], where C is the number
     // of connected components.
     lemon::SmartGraph::NodeMap<unsigned> compMap(theGraph);

     // new random labels of the conn. comp. 
     std::vector<unsigned short > newLabelSet;

     // Using edgeMap as filter map, get a subgraph with some edges removed.
     lemon::FilterEdges<SmartGraph> subGraph(theGraph, edgeMap);

     // init label map with random labels.
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  labelMap[nodeIt][M-1] = roll_die();
     }
     
     // probability of all lablels for one connected component.
     std::vector< std::vector<double> > labelWeightSet(par.numClusters);

     // sampling starts here.
     for (unsigned scanIdx = 0; scanIdx < par.burnin + par.numSamples; scanIdx ++) {
	  if (par.verbose >= 1 && (scanIdx%10 == 0) ) {
	       printf("SWSampling(), scan %i -->\n", scanIdx);
	  }
	  // sampling edge
	  for (SmartGraph::EdgeIt edgeIt(theGraph); edgeIt !=INVALID; ++ edgeIt) {
	       if (labelMap[theGraph.u(edgeIt)][M-1] ==  labelMap[theGraph.v(edgeIt)][M-1]) {
		    edgeMap[edgeIt] = Ber();
	       }
	       else {
		    edgeMap[edgeIt] = 0;
	       }
	  } // edgeIt

	  // Sampling labels.

	  // Put the labels of the subgraph into clusters (connected component)
	  numComp = connectedComponents(subGraph, compMap);
	  newLabelSet.resize(numComp);
	  labelWeightSet.assign( numComp, std::vector<double> (par.numClusters, 0) );

	  // test
	  for (unsigned compIdx = 0; compIdx < numComp; compIdx ++) {
	       labelWeightSet[compIdx][par.numClusters - 1] = 0;
	  }

	  // compute probability for each label at each cc 
	  ComputeExtField(labelWeightSet, compMap, theGraph, labelMap, dataMap, par);
	  
	  // Generate samples for all connected components.

	  for (unsigned compIdx = 0; compIdx < numComp; compIdx ++) {
	       boost::random::discrete_distribution<>::param_type myparam(labelWeightSet[compIdx].begin(), labelWeightSet[compIdx].end() );
	       newLabelSet[compIdx] = disc_dist(mygenerator, myparam);
	  }

	  // assign new labels to all connected components.
	  for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	       labelMap[nodeIt][M-1] = newLabelSet[compMap[nodeIt]];
	  }

	  // after burnin pariod, save it to correct place.
	  if (scanIdx >= par.burnin) {
	       for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
		    labelMap[nodeIt][scanIdx - par.burnin] = labelMap[nodeIt][M-1];
	       }	       
	  }

     } // scanIdx

     return 0;
}


int ComputeExtField(std::vector< std::vector<double> > & labelWeightSet,
		    lemon::SmartGraph::NodeMap<unsigned> & compMap,
		    lemon::SmartGraph & theGraph, 
		    lemon::SmartGraph::NodeMap< std::vector< unsigned char > > & labelMap,
		    lemon::SmartGraph::NodeMap<float > & dataMap,
		    SWParType & par) 
{
     unsigned clsIdx = 0;
     double alpha = 0;
     
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       alpha = - ( pow( (dataMap[nodeIt] - par.comp[clsIdx].mu), 2) ) / (2 * pow(par.comp[clsIdx].sigma, 2) ) - (log(par.comp[clsIdx].sigma) );
	       labelWeightSet[compMap[nodeIt]][clsIdx] += alpha;
	  }
     }

     // substract mean of each weight vector for numerical stability.
     double maxWeight = 0;
     for (unsigned ccIdx = 0; ccIdx < labelWeightSet.size(); ccIdx ++) {
	  maxWeight = -1000;
	  for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       if (labelWeightSet[ccIdx][clsIdx] > maxWeight) {
		    maxWeight = labelWeightSet[ccIdx][clsIdx];
	       }
	  }
	  for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       labelWeightSet[ccIdx][clsIdx] =  exp(labelWeightSet[ccIdx][clsIdx] - maxWeight);
	  }
     }

     return 0;
}
		    

int PrintParam(SWParType & par, unsigned plevel)
{
     if (plevel >= 0) {
	  printf("numClusters = %i,    numSamples = %i,    burnin = %i, initTemp = %f, finalTemp = %f, temperature = %f, \nbeta = %f,    gamma = %f,    verbose = %i,    totalPts = %i\n",par.numClusters, par.numSamples, par.burnin, par.initTemp, par.finalTemp, par.temperature, par.beta, par.gamma, par.verbose, par.totalPts);
     }
     if (plevel >= 1) {
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++ ) {
	       printf("cls[%i]: mu = %f, sigma = %f, numPts = %i, prop = %f\n", clsIdx, par.comp[clsIdx].mu, par.comp[clsIdx].sigma, par.comp[clsIdx].numPts, par.comp[clsIdx].prop);
	  }
     }
     
     return 0;
}



// int SaveGraphToImage(lemon::SmartGraph & theGraph, 
// 		     lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coordMap,
// 		     lemon::SmartGraph::NodeMap<std::vector<unsigned char> > & labelMap,
// 		     ImageType3DChar::Pointer maskPtr,
// 		     std::string outFile)
// {
//      unsigned M = labelMap[theGraph.nodeFromId( 0 )].size();

//      ImageType3DChar::IndexType outImageIdx;
//      outImageIdx.Fill(0);
//      ImageType3DChar::SizeType outImageSize = maskPtr->GetLargestPossibleRegion().GetSize();;

//      ImageType3DChar::RegionType outImageRegion;
//      outImageRegion.SetSize(outImageSize);
//      outImageRegion.SetIndex(outImageIdx);
//      ImageType3DChar::Pointer outImagePtr = ImageType3DChar::New();
//      outImagePtr->SetRegions(outImageRegion);
//      outImagePtr->Allocate();
//      outImagePtr->FillBuffer(0);
     
//      // write labels (plus 1) at each node to image.
//      for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
// 	  outImageIdx = coordMap[nodeIt];
// 	  outImagePtr->SetPixel(outImageIdx, labelMap[nodeIt][M-1]+1);
//      }

//      WriterType3DChar::Pointer writer = WriterType3DChar::New();
//      writer->SetInput( outImagePtr );
//      writer->SetFileName(outFile);

//      try 
//      { 
//      	  writer->Update(); 
//      } 
//      catch( itk::ExceptionObject & err ) 
//      { 
//      	  std::cerr << "ExceptionObject caught !" << std::endl; 
//      	  std::cerr << err << std::endl; 
//      	  return EXIT_FAILURE;
//      } 

//      std::cout << "SaveGraphToImage: File " << outFile << " saved.\n";

//      return 0;
// }
