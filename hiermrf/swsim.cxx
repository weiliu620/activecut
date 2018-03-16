#include <common.h>
#include <utility.h>
twister_base_gen_type mygenerator(42u);

int BuildGraph(lemon::SmartGraph & theGraph, 
	       lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coordMap,
	       ImageType3DChar::Pointer maskPtr);

int Sampling(lemon::SmartGraph & theGraph, 
	     lemon::SmartGraph::NodeMap<std::vector<unsigned char> > & labelMap,
	     SWParType & par);

int SWSampling(lemon::SmartGraph & theGraph, 
	       lemon::SmartGraph::NodeMap<std::vector<unsigned char> > & labelMap,
	       lemon::SmartGraph::EdgeMap<bool > & edgeMap,
	       SWParType & par);

int SWSamplingMaster(lemon::SmartGraph & theGraph, 
		     lemon::SmartGraph::NodeMap< unsigned char > & labelMap,
		     lemon::SmartGraph::EdgeMap<bool > & edgeMap,
		     SWParType & par);

int SaveGraphToImage(lemon::SmartGraph & theGraph, 
		     lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coordMap,
		     lemon::SmartGraph::NodeMap<std::vector<unsigned char> > & labelMap,
		     ImageType3DChar::Pointer maskPtr,
		     std::string outFile);

int SaveSamplesToImage(lemon::SmartGraph & theGraph, 
		       lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coordMap,
		       lemon::SmartGraph::NodeMap<std::vector<unsigned char> > & labelMap,
		       ImageType3DChar::Pointer maskPtr,
		       std::string outFile);

unsigned short Delta(unsigned xi, unsigned xj);

void * SWSamplingSlave(void * argList);

using namespace lemon;
int main(int argc, char* argv[])
{
     SWParType par;
     unsigned seed = 0;
     std::string  maskFile, oLabelFile, sampleFile;

     bool usesw = false;

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Simulate Potts model by Swendsen-Wang algorithm.")
	  ("burnin,b", po::value<unsigned>(&par.burnin)->default_value(20),
	   "number of scans for burn-in period. ")
	  ("numSamples,n", po::value<unsigned int>(&par.numSamples)->default_value(20),
	   "Number of Monte Carlo samples saved at the end of burnin. ")
	  ("inittemp", po::value<double>(&par.initTemp)->default_value(1),
	   "Initial temperature for annealing.")
	  ("finaltemp", po::value<double>(&par.finalTemp)->default_value(1),
	   "Final temperature for annealing to stop")

	  ("beta", po::value<double>(&par.beta)->default_value(0.5),
	   "pairwise interation term")

	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("numClusters,k", po::value<unsigned>(&par.numClusters)->default_value(6),
	   "Number of clusters. Default is 6.")

	   ("maskfile,m", po::value<std::string>(&maskFile)->default_value("maskFile.nii.gz"), "mask file.")

	  ("outlabel,g", po::value<std::string>(&oLabelFile)->default_value("oLabelFile.nii.gz"), "output label file. must be nii or nii.gz.")

	  ("samplefile,a", po::value<std::string>(&sampleFile)->default_value("samples.nii.gz"), "output sample file. must be nii or nii.gz.")

	  ("usesw,s", po::bool_switch(&usesw), "whether to use Swendsen-Wang for sampling.")

	  ("verbose,v", po::value<unsigned short>(&par.verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: swsim [options]\n";
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

     // read in mask file.
     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(maskFile);
     maskReader->Update();
     ImageType3DChar::Pointer maskPtr = maskReader->GetOutput();

     // Construct graph.
     lemon::SmartGraph theGraph;
     lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> coordMap(theGraph); // node --> coordinates.
     BuildGraph(theGraph, coordMap, maskPtr);

     // define label map.
     lemon::SmartGraph::NodeMap<std::vector<unsigned char> > labelMap(theGraph);
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  labelMap[nodeIt].resize(par.numSamples);
     }

     // Define edge map.
     lemon::SmartGraph::EdgeMap<bool > edgeMap(theGraph);

     par.temperature = 1;


     // sampling.
     if (usesw) {
	  SWSampling(theGraph, labelMap, edgeMap, par);
     }
     else {
	  Sampling(theGraph, labelMap, par);
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
     unsigned int neiIdxSet[] = {14, 12, 16, 10, 22, 4}; 

     int curNodeId = 0, neiNodeId = 0;

     for (maskNeiIt.GoToBegin(); !maskNeiIt.IsAtEnd(); ++ maskNeiIt) {
	  if (maskNeiIt.GetCenterPixel() > 0) {
	       maskIdx = maskNeiIt.GetIndex();
	       nodeMapIdx = maskIdx;
	       curNodeId = nodeMapPtr->GetPixel(nodeMapIdx);

	       curNode = theGraph.nodeFromId( curNodeId );

	       for (unsigned nid = 0; nid < 6; nid++) {
		    
		    if (maskNeiIt.GetPixel(neiIdxSet[nid]) > 0) {
			 neiIdx = maskNeiIt.GetIndex(neiIdxSet[nid]);
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
     
     // for (int nodeId = 0; nodeId < 10; nodeId ++) {
     // 	  maskIdx = coordMap[theGraph.nodeFromId( nodeId )];
     // 	  // printf("nodeId = %i, Idx = [%i, %i, %i]\n", nodeId, maskIdx[0], maskIdx[1], maskIdx[2]);
     // }

     // for (int nodeId = theGraph.maxNodeId() - 10; nodeId <= theGraph.maxNodeId(); nodeId ++) {
     // 	  maskIdx = coordMap[theGraph.nodeFromId( nodeId )];
     // 	  // printf("nodeId = %i, Idx = [%i, %i, %i]\n", nodeId, maskIdx[0], maskIdx[1], maskIdx[2]);
     // }

     // lemon::SmartGraph::Edge edge;     
     // for (int edgeId = 0; edgeId < 10; edgeId ++) {
     // 	  edge = theGraph.edgeFromId(edgeId);
     // 	  // printf("edgeId = %i, btw node %i and %i\n", edgeId, theGraph.id(theGraph.u(edge)), theGraph.id(theGraph.v(edge)));
     // }
     // for (int edgeId = theGraph.maxEdgeId() - 10; edgeId <= theGraph.maxEdgeId(); edgeId ++) {
     // 	  edge = theGraph.edgeFromId(edgeId);
     // 	  // printf("edgeId = %i, btw node %i and %i\n", edgeId, theGraph.id(theGraph.u(edge)), theGraph.id(theGraph.v(edge)));
     // }
     return 0;
}

int Sampling(lemon::SmartGraph & theGraph, 
	     lemon::SmartGraph::NodeMap<std::vector<unsigned char> > & labelMap,
	     SWParType & par)
{
     unsigned M = par.numSamples;

     // Uniform integer generator.
     boost::random::uniform_int_distribution<> uni_int(0, par.numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::random::uniform_int_distribution<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::random::uniform_real_distribution<> uni_real(0, 1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::random::uniform_real_distribution<> > uni(mygenerator, uni_real);

     unsigned scanIdx = 0;
     float oldPriorEngy = 0, newPriorEngy = 0, denergy;
     unsigned short cand = 0;
     float p_acpt = 0;

     // init label map with random labels.
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  labelMap[nodeIt][M-1] = roll_die();
	  // labelMap[nodeIt][M-1] = 0;
     }

     // sampling starts here.
     printf("par.burnin = %i", par.burnin);
     for (scanIdx = 0 ;scanIdx < par.burnin + par.numSamples; scanIdx ++) {
	  if (par.verbose >= 1 && (scanIdx%10 == 0) ) {
	       printf("Sampling(), scan %i -->\n", scanIdx);
	  }
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
		    if (theGraph.u(edgeIt) == nodeIt) {
			 // v is neighbors.
			 newPriorEngy += Delta( cand, labelMap[theGraph.v(edgeIt)][M-1] );
		    }
		    else {
			 // u is neighbors.
			 newPriorEngy += Delta( cand, labelMap[theGraph.u(edgeIt)][M-1] );
		    }
	       }
	       newPriorEngy = - par.beta * newPriorEngy;
	       
	       // change of energy.
	       denergy = newPriorEngy - oldPriorEngy;
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
	       lemon::SmartGraph::NodeMap<std::vector<unsigned char> > & labelMap,
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

     // Define component Map. The values is in [0, C-1], where C is the number
     // of connected components.
     lemon::SmartGraph::NodeMap<unsigned> compMap(theGraph);

     // new random labels of the conn. comp. 
     std::vector<unsigned short> newLabelSet;

     // Using edgeMap as filter map, get a subgraph with some edges removed.
     lemon::FilterEdges<SmartGraph> subGraph(theGraph, edgeMap);

     // init label map with random labels.
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  // labelMap[nodeIt][M-1] = roll_die();
	  labelMap[nodeIt][M-1] = 0;
     }

     // sampling starts here.
     for (unsigned scanIdx = 0 ;scanIdx < par.burnin + par.numSamples; scanIdx ++) {
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

	  // Generate samples for all connected components.
	  for (unsigned compIdx = 0; compIdx < numComp; compIdx ++) {
	       newLabelSet[compIdx] = roll_die();
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



// int SWSamplingMaster(lemon::SmartGraph & theGraph, 
// 			lemon::SmartGraph::NodeMap< unsigned char > & labelMap,
// 			lemon::SmartGraph::EdgeMap<bool > & edgeMap,
// 			SWParType & par)
// {
     
//      unsigned numCPU = sysconf( _SC_NPROCESSORS_ONLN );
//      unsigned numThreads = int(numCPU/2);
//      std::vector<pthread_t> threadIDSet( numThreads );
//      // std::vector<pthread_t> threadReSet( int(numCPU/2));
     
//      ArgList argList;
//      argList.par = par;
//      argList.numThreads = numThreads;
//      argList.numPtsPerThreads = ceil(double(lemon::countNodes(theGraph)) / double(numThreads));
//      argList.theGraph = & theGraph;
//      argList.labelMap = & labelMap;
//      argList.edgeMap = & edgeMap;
	  
//      for (unsigned threadIdx = 0; threadIdx < threadIDSet.size(); threadIdx ++) {
// 	  pthread_create(& threadIDSet[threadIdx], NULL, SWSamplingSingleThread, (void*) & argList);
//      }
// };

// void * SWSamplingSlave(void * argList)
// {
//      unsigned K = ((ArgList *) argList)->numPtsPerThreads;
// }
