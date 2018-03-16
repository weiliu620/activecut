#include <common.h>
#include <lemon/list_graph.h>

twister_base_gen_type mygenerator(42u);

unsigned ComputeNumSubs(std::string srcpath);

int BuildGraph(lemon::ListGraph & theGraph, 
	       lemon::ListGraph::NodeMap<SuperCoordType> & coordMap,
	       unsigned numSubs, 
	       ImageType3DChar::Pointer maskPtr);

int BuildDataMap(lemon::ListGraph & theGraph, 
		 lemon::ListGraph::NodeMap<SuperCoordType> & coordMap,
		 lemon::ListGraph::NodeMap<vnl_vector<float> > & tsMap,
		 std::string srcpath);

using namespace lemon;
int main(int argc, char* argv[])
{
     ParStruct par;
     unsigned seed = 0;
     unsigned int EMIter = 1;

     std::string fmriPath, initGrplabel, outGrpLabel, outGrpProb, outSubBase, 
	  grpSampleName, samplePrefix, initsubpath, groupprob, subbaseprob;

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
	   "initial number of Monte Carlo samples. ")
	  ("emiter", po::value<unsigned int>(&EMIter)->default_value(30),
	   "initial number of Monte Carlo samples. ")

	  ("betag", po::value<float>(&par.betag)->default_value(0.5),
	   "group level pairwise interation term")
	  ("betaz", po::value<float>(&par.betaz)->default_value(0.8),
	   "subject level pairwise interation term")
	  ("alpha", po::value<float>(&par.alpha)->default_value(0.7),
	   "connection between group lebel and individual subject level.")

	  ("alpha0", po::value<float>(&par.alpha0)->default_value(0.7),
	   "Mean of prior on alpha.")

	  ("sigma2", po::value<float>(&par.sigma2)->default_value(1),
	   "variance of prior on alpha")

	  ("gamma", po::value<float>(&par.gamma)->default_value(1),
	   " weights of lilelihood term.")

	  ("estprior", po::bool_switch(&estprior), "whether to estimate alpha and beta, the prior parameter.")
	  ("initsame", po::bool_switch(&initsame), "whether to init subject same to group label. choose yes to init same, otherwise init with subject label map.")

	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("numClusters,k", po::value<unsigned>(&par.numClusters)->default_value(6),
	   "Number of clusters. Default is 6.")

	   ("initgrouplabel,i", po::value<std::string>(&initGrplabel)->default_value("grouplabel.nii"), "Initial group level label map. Also used as mask file.")

	   ("initsubpath,t", po::value<std::string>(&initsubpath)->default_value("subpath"), "Initial subject label maps path")

	  ("fmripath,f", po::value<std::string>(&fmriPath)->default_value("."), "noised image file")

	  ("sampleprefix", po::value<std::string>(&samplePrefix)->default_value("subject"), "Monte Carlo samples file name prefix")

	   ("subbasename", po::value<std::string>(&outSubBase)->default_value("subject"), "Output individual subject's base name (with path).")

	  ("grouplabel,g", po::value<std::string>(&outGrpLabel)->default_value("outgrplabel.nii"), "output group label file. nii format.")

	  ("grpsamplename", po::value<std::string>(&grpSampleName)->default_value("grpSample.nii.gz"), "Monte Carlo group label samples file name.")


	  ("verbose,v", po::value<unsigned short>(&par.verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: hiermarf [options]\n";
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

     // original 1-based group map only used for mask.
     ImageType3DChar::Pointer maskPtr = grpReader->GetOutput();

     // get number of subjects.
     par.numSubs = ComputeNumSubs(fmriPath);

     printf("number of subjects: %i.\n", par.numSubs);

     // Construct graph.
     lemon::ListGraph theGraph;
     lemon::ListGraph::NodeMap<SuperCoordType> coordMap(theGraph); // node --> coordinates.
     BuildGraph(theGraph, coordMap, par.numSubs, maskPtr);

     // Build Data map.
     lemon::ListGraph::NodeMap<vnl_vector<float> > tsMap(theGraph); // node --> time series.
     BuildDataMap(theGraph, coordMap, tsMap, fmriPath);

     // define label map.
     lemon::ListGraph::NodeMap<unsigned char> labelMap(theGraph);
     
     // Define edge map.
     lemon::ListGraph::NodeMap<bool> edgeMap(theGraph);

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


int BuildGraph(lemon::ListGraph & theGraph, 
	       lemon::ListGraph::NodeMap<SuperCoordType> & coordMap,
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

     
     // Add nodes.

     lemon::ListGraph::Node curNode, neiNode;

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

     // Fill in group level nodes.
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
     	  if (maskIt.Get() > 0) {
     	       curNode = theGraph.addNode();
	       coordMap[curNode].idx = maskIt.GetIndex();
	       coordMap[curNode].subid = numSubs;
	       
	       nodeMapIdx[0] = coordMap[curNode].idx[0];
	       nodeMapIdx[1] = coordMap[curNode].idx[1];
	       nodeMapIdx[2] = coordMap[curNode].idx[2];
	       nodeMapIdx[3] = numSubs;
	       nodeMapPtr->SetPixel(nodeMapIdx, theGraph.id(curNode) );
	  }
     }

     printf("BuildGraph(): number of Nodes: %i\n", theGraph.maxNodeId()+1 );

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
     ImageType3DChar::IndexType maskIdx;
     int curNodeId = 0, neiNodeId = 0;
     // std::array<short, 6>::const_iterator neiIdxIt;
     auto neiIdxIt = neiIdxSet.begin();

     // edge within grp.
     for (maskNeiIt.GoToBegin(); !maskNeiIt.IsAtEnd(); ++ maskNeiIt) {
	  if (maskNeiIt.GetCenterPixel() > 0) {
	       maskIdx = maskNeiIt.GetIndex();
	       nodeMapIdx[0] = maskIdx[0];
	       nodeMapIdx[1] = maskIdx[1];
	       nodeMapIdx[2] = maskIdx[2];
	       nodeMapIdx[3] = numSubs;
	       curNodeId = nodeMapPtr->GetPixel(nodeMapIdx);
	       // curNode is group node.
	       curNode = theGraph.nodeFromId( curNodeId );

	       for (neiIdxIt=neiIdxSet.begin(); neiIdxIt < neiIdxSet.end(); neiIdxIt++) {
		    if (maskNeiIt.GetPixel(*neiIdxIt) > 0) {
			 neiIdx = maskNeiIt.GetIndex(*neiIdxIt);
			 nodeMapIdx[0] = neiIdx[0];
			 nodeMapIdx[1] = neiIdx[1];
			 nodeMapIdx[2] = neiIdx[2];
			 nodeMapIdx[3] = numSubs;
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
     printf("BuildGraph(): After adding edges within group, number of edges: %i\n", theGraph.maxEdgeId() + 1);

     // edge between grp and subjects.
     for (maskNeiIt.GoToBegin(); !maskNeiIt.IsAtEnd(); ++ maskNeiIt) {
	  if (maskNeiIt.GetCenterPixel() > 0) {
	       maskIdx = maskNeiIt.GetIndex();
	       nodeMapIdx[0] = maskIdx[0];
	       nodeMapIdx[1] = maskIdx[1];
	       nodeMapIdx[2] = maskIdx[2];

	       // get current node in group.
	       nodeMapIdx[3] = numSubs;
	       curNodeId = nodeMapPtr->GetPixel(nodeMapIdx);
	       curNode = theGraph.nodeFromId( curNodeId );
	       
	       // get subject node and connect them
	       for (unsigned short subIdx = 0; subIdx < numSubs; subIdx ++) {
		    nodeMapIdx[3] = subIdx;
		    neiNodeId = nodeMapPtr->GetPixel(nodeMapIdx);
		    neiNode = theGraph.nodeFromId( neiNodeId );
		    theGraph.addEdge(curNode, neiNode);
	       } // subIdx
	  } // maskNeiIt > 0
     } // MaskNeiIt

     printf("BuildGraph(): After adding edges btw grp and subs, number of edges: %i\n", theGraph.maxEdgeId() + 1);

     // edge within each subject.
     for (unsigned short subIdx = 0; subIdx < numSubs; subIdx ++) {
	  nodeMapIdx[3] = subIdx;

	  for (maskNeiIt.GoToBegin(); !maskNeiIt.IsAtEnd(); ++ maskNeiIt) {
	       if (maskNeiIt.GetCenterPixel() > 0) {
		    maskIdx = maskNeiIt.GetIndex();
		    nodeMapIdx[0] = maskIdx[0];
		    nodeMapIdx[1] = maskIdx[1];
		    nodeMapIdx[2] = maskIdx[2];
		    curNodeId = nodeMapPtr->GetPixel(nodeMapIdx);
		    // curNode is in subject map.
		    curNode = theGraph.nodeFromId( curNodeId );

		    for (neiIdxIt=neiIdxSet.begin(); neiIdxIt < neiIdxSet.end(); neiIdxIt++) {
			 if (maskNeiIt.GetPixel(*neiIdxIt) > 0) {
			      neiIdx = maskNeiIt.GetIndex(*neiIdxIt);
			      nodeMapIdx[0] = neiIdx[0];
			      nodeMapIdx[1] = neiIdx[1];
			      nodeMapIdx[2] = neiIdx[2];
			      nodeMapIdx[3] = numSubs;
			      neiNodeId = nodeMapPtr->GetPixel(nodeMapIdx);
			      // make sure each edge is added only once.
			      if (neiNodeId > curNodeId) {
				   neiNode = theGraph.nodeFromId( neiNodeId );
				   theGraph.addEdge(curNode, neiNode);
			      } // if
			 } // maskNeiIt
		    } // neiIdxIt
	       } // maskNeiIt >0 
	  }
     } // subIdx

     printf("BuildGraph(): number of edges: %i\n", theGraph.maxEdgeId() + 1);

     // test code. Test nodeMap and coordMap.
     // SuperCoordType superCoord;
     // for (int nodeId = 0; nodeId < 10; nodeId ++) {
     // 	  superCoord = coordMap[theGraph.nodeFromId( nodeId )];
     // 	  printf("nodeId = %i, superCoord.idx = [%i, %i, %i], superCoord.subid = %i\n", nodeId, superCoord.idx[0], superCoord.idx[1], superCoord.idx[2], superCoord.subid);
     // }

     // for (int nodeId = theGraph.maxNodeId() - 10; nodeId <= theGraph.maxNodeId(); nodeId ++) {
     // 	  superCoord = coordMap[theGraph.nodeFromId( nodeId )];
     // 	  printf("nodeId = %i, superCoord.idx = [%i, %i, %i], superCoord.subid = %i\n", nodeId, superCoord.idx[0], superCoord.idx[1], superCoord.idx[2], superCoord.subid);
     // }

     // lemon::ListGraph::Edge edge;     
     // for (int edgeId = 0; edgeId < 10; edgeId ++) {
     // 	  edge = theGraph.edgeFromId(edgeId);
     // 	  printf("edgeId = %i, btw node %i and %i\n", edgeId, theGraph.id(theGraph.u(edge)), theGraph.id(theGraph.v(edge)));
     // }
     // for (int edgeId = theGraph.maxEdgeId() - 10; edgeId <= theGraph.maxEdgeId(); edgeId ++) {
     // 	  edge = theGraph.edgeFromId(edgeId);
     // 	  printf("edgeId = %i, btw node %i and %i\n", edgeId, theGraph.id(theGraph.u(edge)), theGraph.id(theGraph.v(edge)));
     // }
     return 0;
}


int BuildDataMap(lemon::ListGraph & theGraph, 
		 lemon::ListGraph::NodeMap<SuperCoordType> & coordMap,
		 lemon::ListGraph::NodeMap<vnl_vector<float> > & tsMap,
		 std::string srcpath)
{
     // some of the code below refers to boost filesystem tutorial tut04.cpp.
     // prepare for sorting directory entries.
     typedef std::vector<boost::filesystem::path> PathVec;
     PathVec  sortedEntries;
     copy(boost::filesystem::directory_iterator(srcpath), boost::filesystem::directory_iterator(), std::back_inserter(sortedEntries) );
     sort(sortedEntries.begin(), sortedEntries.end() );

     // Get number of subjects.
     unsigned numSubs  = sortedEntries.size();

     // init file reader pointers.
     std::vector<ReaderType4DFloat::Pointer> readerVec(numSubs);
     std::vector<ImageType4DF::Pointer> srcPtrVec(numSubs);
     for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	  readerVec[subIdx] = ReaderType4DFloat::New();
	  readerVec[subIdx]->SetFileName( sortedEntries[subIdx].string() );
	  readerVec[subIdx]->Update();
	  srcPtrVec[subIdx]= readerVec[subIdx]->GetOutput();
     }

     ImageType4DF::SizeType srcSize = readerVec[0]->GetOutput()->GetLargestPossibleRegion().GetSize();
     unsigned tsLength = srcSize[3];

     printf("BuildDataMap(): Time series length: %i.\n", tsLength);

     SuperCoordType superCoord;
     ImageType4DF::IndexType srcIdx;
     for (ListGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
     	  superCoord = coordMap[nodeIt];
	  // when not in group map.
	  if (superCoord.subid < numSubs) {
	       tsMap[nodeIt].set_size(tsLength);

	       srcIdx[0] = superCoord.idx[0];
	       srcIdx[1] = superCoord.idx[1];
	       srcIdx[2] = superCoord.idx[2];
	       for (srcIdx[3] = 0; srcIdx[3] < tsLength; srcIdx[3] ++) {
		    tsMap[nodeIt][srcIdx[3]] = srcPtrVec[superCoord.subid]->GetPixel(srcIdx);
	       }
	  }
     }

     return 0;
}
