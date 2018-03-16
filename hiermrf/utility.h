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
