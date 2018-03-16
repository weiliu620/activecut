double logBesselI(float nu, double x);

int SaveCumuSample(lemon::SmartGraph & theGraph, 
		   lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coord,
		    lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
		    ParStruct & par,
		   std::string cumuSampleFile);

int PrintPar(unsigned prlevel, ParStruct & par);


