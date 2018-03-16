#include <common.h>
#include <lemon/list_graph.h>
#include <lemon/adaptors.h>
#include <lemon/connectivity.h>

twister_base_gen_type mygenerator(42u);
using namespace lemon;

struct Derivatives {
     double first;
     double second;
};

int BuildGraph(lemon::SmartGraph & theGraph, 
	       lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coordMap,
	       ImageType3DChar::Pointer maskPtr);

int EstimateBeta(lemon::SmartGraph & theGraph,
		 lemon::SmartGraph::NodeMap<unsigned char > & labelMap,
		 SWParType & par);

int ComputeDeriv(lemon::SmartGraph & theGraph,
		 lemon::SmartGraph::NodeMap<unsigned char > & labelMap,
		 SWParType & par,
		 Derivatives & deriv);

unsigned short Delta(unsigned xi, unsigned xj);
     
int main(int argc, char* argv[])
{
     SWParType par;
     std::string  labelFile;

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Inference on group mrf.")

	  ("numClusters,k", po::value<unsigned>(&par.numClusters)->default_value(6),
	   "Number of clusters. Default is 6.")

	  ("beta", po::value<double>(&par.beta)->default_value(0.5),
	   "pairwise interation term")

	   ("labelfile,i", po::value<std::string>(&labelFile)->default_value("labelFile.nii.gz"), "label file.")

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

     // read in label file and minus 1.
     ReaderType3DChar::Pointer labelReader = ReaderType3DChar::New();
     labelReader->SetFileName(labelFile);

     AddConstantToImageFilterType::Pointer myfilter = AddConstantToImageFilterType::New();
     myfilter->SetInput( labelReader->GetOutput());
     myfilter->SetConstant(-1);
     myfilter->InPlaceOff();
     myfilter->Update();
     ImageType3DChar::Pointer labelPtr = myfilter->GetOutput();

     // original 1-based group map only used for mask.
     ImageType3DChar::Pointer maskPtr = labelReader->GetOutput();


     // Construct graph.
     lemon::SmartGraph theGraph;
     lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> coordMap(theGraph); // node --> coordinates.
     BuildGraph(theGraph, coordMap, maskPtr);

     // Build label map.
     lemon::SmartGraph::NodeMap<unsigned char > labelMap(theGraph);
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
  	  labelMap[nodeIt] = labelPtr->GetPixel(coordMap[nodeIt]);
     }

     // Estimate Beta.
     EstimateBeta(theGraph, labelMap, par);
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

     for (maskNeiIt.GoToBegin(); !maskNeiIt.IsAtEnd(); ++ maskNeiIt) {
	  if (maskNeiIt.GetCenterPixel() > 0) {
	       maskIdx = maskNeiIt.GetIndex();
	       nodeMapIdx = maskIdx;
	       curNodeId = nodeMapPtr->GetPixel(nodeMapIdx);

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

     return 0;
}


int EstimateBeta(lemon::SmartGraph & theGraph,
		 lemon::SmartGraph::NodeMap<unsigned char > & labelMap,
		 SWParType & par)
{
     Derivatives deriv;
     do {
	  ComputeDeriv(theGraph, labelMap, par, deriv);
	  printf("1st deriv = %f, 2nd deriv = %f, beta = %f\n", deriv.first, deriv.second, par.beta);
	  par.beta = par.beta - deriv.first/deriv.second;
     } while (fabs(deriv.first) > 1e-5);
     return 0;
}

int ComputeDeriv(lemon::SmartGraph & theGraph,
		 lemon::SmartGraph::NodeMap<unsigned char > & labelMap,
		 SWParType & par,
		 Derivatives & deriv)
{
     double Si = 0, Sk = 0;
     double M0 = 0, M1 = 0, M2 = 0;
     unsigned short clsIdx = 0;
     lemon::SmartGraph::Node neiNode;
     deriv.first = 0;
     deriv.second = 0;
     for (SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=INVALID; ++ nodeIt) {
	  // compute Si
	  Si = 0;
	  for (SmartGraph::IncEdgeIt edgeIt(theGraph, nodeIt); edgeIt != INVALID; ++ edgeIt) {
	       Si += Delta( labelMap[theGraph.u(edgeIt)], labelMap[theGraph.v(edgeIt)] );
	  }

	  // compute M0, M1, M2
	  M0 = 0, M1 = 0, M2 = 0;
	  for (clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       Sk = 0;
	       for (SmartGraph::IncEdgeIt edgeIt(theGraph, nodeIt); edgeIt != INVALID; ++ edgeIt) {
		    // tell which one of u or v is neighbors.
		    if (theGraph.u(edgeIt) == nodeIt) {
			 neiNode = theGraph.v(edgeIt);
		    }
		    else {
			 neiNode = theGraph.u(edgeIt);
		    }
		    Sk += Delta( clsIdx, labelMap[neiNode] );		    
	       } // edgeIt
	       M0 += exp(par.beta * Sk);
	       M1 += exp(par.beta * Sk) * Sk;
	       M2 += exp(par.beta * Sk) * Sk * Sk;
	  } // clsIdx
	  deriv.first += (Si - M1/M0);
	  deriv.second += (M1*M1 - M0*M2) / (M0*M0);
     } // nodeIt
     
     return 0;
}



unsigned short Delta(unsigned xi, unsigned xj)
{
     return ((xi == xj) ? 1:0);
}
