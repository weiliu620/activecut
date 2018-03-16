#include <common.h>
using namespace lemon;
double logBesselI(float nu, double x)
{
     // From
     // http://people.kyb.tuebingen.mpg.de/suvrit/work/progs/movmf/. See
     // matlab code logbesseli.m, and
     // http://people.math.sfu.ca/~cbm/aands/page_378.htm (eq 9.7.7)
     double const Pi = 4 * atan(1);

     double frac = x/nu;
     double square = 1 + frac * frac;
     double root = sqrt(square);
     double eta = root + log(frac) - log(1+root);
     double logb = - log(sqrt(2 * Pi * nu)) + nu * eta - 0.25 * log(square);
     return (logb);
}


int SaveCumuSample(lemon::SmartGraph & theGraph, 
		   lemon::SmartGraph::NodeMap<ImageType3DChar::IndexType> & coord,
		   lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
		   ParStruct & par,
		   std::string cumuSampleFile)
{
     ImageType3DChar::SizeType maskSize = par.maskSize;
     ImageType4DS::IndexType outImageIdx;
     outImageIdx.Fill(0);
     ImageType4DS::SizeType outImageSize;
     outImageSize[0] = maskSize[0];
     outImageSize[1] = maskSize[1];
     outImageSize[2] = maskSize[2];
     outImageSize[3] = par.numClusters;

     ImageType4DS::RegionType outImageRegion;
     outImageRegion.SetSize(outImageSize);
     outImageRegion.SetIndex(outImageIdx);
     ImageType4DS::Pointer outImagePtr = ImageType4DS::New();
     outImagePtr->SetRegions(outImageRegion);
     outImagePtr->Allocate();
     outImagePtr->FillBuffer(0);

     // Save group sample file.
     for (lemon::SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=lemon::INVALID; ++ nodeIt) {

	  outImageIdx[0] = coord[nodeIt][0];
	  outImageIdx[1] = coord[nodeIt][1];
	  outImageIdx[2] = coord[nodeIt][2];
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       outImageIdx[3] = clsIdx;
	       outImagePtr->SetPixel(outImageIdx, cumuSampleMap[nodeIt][clsIdx]);
	  }
     } // nodeIt

     WriterType4DS::Pointer writer = WriterType4DS::New();
     writer->SetInput( outImagePtr );
     writer->SetFileName(cumuSampleFile);

     try 
     { 
	  writer->Update(); 
     } 
     catch( itk::ExceptionObject & err ) 
     { 
	  std::cerr << "ExceptionObject caught !" << std::endl; 
	  std::cerr << err << std::endl; 
	  return EXIT_FAILURE;
     }      
     std::cout << "SaveGrpCumuSample(): File  " << cumuSampleFile << " saved.\n";

     return 0;
}

int PrintPar(unsigned prlevel, ParStruct & par)
{
     printf(" numClusters: %d\n numSamples: %d\n burnin: %d\n initTemp: %4.3f\n finalTemp: %4.3f\n temp: %4.3f\n tsLength:%d\n", par.numClusters, par.numSamples, par.burnin, par.initTemp, par.finalTemp, par.temperature, par.tsLength);

     if (prlevel >= 1) {
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       printf("cls[%i]: mu=[%2.2f %2.2f], kappa=%3.1f, numPts=%ld\n", clsIdx+1, par.vmm.comp[clsIdx].mu[0], par.vmm.comp[clsIdx].mu[1], par.vmm.comp[clsIdx].kappa, par.vmm.comp[clsIdx].numPts);
	  }
     }
     return 0;
}

