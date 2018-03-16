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


int SaveCumuSamples(lemon::SmartGraph & theGraph, 
		    lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		    lemon::SmartGraph::NodeMap<std::vector<unsigned short> > & cumuSampleMap,
		    ParStruct & par)
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
     
     for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  std::string thisSubFilename = par.cumusampledir + "/" + par.vmm.sub[subIdx].name + ".nii.gz";
	  outImagePtr->FillBuffer(0);

	  for (lemon::SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=lemon::INVALID; ++ nodeIt) {
	       // if this is subject node.
	       if (coordMap[nodeIt].subid == subIdx) {
		    outImageIdx[0] = coordMap[nodeIt].idx[0];
		    outImageIdx[1] = coordMap[nodeIt].idx[1];
		    outImageIdx[2] = coordMap[nodeIt].idx[2];
		    for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
			 outImageIdx[3] = clsIdx;
			 outImagePtr->SetPixel(outImageIdx, cumuSampleMap[nodeIt][clsIdx]);
		    }
	       }
	  } // nodeIt

	  WriterType4DS::Pointer writer = WriterType4DS::New();
	  writer->SetInput( outImagePtr );
	  writer->SetFileName(thisSubFilename);

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
     } // subIdx
     std::cout << "SaveSubCumuSamples(): Subject samples saved into folder " << par.cumusampledir << std::endl;

     // Save group sample file.
     outImagePtr->FillBuffer(0);
     for (lemon::SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=lemon::INVALID; ++ nodeIt) {
	  // if this is group node.
	  if (coordMap[nodeIt].subid == par.numSubs) {
	       outImageIdx[0] = coordMap[nodeIt].idx[0];
	       outImageIdx[1] = coordMap[nodeIt].idx[1];
	       outImageIdx[2] = coordMap[nodeIt].idx[2];
	       for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
		    outImageIdx[3] = clsIdx;
		    outImagePtr->SetPixel(outImageIdx, cumuSampleMap[nodeIt][clsIdx]);
	       }
	  }
     } // nodeIt

     WriterType4DS::Pointer writer = WriterType4DS::New();
     writer->SetInput( outImagePtr );
     std::string outGrpFile = par.cumusampledir + "/" + "grp.nii.gz";
     writer->SetFileName(outGrpFile);

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
     std::cout << "SaveGrpCumuSamples(): File  " << outGrpFile << " saved.\n";

     return 0;
}

int SaveRunningSamples(lemon::SmartGraph & theGraph, 
		       lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
		       lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
		       ParStruct & par)

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

     std::string thisSubFilename;
     for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  thisSubFilename = par.rsampledir + "/" + par.vmm.sub[subIdx].name + ".nii.gz";
	  outImagePtr->FillBuffer(0);
	  for (lemon::SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=lemon::INVALID; ++ nodeIt) {
	       // if this is subject node.
	       if (coordMap[nodeIt].subid == subIdx) {
		    outImageIdx[0] = coordMap[nodeIt].idx[0];
		    outImageIdx[1] = coordMap[nodeIt].idx[1];
		    outImageIdx[2] = coordMap[nodeIt].idx[2];
		    outImageIdx[3] = rSampleMap[nodeIt].find_first();
		    outImagePtr->SetPixel(outImageIdx, 1);
	       }
	  } // nodeIt

	  WriterType4DS::Pointer writer = WriterType4DS::New();
	  writer->SetInput( outImagePtr );
	  writer->SetFileName(thisSubFilename);
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
     } // subIdx

     std::cout << "SaveRunningSample(): Subject samples saved into folder " << par.rsampledir << std::endl;

     // Now it's for group 
     outImagePtr->FillBuffer(0);
     for (lemon::SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=lemon::INVALID; ++ nodeIt) {
	  // if this is group node.
	  if (coordMap[nodeIt].subid == par.numSubs) {
	       outImageIdx[0] = coordMap[nodeIt].idx[0];
	       outImageIdx[1] = coordMap[nodeIt].idx[1];
	       outImageIdx[2] = coordMap[nodeIt].idx[2];
	       outImageIdx[3] = rSampleMap[nodeIt].find_first();
	       outImagePtr->SetPixel(outImageIdx, 1);
	  }
     } // nodeIt

     WriterType4DS::Pointer writer = WriterType4DS::New();
     writer->SetInput( outImagePtr );
     std::string outGrpFile = par.rsampledir + "/" + "grp.nii.gz";
     writer->SetFileName(outGrpFile);

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
     std::cout << "SaveRunnignSample(): File  " << outGrpFile << " saved.\n";
     return 0;
}

int PrintPar(unsigned prlevel, ParStruct & par)
{
     if (prlevel >= 1) {
	  // print VMF's mu.
	  printf("%-11s", "mu");
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       printf("%10s%02d", "mu_", clsIdx+1);
	  }
	  printf("\n");
	  
	  for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	       printf("%s%-9d", "sub", subIdx+1);
	       for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {

		    printf("[%3.2f %3.2f] ", par.vmm.sub[subIdx].comp[clsIdx].mu[0], par.vmm.sub[subIdx].comp[clsIdx].mu[1]);
	       }
	       printf("\n");
	  } // subIdx
	  printf("\n");

	  // print VMF's kappa.
	  printf("%-11s", "kappa");
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       printf("%10s%02d", "kappa_", clsIdx+1);
	  }
	  printf("\n");

	  
	  for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	       printf("%s%-8d", "sub", subIdx+1);
	       for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {

		    printf("%12.4f", par.vmm.sub[subIdx].comp[clsIdx].kappa);
	       }
	       printf("\n");
	  } // subIdx
	  printf("\n");

	  // print the number of points in each sub's each clusters.
	  printf("%-11s", "numPts");
	  for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
	       printf("%10s%02d", "numPts_", clsIdx+1);
	  }
	  printf("\n");

	  for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	       printf("%s%-8d", "sub", subIdx+1);
	       for (unsigned clsIdx = 0; clsIdx < par.numClusters; clsIdx ++) {
		    printf("%12ld", par.vmm.sub[subIdx].comp[clsIdx].numPts);
	       }
	       printf("\n");
	  } // subIdx
	  printf("\n");
	  
     }

     if (prlevel >=2) {
	  printf("numClusters: %d,   numSamples: %d,   numSubs: %d,   burnin: %d,   \ninitTemp: %4.3f,   finalTemp: %4.3f, temperature: %4.3f, tsLength: %d\n", par.numClusters, par.numSamples, par.numSubs, par.burnin, par.initTemp, par.finalTemp, par.temperature, par.tsLength);
	  printf("alpha: %4.3f,   beta: %4.3f,   gamma:%4.3f,   alpha0: %4.3f,   sigma2: %4.3f,   nthreads: %d,   sweepPerThread: %d\n", par.alpha, par.beta, par.gamma, par.alpha0, par.sigma2, par.nthreads, par.sweepPerThread);
     }
     if (prlevel >=4 ) {
	  printf("subject names:\n");
	  for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	       printf("%s%02d: %s\n", "sub", subIdx+1, par.vmm.sub[subIdx].name.c_str());
	  }	  
     }


     return 0;
}


int save3dcharInc(ImageType3DChar::Pointer ip, std::string fname)


{
     // use this filter to Add 1 to all labels.
     AddConstantToImageFilterType::Pointer myfilter = AddConstantToImageFilterType::New();
     myfilter->SetInput(ip);
     myfilter->SetConstant(1);
     myfilter->Update();
     
     // use this filter to cast char to short int for saving file. Original data
     // should not be changed.
     typedef itk::CastImageFilter< ImageType3DChar, ImageType3DShort > CastFilterType;
     CastFilterType::Pointer castFilter = CastFilterType::New();
     castFilter->SetInput(myfilter->GetOutput());

     WriterType3DShort::Pointer writer = WriterType3DShort::New();
     writer->SetInput(castFilter->GetOutput() );
     writer->SetFileName(fname);

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

     std::cout << "save3dcharInc(): File " << fname << " saved.\n";

     return 0;
}


int printVnlVector(vnl_vector<float> vec, unsigned numElements)
{
     unsigned short eleIdx = 0;
     unsigned vecSize = vec.size();
     unsigned numPrint = vnl_math_min(vecSize, numElements);
     printf("[ ");
     for (eleIdx = 0; eleIdx < numPrint; eleIdx ++) {
	  // cout << vec[eleIdx] << " ";
	  printf("%.4f ", vec[eleIdx]);
     }
     printf("]\n");
     return (0);
}

int SaveRunningSamples2(lemon::SmartGraph & theGraph, 
			lemon::SmartGraph::NodeMap<SuperCoordType> & coordMap,
			lemon::SmartGraph::NodeMap< boost::dynamic_bitset<> > & rSampleMap,
			ParStruct & par,
			unsigned scanIdx)

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

     std::string thisSubFilename;
     std::ostringstream str_beta, str_scan;
     str_beta << par.beta;
     str_scan << std::setw(5) << std::setfill('0') << scanIdx;
     
     for (unsigned subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  thisSubFilename = par.rsampledir + "/" +  par.vmm.sub[subIdx].name + "_cls_" +  std::to_string(static_cast<long long>(par.numClusters)) +"_beta_" + str_beta.str() +"_scan_" + str_scan.str() + ".nii.gz";
	  outImagePtr->FillBuffer(0);
	  for (lemon::SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=lemon::INVALID; ++ nodeIt) {
	       // if this is subject node.
	       if (coordMap[nodeIt].subid == subIdx) {
		    outImageIdx[0] = coordMap[nodeIt].idx[0];
		    outImageIdx[1] = coordMap[nodeIt].idx[1];
		    outImageIdx[2] = coordMap[nodeIt].idx[2];
		    outImageIdx[3] = rSampleMap[nodeIt].find_first();
		    outImagePtr->SetPixel(outImageIdx, 1);
	       }
	  } // nodeIt

	  WriterType4DS::Pointer writer = WriterType4DS::New();
	  writer->SetInput( outImagePtr );
	  writer->SetFileName(thisSubFilename);
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
	  std::cout << "SaveRunnignSample2(): File  " << thisSubFilename << " saved.\n";
     } // subIdx

     // Now it's for group 
     outImagePtr->FillBuffer(0);
     for (lemon::SmartGraph::NodeIt nodeIt(theGraph); nodeIt !=lemon::INVALID; ++ nodeIt) {
	  // if this is group node.
	  if (coordMap[nodeIt].subid == par.numSubs) {
	       outImageIdx[0] = coordMap[nodeIt].idx[0];
	       outImageIdx[1] = coordMap[nodeIt].idx[1];
	       outImageIdx[2] = coordMap[nodeIt].idx[2];
	       outImageIdx[3] = rSampleMap[nodeIt].find_first();
	       outImagePtr->SetPixel(outImageIdx, 1);
	  }
     } // nodeIt

     WriterType4DS::Pointer writer = WriterType4DS::New();
     writer->SetInput( outImagePtr );
     std::string outGrpFile = par.rsampledir + "/" + "grp" "_cls_" +  std::to_string(static_cast<long long>(par.numClusters))+ "_beta_" + str_beta.str()  +"_scan_" + str_scan.str() + ".nii.gz";
     writer->SetFileName(outGrpFile);

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
     std::cout << "SaveRunnignSample2(): File  " << outGrpFile << " saved.\n";
     return 0;
}
