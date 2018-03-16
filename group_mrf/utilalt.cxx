#include "commonalt.h"
using std::cout;
using std::endl;

int saveimage2d(ImageType2D::Pointer ip, std::string fname)

{
     WriterType::Pointer writer = WriterType::New();
     ImageType2D::Pointer imP = ImageType2D::New();
     
     writer->SetInput(ip);
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
     return 0;
}


int saveimage3d(ImageType3D::Pointer ip, std::string fname)

{
     WriterType3D::Pointer writer = WriterType3D::New();
     ImageType3D::Pointer imP = ImageType3D::New();
     
     writer->SetInput(ip);
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

     std::cout << "saveimage3d(): File " << fname << " saved.\n";
     return 0;
}

int saveimage4dchar(ImageType4DChar::Pointer ip, std::string fname)

{
     WriterType4DChar::Pointer writer = WriterType4DChar::New();
     ImageType4DChar::Pointer imP = ImageType4DChar::New();
     
     writer->SetInput(ip);
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

     std::cout << "saveimage4dchar(): File " << fname << " saved.\n";
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


int saveExtractSlice(ImageType3D::Pointer inPtr,
		     unsigned sliceIdx, std::string fname)
{
     WriterType::Pointer writer = WriterType::New();
     writer->SetFileName(fname);
     typedef itk::ExtractImageFilter< ImageType3D, ImageType2D > FilterType;
     FilterType::Pointer filter = FilterType::New();
     ImageType3D::RegionType inputRegion =
           inPtr->GetLargestPossibleRegion();
     
     ImageType3D::SizeType size = inputRegion.GetSize();
     size[2] = 0;

     ImageType3D::IndexType start = inputRegion.GetIndex();
     start[2] = sliceIdx;

     ImageType3D::RegionType desiredRegion;
     desiredRegion.SetSize(  size  );
     desiredRegion.SetIndex( start );
     filter->SetExtractionRegion( desiredRegion );
     filter->SetInput( inPtr);
     writer->SetInput( filter->GetOutput() );

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
     std::cout << "saveExtractSlice(): File " << fname << " saved.\n";

     return (0);
}


int saveExtractVolume(ImageType4DChar::Pointer inPtr,
		     unsigned volumeIdx, std::string fname)
{
     WriterType3DChar::Pointer writer = WriterType3DChar::New();
     writer->SetFileName(fname);
     typedef itk::ExtractImageFilter< ImageType4DChar, ImageType3DChar > FilterType;
     FilterType::Pointer filter = FilterType::New();
     ImageType4DChar::RegionType inputRegion =
           inPtr->GetLargestPossibleRegion();
     
     ImageType4DChar::SizeType size = inputRegion.GetSize();
     size[3] = 0;

     ImageType4DChar::IndexType start = inputRegion.GetIndex();
     start[3] = volumeIdx;

     ImageType4DChar::RegionType desiredRegion;
     desiredRegion.SetSize(  size  );
     desiredRegion.SetIndex( start );
     filter->SetExtractionRegion( desiredRegion );
     filter->SetInput( inPtr);
     writer->SetInput( filter->GetOutput() );

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
     std::cout << "saveExtractVolume(): File " << fname << " saved.\n";

     return (0);
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

int printVnlMatrix(vnl_matrix<float> mat, unsigned numElements)
{

     unsigned maxRow = vnl_math_min(numElements, mat.rows());
     unsigned maxCol = vnl_math_min(numElements, mat.columns());

     unsigned short rowIdx = 0, colIdx = 0;

     for (rowIdx = 0; rowIdx < maxRow; rowIdx ++) {
	  for (colIdx = 0; colIdx < maxCol; colIdx ++) {
	       cout << mat(rowIdx, colIdx) << " ";
	  }
	  cout << endl;
     }
     cout << endl;
     return (0);
}

int printVnlMatrixSci(vnl_matrix<float> mat, unsigned numElements)
{

     unsigned maxRow = vnl_math_min(numElements, mat.rows());
     unsigned maxCol = vnl_math_min(numElements, mat.columns());

     unsigned short rowIdx = 0, colIdx = 0;

     for (rowIdx = 0; rowIdx < maxRow; rowIdx ++) {
	  for (colIdx = 0; colIdx < maxCol; colIdx ++) {
	       printf("%E  ", mat(rowIdx, colIdx));
	  }
	  cout << endl;
     }
     cout << endl;
     return (0);
}


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


int printvmm(std::vector <VMM> & vmm)
{
     unsigned numSubs = vmm.size();
     unsigned numClusters = vmm[0].comp.size();

     unsigned int clsIdx = 0, subIdx = 0;

     for (subIdx = 0; subIdx < numSubs; subIdx ++) {
	  printf("subject %i: \n", subIdx);
	  for (clsIdx = 0; clsIdx < numClusters; clsIdx++) {
	       printf("       cluster[%i]: # of points: %ld, meanNorm = %f, kappa = %4.2f, label = %i. ", 
		      clsIdx, 
		      vmm[subIdx].comp[clsIdx].numPoints, 
		      vmm[subIdx].comp[clsIdx].meanNorm,
		      vmm[subIdx].comp[clsIdx].kappa,
		      vmm[subIdx].comp[clsIdx].label); 
	       printf("mu = ");
	       printVnlVector(vmm[subIdx].comp[clsIdx].mu, 3);
	  } // clsIdx
     } // subIdx.
     return 0;
}

int SaveSample(ImageType3DChar::Pointer samplePtr, 		 
	       std::string filename)
{
     AddConstantToImageFilterType::Pointer myfilter = AddConstantToImageFilterType::New();
     myfilter->SetInput(samplePtr);
     myfilter->SetConstant(1);
     myfilter->Update();
     
     WriterType3DChar::Pointer writer = WriterType3DChar::New();
     writer->SetInput(myfilter->GetOutput() );
     writer->SetFileName(filename);
     
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

     std::cout << "SaveSample(): File " << filename << " saved.\n";
     return 0;
}

int SaveSamples(std::vector< std::vector<ImageType3DChar::Pointer> > &  sampleVec,
		ImageType3DChar::Pointer maskPtr,
		std::string basename)
{

     ImageType3DChar::SizeType sampleSize = 
	  sampleVec[0][0]->GetLargestPossibleRegion().GetSize();
     unsigned numSubs = sampleVec.size();
     unsigned numSamples = sampleVec[0].size();

     // mask.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() ) ;
     ImageType3DChar::IndexType maskIdx;     

     // output.
     ImageType4DFloat::Pointer outPtr = ImageType4DFloat::New();
     ImageType4DFloat::IndexType outIdx;
     outIdx.Fill(0);
     ImageType4DFloat::SizeType outSize;

     outSize[0] = sampleSize[0];
     outSize[1] = sampleSize[1];
     outSize[2] = sampleSize[2];
     outSize[3]= numSamples;
     ImageType4DFloat::RegionType outRegion;
     outRegion.SetSize(outSize);
     outRegion.SetIndex(outIdx);
     outPtr->SetRegions(outRegion);
     outPtr->Allocate();
     outPtr->FillBuffer(0);

     ImageType4DFloat::PointType outOrigin;
     ImageType4DFloat::SpacingType outSpacing;
     ImageType4DFloat::DirectionType outDirection;

     for (unsigned subIdx = 0; subIdx < numSubs; subIdx ++) {
	  std::string subFileNum = boost::lexical_cast<std::string> (subIdx+1);
	  if (subIdx < 9) {
	       subFileNum.insert(0,"0");
	  }

	  std::string thisSubFilename(basename);
	  // add sub number.
	  thisSubFilename.append(subFileNum);
	  thisSubFilename.append(".nii.gz");

	  outPtr->FillBuffer(0);
	  for (unsigned mcIdx = 0; mcIdx < numSamples; mcIdx ++) {

	       ConstIteratorType3DChar sampleIt(
		    sampleVec[subIdx][mcIdx],
		    sampleVec[subIdx][mcIdx]->GetLargestPossibleRegion() ) ;

	       outIdx[3] = mcIdx;
	       
	       for (maskIt.GoToBegin(), sampleIt.GoToBegin(); !sampleIt.IsAtEnd(); ++ sampleIt, ++ maskIt) {
		    if (maskIt.Get() > 0) {
			 maskIdx = maskIt.GetIndex();
			 outIdx[0] = maskIdx[0];
			 outIdx[1] = maskIdx[1];
			 outIdx[2] = maskIdx[2];
			 
			 outPtr->SetPixel(outIdx, sampleIt.Get() + 1 );
		    } // masIt > 0
	       } // maskIt
	  } // mcIdx.

	  // Write to file
	  WriterType4DFloat::Pointer writer = WriterType4DFloat::New();
	  
	  writer->SetInput(outPtr);
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
     } // subIdx.
     std::cout << "SaveSamples(): File with basename " << basename << " saved.\n";

     return 0;
}


int SaveGrpSamples(std::vector<ImageType3DChar::Pointer>  &  grpVec,
		ImageType3DChar::Pointer maskPtr,
		std::string outFilename)
{

     ImageType3DChar::SizeType sampleSize = 
	  grpVec[0]->GetLargestPossibleRegion().GetSize();
     unsigned numSamples = grpVec.size();

     // mask.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() ) ;
     ImageType3DChar::IndexType maskIdx;     

     // output.
     ImageType4DFloat::Pointer outPtr = ImageType4DFloat::New();
     ImageType4DFloat::IndexType outIdx;
     outIdx.Fill(0);
     ImageType4DFloat::SizeType outSize;

     outSize[0] = sampleSize[0];
     outSize[1] = sampleSize[1];
     outSize[2] = sampleSize[2];
     outSize[3]= numSamples;
     ImageType4DFloat::RegionType outRegion;
     outRegion.SetSize(outSize);
     outRegion.SetIndex(outIdx);
     outPtr->SetRegions(outRegion);
     outPtr->Allocate();
     outPtr->FillBuffer(0);

     for (unsigned mcIdx = 0; mcIdx < numSamples; mcIdx ++) {

	  ConstIteratorType3DChar srcIt(
	       grpVec[mcIdx], grpVec[mcIdx]->GetLargestPossibleRegion() ) ;

	  outIdx[3] = mcIdx;
	       
	  for (maskIt.GoToBegin(), srcIt.GoToBegin(); !maskIt.IsAtEnd(); ++ srcIt, ++ maskIt) {
	       if (maskIt.Get() > 0) {
		    maskIdx = maskIt.GetIndex();
		    outIdx[0] = maskIdx[0];
		    outIdx[1] = maskIdx[1];
		    outIdx[2] = maskIdx[2];
			 
		    outPtr->SetPixel(outIdx, srcIt.Get() + 1 );
	       } // masIt > 0
	  } // maskIt
     } // mcIdx.

     // Write to file
     WriterType4DFloat::Pointer writer = WriterType4DFloat::New();
	  
     writer->SetInput(outPtr);
     writer->SetFileName(outFilename);
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

     std::cout << "SaveGrpSamples(): File " << outFilename << " saved.\n";

     return 0;
}

