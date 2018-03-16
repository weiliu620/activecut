#include "common_mcem.h"
using std::cout;
using std::endl;

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

int save3dchar(ImageType3DChar::Pointer ip, std::string fname)


{

     AddConstantToImageFilterType::Pointer myfilter = AddConstantToImageFilterType::New();
     myfilter->SetInput(ip);
     myfilter->SetConstant(1);
     myfilter->Update();
     
     WriterType3DChar::Pointer writer = WriterType3DChar::New();
     writer->SetInput(myfilter->GetOutput() );
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

     std::cout << "save3dchar(): File " << fname << " saved.\n";

     return 0;
}

int printVnlVector(vnl_vector<float> vec, unsigned numElements)
{
     unsigned short eleIdx = 0;
     unsigned vecSize = vec.size();
     unsigned numPrint = vnl_math_min(vecSize, numElements);
     for (eleIdx = 0; eleIdx < numPrint; eleIdx ++) {
	  cout << vec[eleIdx] << " ";
     }
     cout << endl;
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

     std::cout << "SaveSamples(): File " << basename << " saved.\n";
     return 0;
}


int ReadInitGrpLabel(ImageType3DVecF::Pointer outPtr, 
		     ImageType3DChar::Pointer initGrpPtr)
{
     ImageType3DChar::RegionType grpRegion = initGrpPtr->GetLargestPossibleRegion();
     ConstIteratorType3DChar grpIt(initGrpPtr, initGrpPtr->GetLargestPossibleRegion() );

     IteratorType3DVecF outIt(outPtr, outPtr->GetLargestPossibleRegion() );
     unsigned int numClusters = outPtr->GetNumberOfComponentsPerPixel();

     itk::VariableLengthVector<float> grpVector(numClusters);
     unsigned short label = 0;

     // before read in value. Fill 0 vector to all voxels.
     grpVector.Fill(0);
     outPtr->FillBuffer(grpVector);

     // read data from file.
     for (grpIt.GoToBegin(), outIt.GoToBegin(); !grpIt.IsAtEnd(); 
	  ++grpIt, ++outIt) {
	  if (grpIt.Get() > 0) {
	       label = grpIt.Get();
	       grpVector.Fill(0);
	       grpVector[label-1] = 1;
	       outIt.Set(grpVector); 
	  }
     }
     return 0;
}



unsigned short GetGrpLabel(ImageType3DVecF::Pointer groupPtr,      
			   ImageType3DVecF::IndexType groupIdx)
{
     ImageType3DVecF::PixelType pixelValue = groupPtr->GetPixel(groupIdx);
     unsigned short clsIdx = 0;
     float maxPr = 0;
     for (int k = 0; k < pixelValue.GetSize(); k ++) {
	  if (pixelValue[k] > maxPr) {
	       maxPr = pixelValue[k];
	       clsIdx = k;
	  }
     }
     return clsIdx;
}

int SetGrpLabel(ImageType3DVecF::Pointer groupPtr,      
		ImageType3DVecF::IndexType groupIdx,
		unsigned short label)
{
     ImageType3DVecF::PixelType pixelValue = groupPtr->GetPixel(groupIdx);
     for (int k = 0; k < pixelValue.GetSize(); k ++) {
	  pixelValue[k] = 0;
     }
     pixelValue[label] = 1;
     groupPtr->SetPixel(groupIdx, pixelValue);
     return 0;
}

int SaveGrpLabel(ImageType3DVecF::Pointer groupPtr, 
		 ImageType3DChar::Pointer maskPtr,
		 std::string filename)
{
     ImageType3DChar::Pointer outPtr = ImageType3DChar::New();
     ImageType3DChar::IndexType outIdx;
     outIdx.Fill(0);
     ImageType3DChar::SizeType outSize;
     ImageType3DVecF::SizeType groupSize = groupPtr->GetLargestPossibleRegion().GetSize();
     outSize = groupSize;
     ImageType3DChar::RegionType outRegion;
     outRegion.SetSize(outSize);
     outRegion.SetIndex(outIdx);
     outPtr->SetRegions(outRegion);
     outPtr->Allocate();
     outPtr->FillBuffer(0);

     ImageType3DChar::PointType outOrigin;
     ImageType3DChar::SpacingType outSpacing;
     ImageType3DChar::DirectionType outDirection;
     ImageType3DVecF::IndexType groupIdx;

     ImageType3DVecF::PointType groupOrigin = groupPtr->GetOrigin();
     ImageType3DVecF::SpacingType groupSpacing = groupPtr->GetSpacing();
     ImageType3DVecF::DirectionType groupDirection = groupPtr->GetDirection();


     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     ImageType3DChar::RegionType maskRegion = maskPtr->GetLargestPossibleRegion();
     IteratorType3DChar maskIt(maskPtr, maskRegion);
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();

     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       maskIdx = maskIt.GetIndex();
	       groupIdx[0] = maskIdx[0];
	       groupIdx[1] = maskIdx[1];
	       groupIdx[2] = maskIdx[2];

	       outIdx[0] = maskIdx[0];
	       outIdx[1] = maskIdx[1];
	       outIdx[2] = maskIdx[2];

	       outPtr->SetPixel(outIdx, GetGrpLabel(groupPtr, groupIdx) + 1 );
	  }
     }

     // Write to file
     WriterType3DChar::Pointer writer = WriterType3DChar::New();
     
     writer->SetInput(outPtr);
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

     std::cout << "SaveGrpLabel(): File " << filename << " saved.\n";
     return 0;
     
}


int SaveGrpProb(ImageType3DVecF::Pointer groupPtr, 
		ImageType3DChar::Pointer maskPtr,
		std::string filename)
{
     ImageType3DVecF::IndexType groupIdx;
     itk::VariableLengthVector<float> groupValue;
     ImageType3DVecF::SizeType groupSize 
	  = groupPtr->GetLargestPossibleRegion().GetSize();
     groupIdx.Fill(0);
     unsigned short numClusters = groupPtr->GetNumberOfComponentsPerPixel();

     ImageType4DFloat::Pointer outPtr = ImageType4DFloat::New();
     ImageType4DFloat::IndexType outIdx;
     outIdx.Fill(0);
     ImageType4DFloat::SizeType outSize;

     outSize[0] = groupSize[0];
     outSize[1] = groupSize[1];
     outSize[2] = groupSize[2];
     outSize[3]= numClusters;
     ImageType4DFloat::RegionType outRegion;
     outRegion.SetSize(outSize);
     outRegion.SetIndex(outIdx);
     outPtr->SetRegions(outRegion);
     outPtr->Allocate();
     outPtr->FillBuffer(0);

     ImageType4DFloat::PointType outOrigin;
     ImageType4DFloat::SpacingType outSpacing;
     ImageType4DFloat::DirectionType outDirection;

     ImageType3DVecF::PointType groupOrigin = groupPtr->GetOrigin();
     ImageType3DVecF::SpacingType groupSpacing = groupPtr->GetSpacing();
     ImageType3DVecF::DirectionType groupDirection = groupPtr->GetDirection();

     // origin, orientation, and spacing.
     outOrigin[0] = groupOrigin[0];
     outOrigin[1] = groupOrigin[1];
     outOrigin[2] = groupOrigin[2];
     outSpacing[0] = groupSpacing[0];
     outSpacing[1] = groupSpacing[1];
     outSpacing[2] = groupSpacing[2];
     outSpacing[3] = 3;

     outDirection(0,0) = groupDirection(0,0);
     outDirection(0,1) = groupDirection(0,1);
     outDirection(0,2) = groupDirection(0,2);
     outDirection(1,0) = groupDirection(1,0);
     outDirection(1,1) = groupDirection(1,1);
     outDirection(1,2) = groupDirection(1,2);
     outDirection(2,0) = groupDirection(2,0);
     outDirection(2,1) = groupDirection(2,1);
     outDirection(2,2) = groupDirection(2,2);

     outDirection(3,3) = 1;
     
     outPtr->SetOrigin(outOrigin);
     outPtr->SetSpacing(outSpacing);
     outPtr->SetDirection(outDirection);

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     ImageType3DChar::RegionType maskRegion = maskPtr->GetLargestPossibleRegion();
     IteratorType3DChar maskIt(maskPtr, maskRegion);
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();

     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       maskIdx = maskIt.GetIndex();

		    if (maskIdx[0] == 7 && maskIdx[1] == 23 && maskIdx[2] == 10) {
			 printf(" ");
		    }


	       groupIdx[0] = maskIdx[0];
	       groupIdx[1] = maskIdx[1];
	       groupIdx[2] = maskIdx[2];

	       outIdx[0] = maskIdx[0];
	       outIdx[1] = maskIdx[1];
	       outIdx[2] = maskIdx[2];

	       groupValue = groupPtr->GetPixel(groupIdx);
	       for (outIdx[3] = 0; outIdx[3] < numClusters; outIdx[3] ++) {
		    outPtr->SetPixel(outIdx, groupValue[outIdx[3]]);
	       }
	  }
     }

     // Write to file
     WriterType4DFloat::Pointer writer = WriterType4DFloat::New();
     
     writer->SetInput(outPtr);
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

     std::cout << "SaveGrpProb(): File " << filename << " saved.\n";
     return 0;
     
}



unsigned ReadObs(std::string obspath,
	    VnlVectorImageType::Pointer & obsPtr,
	    ImageType3DChar::Pointer maskPtr)
{
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     boost::filesystem::path obspathVar(obspath);
     boost::filesystem::directory_iterator obspathEnd;

     ReaderType3DFloat::Pointer obsReader = ReaderType3DFloat::New();

     // Get number of subjects.
     unsigned numSubs = 0;
     for (boost::filesystem::directory_iterator obspathIt(obspathVar); obspathIt != obspathEnd; ++obspathIt) {
	  numSubs ++;
     }

     vnl_vector<float> obsVector(numSubs, 0);
     ImageType3DFloat::Pointer singleObsPtr = ImageType3DFloat::New();

     // allocate memory for obsPtr.
     VnlVectorImageType::IndexType obsIdx;
     obsIdx.Fill ( 0 );
     VnlVectorImageType::SizeType obsSize;

     // get image size.
     boost::filesystem::directory_iterator obspathIt(obspathVar);
     obsReader->SetFileName( (*obspathIt).path().string() );
     obsReader->Update();
     singleObsPtr = obsReader->GetOutput();
     obsSize = singleObsPtr->GetLargestPossibleRegion().GetSize();

     // region
     VnlVectorImageType::RegionType obsRegion;
     obsRegion.SetSize (obsSize);
     obsRegion.SetIndex( obsIdx );
     obsPtr->SetRegions( obsRegion );
     obsPtr->Allocate();
     obsPtr->FillBuffer ( obsVector );

     IteratorTypeVnlVector obsIt(obsPtr, obsPtr->GetLargestPossibleRegion() );
     unsigned subIdx = 0;

     for (boost::filesystem::directory_iterator obspathIt(obspathVar); obspathIt != obspathEnd; ++obspathIt) {
	  obsReader->SetFileName( (*obspathIt).path().string() );
	  obsReader->Update();
	  singleObsPtr = obsReader->GetOutput();
	  IteratorType3DFloat singleObsIt(singleObsPtr, singleObsPtr->GetLargestPossibleRegion() );
	  std::cout <<  "add " << (*obspathIt).path().string() << "\n";

	  // read in data for this subject.

	  for (maskIt.GoToBegin(), obsIt.GoToBegin(); !obsIt.IsAtEnd(); ++ obsIt, ++maskIt){
	       if (maskIt.Get() > 0) {

		    obsIdx = maskIt.GetIndex();
		    obsVector = obsIt.Get();
		    obsVector[subIdx] = singleObsPtr->GetPixel(obsIdx);		    
		    obsIt.Set( obsVector );
	       } // maskIt > 0
	  } // maskIt
	  
	  // done with read data for this sub. 
	  subIdx ++;
     }
     return numSubs;
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

