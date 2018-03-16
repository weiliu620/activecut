#include "commonalt.h"

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

     printf("Image saved.\n");
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
}


int printcls(CompType* cls, int nLabels)
{
     for (int k = 0; k < nLabels; k++) {
	  printf("cluster[%i]: mum of points: %ld, mu = %f, sigma = %f.\n", k, cls[k].numPoints, cls[k].mu, cls[k].sigma); 
     }
     
     return 0;
}

