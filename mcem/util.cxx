#include "common.h"
#include "itkRescaleIntensityImageFilter.h"

int SaveToy(MatrixXf image, string fname)

{
     WriterType::Pointer writer = WriterType::New();
     ImageType2D::Pointer imP = ImageType2D::New();
     ImageType2D::IndexType start;               
     
     start[0] =   0;  // first index on X
     start[1] =   0;  // first index on Y

     ImageType2D::SizeType  size;
     size[0]  = image.rows();  // size along X
     size[1]  = image.cols();  // size along Y
//     printf("in Savetoy, image.rows = %d, image.cols = %d\n", image.rows(), image.cols());
     printf("Image saved.\n");

     ImageType2D::RegionType region;
     region.SetSize( size );
     region.SetIndex( start );

     imP->SetRegions( region );
     imP->Allocate();

     // The image buffer is initialized to a particular value
     imP->FillBuffer(0);

     ImageType2D::IndexType pixelIndex;

     for (pixelIndex[0] = 0; pixelIndex[0] < image.rows(); pixelIndex[0]++){
          for (pixelIndex[1] = 0; pixelIndex[1] < image.cols(); pixelIndex[1]++){
	       imP->SetPixel(pixelIndex, image(pixelIndex[0], pixelIndex[1]));
	  }
     }

     writer->SetInput(imP);
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


     // Write to a png image.
     typedef unsigned char WriterPixelType;
     typedef itk::Image<WriterPixelType, 2> WriteImageType;
     typedef itk::ImageFileWriter<WriteImageType> PNGWriterType;
     typedef itk::RescaleIntensityImageFilter<ImageType2D, WriteImageType> RescaleFilterType;
     RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
     
     rescaler->SetOutputMinimum(0);
     rescaler->SetOutputMaximum(255);
     rescaler->SetInput(imP);
     PNGWriterType::Pointer pngWriter = PNGWriterType::New();
     string pngfilename = fname.append(".png");
     pngWriter->SetFileName(pngfilename);
     pngWriter->SetInput(rescaler->GetOutput() );

     try 
     { 
     	  pngWriter->Update(); 
     } 
     catch( itk::ExceptionObject & err ) 
     { 
     	  std::cerr << "ExceptionObject caught !" << std::endl; 
     	  std::cerr << err << std::endl; 
     	  return EXIT_FAILURE;
     } 



     return 0;
}

int saveimage2d(ImageType2D::Pointer ip, string fname)

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


int saveimage3d(ImageType3D::Pointer ip, string fname)

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


int printcls(CompType* cls, int nLabels)
{
     for (int k = 0; k < nLabels; k++) {
	  printf("cluster[%i]: mum of points: %ld, mu = %f, sigma = %f.\n", k, cls[k].numPoints, cls[k].mu, cls[k].sigma); 
     }
     
     return 0;
}
