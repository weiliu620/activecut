#include <common.h>
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

int saveimage(ImageType3DS::Pointer ip, std::string fname)
{
     WriterType3DS::Pointer writer = WriterType3DS::New();
     ImageType3DS::Pointer imP = ImageType3DS::New();
     
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

int just_test(int a)
{
     return a;
}
