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

int save3dchar(ImageType3DChar::Pointer ip, std::string fname)


{

     WriterType3DChar::Pointer writer = WriterType3DChar::New();
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

     std::cout << "save3dchar(): File " << fname << " saved.\n";

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

