/*=========================================================================
 * test to do supervoxel for multimodality volumes
 *=========================================================================*/

// add common libraries for super pixel
#include <commonSupervoxel.h>

// add library for super pixel
#include <string>
#include <MSLIC.h>

// for using boost
//namespace po = boost::program_options;

int main(int argc, char * argv[])
{

  if( argc < 8 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputT1Volume inputT2Volume inputGREVolume inputFLAIRVolume inputSWIVolume outputImageFile SuperVoxelCubeSize CompactParameter " << std::endl;
    return EXIT_FAILURE;
    }

  // define input & output pixel type
  typedef  unsigned short  InputPixelType;
  typedef  unsigned short  OutputPixelType;

  // define input & output image type
  typedef itk::Image< InputPixelType,  3 >   InputImageType;
  typedef itk::Image< OutputPixelType, 3 >   OutputImageType;

  // define reader and writer type
  typedef itk::ImageFileReader< InputImageType >  ReaderType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  // read the first volume
  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();
  ReaderType::Pointer reader3 = ReaderType::New();
  ReaderType::Pointer reader4 = ReaderType::New();
  ReaderType::Pointer reader5 = ReaderType::New();
  
  // First the image type should be declared
  typedef itk::Image< unsigned short, 3 > ImageType;

  // Then the image object can be created
  ImageType::Pointer image = ImageType::New();
  ImageType::Pointer imageInput2 = ImageType::New();
  ImageType::Pointer imageInput3 = ImageType::New();
  ImageType::Pointer imageInput4 = ImageType::New();
  ImageType::Pointer imageInput5 = ImageType::New();
  ImageType::Pointer OutputImage = ImageType::New();

  // asign image
  reader->SetFileName( argv[1] );
  reader->Update();
  image = reader->GetOutput();
  reader2->SetFileName( argv[2] );
  reader2->Update();
  imageInput2 = reader2->GetOutput();
  reader3->SetFileName( argv[3] );
  reader3->Update();
  imageInput3 = reader3->GetOutput();
  reader4->SetFileName( argv[4] );
  reader4->Update();
  imageInput4 = reader4->GetOutput();
  reader5->SetFileName( argv[5] );
  reader5->Update();
  imageInput5 = reader5->GetOutput();

  // Read information
  const ImageType::SizeType sizeOfImage = image->GetLargestPossibleRegion().GetSize();

  int width = sizeOfImage[0];
  int height = sizeOfImage[1];
  int slice = sizeOfImage[2];

  // print the image size
  std::cout << "Volume Size:";
  std::cout << "width: " <<  sizeOfImage[0] << ", height: " << sizeOfImage[1] << ", slice: " << sizeOfImage[2] << std::endl;

  int SliceSize = width*height;
  int** pbuffT1;
  pbuffT1 = new int*[slice];
  int** pbuffT2;
  pbuffT2 = new int*[slice];
  int** pbuffGRE;
  pbuffGRE = new int*[slice];
  int** pbuffFLAIR;
  pbuffFLAIR = new int*[slice];
  int** pbuffSWI;
  pbuffSWI = new int*[slice];
  
  int i,j,k,s;
  ImageType::IndexType pixelIndex;
  ImageType::PixelType pixelValue;

  std::cout<<"Input origin is " << image->GetOrigin()<<std::endl;

  // ---------------------->>>>
  // TO-DO: load all five modalities
  // ---------------------->>>>
  for(k=0; k<slice; k++)
  {
       //std::cout << "current k: " << k << std::endl;

       pbuffT1[k] = new int[SliceSize];
       pbuffT2[k] = new int[SliceSize];
       pbuffGRE[k] = new int[SliceSize];
       pbuffFLAIR[k] = new int[SliceSize];
       pbuffSWI[k] = new int[SliceSize];
       for(j=0; j< height; j++)
       {
	    for(i=0; i< width; i++)
	    {
		 pixelValue = 0;
		 pixelIndex[0] = i;
		 pixelIndex[1] = j;
		 pixelIndex[2] = k;
		 int CurIndex = j*width + i;
		 // for T1
		 pixelValue = image->GetPixel( pixelIndex );
		 pbuffT1[k][CurIndex] = pixelValue;
		 // for T2
		 pixelValue = imageInput2->GetPixel( pixelIndex );
		 pbuffT2[k][CurIndex] = pixelValue;
		 // for GRE
		 pixelValue = imageInput3->GetPixel( pixelIndex );
		 pbuffGRE[k][CurIndex] = pixelValue;
		 // for FLAIR
		 pixelValue = imageInput4->GetPixel( pixelIndex );
		 pbuffFLAIR[k][CurIndex] = pixelValue;
		 // for SWI
		 pixelValue = imageInput5->GetPixel( pixelIndex );
		 pbuffSWI[k][CurIndex] = pixelValue;

		 if(k==176)
		 {
		      std::cout << "Index I: " << i << " Index J: " << j << std::endl;
		      std::cout << "current index: " << CurIndex << std::endl;
		 }
	    }
       }
  }
  
  
  // Call the Multi-channel supervoxel function
  
  //----------------------------------
  // Initialize parameters
  //----------------------------------
  //int SizeSuperVoxelCube = 500;//Desired number of superpixels.
  int SizeSuperVoxelCube = atoi(argv[7]);
  //double VarCompact = 35;//Compactness factor. use a value ranging from 10 to 40 depending on your needs. Default is 10
  double VarCompact = atof(argv[8]);

  std::cout << "Input Number of supervoxels: " << SizeSuperVoxelCube << std::endl;
  std::cout << "Input Compactness parameter: " << VarCompact << std::endl;

  int** klabels = new int*[slice];
  for( int d = 0; d < slice; d++ )
  {
       klabels[d] = new int[SliceSize];
  }

  int numlabels(0);
  
  //----------------------------------
  // Perform MSLIC on the image buffer
  //----------------------------------

  // ---------------------->>>>
  // TO-DO: pass all five modalities to super voxel volume
  // ---------------------->>>>
  std::cout << "Do multimodality super voxel! " << std::endl;
  DoSupervoxelSegmentation(pbuffT1, pbuffT2, pbuffGRE, pbuffFLAIR, pbuffSWI, width, height, slice, klabels, numlabels, SizeSuperVoxelCube, VarCompact);

  // set the value of output image
  OutputImageType::RegionType outRegion = image->GetLargestPossibleRegion();
  OutputImage->SetRegions(outRegion);
  OutputImage->Allocate();
  OutputImage->FillBuffer( 0 );

  OutputImage->SetOrigin(image->GetOrigin() );
  OutputImage->SetSpacing(image->GetSpacing() );
  OutputImage->SetDirection(image->GetDirection() );

  int contNumVolume = 0;
  ImageType::IndexType NewPixelIndex;
  ImageType::PixelType NewPixelValue;

  for(k=0; k<slice; k++)
  {
       for(j=0; j< height; j++)
       {
	    for(i=0; i< width; i++)
	    {
		 

		 int CurIndex = j*width + i;
		 NewPixelValue = klabels[k][CurIndex];

		 NewPixelIndex[0] = i;
		 NewPixelIndex[1] = j;
		 NewPixelIndex[2] = k;
		 OutputImage->SetPixel( NewPixelIndex, NewPixelValue );
	    }
       }
  }

  // specify the name of the output 
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[6] );
  
  writer->SetInput( OutputImage );

  // pass the arg file name to writer
  
  writer->Update();

  std::cout << "SuperVoxel is done!! " << std::endl;

  //if(pbuffT1) delete [] pbuffT1;
  //if(klabels) delete [] klabels;

  return EXIT_SUCCESS;

}

