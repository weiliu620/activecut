/*=========================================================================
 * main function to do supervoxel for multimodality volumes
 *=========================================================================*/

// add common libraries for super pixel
#include <common.h>

// add library for super pixel
#include <MSLIC.h>

namespace po = boost::program_options;


int main(int argc, char * argv[])
{

     /*
     if( argc < 8 )
     {
	  std::cerr << "Usage: " << argv[0];
	  std::cerr << " inputT1Volume inputT2Volume inputGREVolume inputFLAIRVolume inputSWIVolume outputImageFile SuperVoxelCubeSize CompactParameter " << std::endl;
	  return EXIT_FAILURE;
     }
     */

     std::string input_image_file, output_image_file;

     int SuperVoxelCubeSize = 0;
     int CompactParameter = 0;

     // program options.
     po::options_description svdesc("Options (to use supervoxel) can only used at commandline");
     svdesc.add_options()
	  
	  ("help,h", "Supervoxel function for multi-modal 3D volume data")

	  ("data,d", po::value<std::string>(&input_image_file),
	   "Input all-channel image file. A 4D gipl or nii or nii.gz file.")
	  
	  ("output,o", po::value<std::string>(&output_image_file),
	   "Output label image file. A 3D nii file.")
	  
	  ("supervoxelsize,s", po::value<int>(&SuperVoxelCubeSize)->default_value(500),
	   "Volume size of the super voxel cube, e.g., 3*3*3 = 27.")

	  ("compact,c", po::value<int>(&CompactParameter)->default_value(35),
	   "Compactness parameter.")
	  ;

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, svdesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: Supervoxel [options]\n";
	       std::cout << svdesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }

     std::cout << "came to here!!" << "\n";

     // input image
     typedef itk::ImageFileReader< ImageType4DF >  ReaderType;
     ReaderType::Pointer reader = ReaderType::New();
     reader->SetFileName( input_image_file );
     reader->Update();
     
     ImageType4DF::Pointer inputimage = ImageType4DF::New();
     inputimage = reader->GetOutput();

     const ImageType4DF::SizeType sizeOfImage = inputimage->GetLargestPossibleRegion().GetSize();

     int width = sizeOfImage[0];
     int height = sizeOfImage[1];
     int slice = sizeOfImage[2];
     int numchannel = sizeOfImage[3];
     int SliceSize = width*height;
     
     std::cout << "width: " << width << " ; height: " << height << "\n";
     std::cout << "slice: " << slice << " ; num of modalities: " << numchannel << "\n";
     std::cout << "Supervoxel cube size: " << SuperVoxelCubeSize << "\n";
     std::cout << "Compactness: " << CompactParameter << "\n";
     
     // create image
     int** pbuffAll[5];
     
     
     // extract 3d vol from 4d.
     ImageType4DF::IndexType vol_idx;
     vol_idx.Fill(0);
     ImageType4DF::SizeType vol_size = sizeOfImage;
     vol_size[3] = 0;

     int i,j,k;
     ImageType3DF::IndexType pixelIndex;
     ImageType3DF::PixelType pixelValue;

     ImageType3DF::Pointer FirstModalityVolume;

     // for each channel
     typedef itk::ExtractImageFilter<ImageType4DF, ImageType3DF> ExtFilter;
     for (unsigned p = 0; p < sizeOfImage[3]; p ++) {
	  vol_idx.Fill(0);
	  vol_idx[3] = p;
	  ImageType4DF::RegionType vol_region(vol_idx, vol_size);
	  ExtFilter::Pointer ext_filter = ExtFilter::New();
	  ext_filter->SetExtractionRegion(vol_region);
	  //ext_filter->SetInput(dataPtr);
	  ext_filter->SetInput(inputimage);
	  ext_filter->Update();

	  ext_filter->GetOutput();

	  ImageType3DF::Pointer outputVolume = ext_filter->GetOutput();
	  if(p == 1)
	  {
	       FirstModalityVolume = outputVolume;
	  }

	  pbuffAll[p] = new int*[slice];

	  for(k=0; k<slice; k++)
	  {
	       pbuffAll[p][k] = new int[SliceSize];
	       
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
			 pixelValue = outputVolume->GetPixel( pixelIndex );
			 pbuffAll[p][k][CurIndex] = int(pixelValue);
		    }
	       }
	  }

     }

     // Call the Multi-channel supervoxel function
  
     //----------------------------------
     // Initialize parameters
     //----------------------------------

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
     DoSupervoxelSegmentation(pbuffAll[0], pbuffAll[1], pbuffAll[2], pbuffAll[3], pbuffAll[4], width, height, slice, klabels, numlabels, SuperVoxelCubeSize, CompactParameter);
     

     // output image
     typedef  unsigned short  OutputPixelType;
     typedef itk::Image< OutputPixelType, 3 >   OutputImageType;
     typedef itk::ImageFileWriter< OutputImageType >  WriterType;
     OutputImageType::Pointer OutputImage = OutputImageType::New();
     
     // set the value of output image
     OutputImageType::RegionType outRegion = FirstModalityVolume->GetLargestPossibleRegion();
     OutputImage->SetRegions(outRegion);
     OutputImage->Allocate();
     OutputImage->FillBuffer( 0 );

     OutputImage->SetOrigin(FirstModalityVolume->GetOrigin() );
     OutputImage->SetSpacing(FirstModalityVolume->GetSpacing() );
     OutputImage->SetDirection(FirstModalityVolume->GetDirection() );

     int contNumVolume = 0;
     OutputImageType::IndexType NewPixelIndex;
     OutputImageType::PixelType NewPixelValue;

     for(k=0; k<slice; k++)
     {
	  for(j=0; j< height; j++)
	  {
	       for(i=0; i< width; i++)
	       {
		 
		    int CurIndex = j*width + i;
		    NewPixelValue = klabels[k][CurIndex];
		    //NewPixelValue = pbuffAll[4][k][CurIndex];

		    NewPixelIndex[0] = i;
		    NewPixelIndex[1] = j;
		    NewPixelIndex[2] = k;
		    OutputImage->SetPixel( NewPixelIndex, NewPixelValue );
	       }
	  }
     }

     // specify the name of the output 
     WriterType::Pointer writer = WriterType::New();
     writer->SetFileName( output_image_file );
  
     writer->SetInput( OutputImage );

     // pass the arg file name to writer
  
     writer->Update();

     std::cout << "supervoxelseg(): save output to: " <<  output_image_file << std::endl;
     
     
     //if(pbuffT1) delete [] pbuffT1;
     //if(klabels) delete [] klabels;

     return EXIT_SUCCESS;

}
