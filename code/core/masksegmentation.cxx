/*=========================================================================
 * load mask to mask supervoxel result
 *=========================================================================*/

// add common libraries for super pixel
#include <common.h>

// add library for super pixel
#include <MSLIC.h>

struct LabelMappingArray{
     int labelValue;
     bool labelExist;
     int newlabelValue;
};


namespace po = boost::program_options;

int main(int argc, char * argv[])
{

     std::string input_image_file, input_mask_file, output_image_file;

     // program options.
     po::options_description svdesc("Options (to use this mask function) can only used at commandline");
     svdesc.add_options()
	  
	  ("help,h", "Mask function for 3D volume data")

	  ("data,d", po::value<std::string>(&input_image_file),
	   "Input 3D volume file. A 3D gipl or nii or nii.gz file.")

	  ("mask,m", po::value<std::string>(&input_mask_file),
	   "Input 3D volume mask file. A 3D gipl or nii or nii.gz file.")
	  
	  ("output,o", po::value<std::string>(&output_image_file),
	   "Output label image file. A 3D nii file.")
	  
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


     ReaderType3DI::Pointer readerVolume = ReaderType3DI::New();
     ReaderType3DI::Pointer readerMask = ReaderType3DI::New();
     

     ImageType3DI::Pointer seglabel = ImageType3DI::New();
     ImageType3DI::Pointer mask = ImageType3DI::New();

     ImageType3DI::Pointer OutputImage = ImageType3DI::New();

     readerVolume->SetFileName( input_image_file );
     readerVolume->Update();
     seglabel = readerVolume->GetOutput();

     readerMask->SetFileName( input_mask_file );
     readerMask->Update();
     mask = readerMask->GetOutput();

     // Read information
     const ImageType3DI::SizeType sizeOfImage = seglabel->GetLargestPossibleRegion().GetSize();
     
     int width = sizeOfImage[0];
     int height = sizeOfImage[1];
     int slice = sizeOfImage[2];

     // print the image size
     std::cout << "Volume Size:";
     std::cout << "width: " <<  width << ", height: " << height << ", slice: " << slice << std::endl;

     int i,j,k;

     ImageType3DI::IndexType pixelIndex;
     ImageType3DI::PixelType pixelValue, maskValue, maxLabelValue;

     ImageType3DI::RegionType outRegion = seglabel->GetLargestPossibleRegion();
     OutputImage->SetRegions(outRegion);
     OutputImage->Allocate();
     OutputImage->FillBuffer( 0 );
     
     // use the mask to mask the supervoxel volume
     // and find the maximum label value
     maxLabelValue = 0;
     for(k=0; k<slice; k++)
     {
	  for(j=0; j< height; j++)
	  {
	       for(i=0; i< width; i++)
	       {
		    pixelValue = 0;
		    maskValue = 0;

		    pixelIndex[0] = i;
		    pixelIndex[1] = j;
		    pixelIndex[2] = k;

		    pixelValue = seglabel->GetPixel( pixelIndex );
		    maskValue = mask->GetPixel( pixelIndex );
		    if(maskValue > 0)
		    {
			 OutputImage->SetPixel( pixelIndex, pixelValue );
			 if(maxLabelValue < pixelValue)
			 {
			      maxLabelValue = pixelValue;
			 }
		    }
	       }
	  }
     }

     // print the maximum label value
     std::cout << "Maximum label value (inside the brain) is: " << maxLabelValue << std::endl;

     // initialize the labe mapping array
     LabelMappingArray alllabels[maxLabelValue];
     int n;
     for(n=0; n<maxLabelValue; n++)
     {
	  alllabels[n].labelValue = n+1;
	  alllabels[n].labelExist = false;
	  alllabels[n].newlabelValue = 0;
     }
     std::cout << "Initialize alllabels is done!" << std::endl;

     // find the mapping of old label value to new ones
     ImageType3DI::PixelType labelValue;
     ImageType3DI::PixelType newlabelValueCumulative = 1;
     for(k=0; k<slice; k++)
     {
	  for(j=0; j< height; j++)
	  {
	       for(i=0; i< width; i++)
	       {
		    labelValue = 0;
		    pixelIndex[0] = i;
		    pixelIndex[1] = j;
		    pixelIndex[2] = k;

		    labelValue = OutputImage->GetPixel( pixelIndex );
		    maskValue = mask->GetPixel( pixelIndex );
		    if(maskValue > 0)
		    {
			 for(n=0; n<maxLabelValue; n++)
			 {
			      if(alllabels[n].labelValue == labelValue && !alllabels[n].labelExist)
			      {
				   alllabels[n].labelExist = true;
				   alllabels[n].newlabelValue = newlabelValueCumulative;
				   newlabelValueCumulative++;
			      }
			 }
		    }
	       }
	  }
     }
     
     // relabel the masked label volume
     ImageType3DI::PixelType newlabelValue;
     for(k=0; k<slice; k++)
     {
	  for(j=0; j< height; j++)
	  {
	       for(i=0; i< width; i++)
	       {
		    labelValue = 0;
		    pixelIndex[0] = i;
		    pixelIndex[1] = j;
		    pixelIndex[2] = k;

		    labelValue = OutputImage->GetPixel( pixelIndex );
		    maskValue = mask->GetPixel( pixelIndex );
		    if(maskValue > 0)
		    {
			 for(n=0; n<maxLabelValue; n++)
			 {
			      if(alllabels[n].labelValue == labelValue && alllabels[n].labelExist)
			      {
				   newlabelValue = alllabels[n].newlabelValue;
				   OutputImage->SetPixel( pixelIndex, newlabelValue );
			      }
			 }
		    }
	       }
	  }
     }

     // test the new maximum label value
     ImageType3DI::PixelType newMaxlabelValue = 0;
     for(k=0; k<slice; k++)
     {
	  for(j=0; j< height; j++)
	  {
	       for(i=0; i< width; i++)
	       {
		    
		    labelValue = 0;
		    pixelIndex[0] = i;
		    pixelIndex[1] = j;
		    pixelIndex[2] = k;

		    labelValue = OutputImage->GetPixel( pixelIndex );
		    if(newMaxlabelValue < labelValue)
		    {
			 newMaxlabelValue = labelValue;
		    }
	       }
	  }
     }

     std::cout << "New maximum label value (inside the brain) is: " << newMaxlabelValue << std::endl;

     WriterType3DI::Pointer writer = WriterType3DI::New();
     writer->SetFileName( output_image_file );
     writer->SetInput( OutputImage );
     writer->Update();

     std::cout << "Mask supervoxel segmentation label is done!! " << std::endl;

     return EXIT_SUCCESS;

}
