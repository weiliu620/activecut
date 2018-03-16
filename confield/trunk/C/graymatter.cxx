#include "cppheader.h"
#include "common.h"
int graymatter(ConfigFile config)
{
     unsigned short x, y, z;
     double count = 0;

     MaskImType::IndexType maskIdx;

     string datadir;
     config.readInto(datadir, "datadir");
     string grayMaskFile(datadir);
     grayMaskFile.append("/mask.nii");

     // Threshhold for gray matter in mask image.
     double grayThresh = 0;
     config.readInto(grayThresh, "GRAYTHRESH");

     // reader.
     maskReaderType::Pointer maskReader = maskReaderType::New();
     maskReader->SetFileName(grayMaskFile);
     maskReader->Update();

     MaskImType::Pointer maskPtr  =  maskReader->GetOutput();

     MaskImType::SizeType size =
	  maskReader->GetOutput()->GetLargestPossibleRegion().GetSize();

     count = 0;
     for (x = 0; x < size[0]; x++){
	  maskIdx[0] = x;
	  for (y = 0; y < size[1]; y++){
	       maskIdx[1] = y;
	       for (z = 0; z < size[2]; z++){
		    maskIdx[2] = z;
		    if (maskPtr->GetPixel(maskIdx) > grayThresh) {
			 count ++;
		    }
	       }
	  }
     }
     printf(" Gray matter voxel above mask threshhold %f is: %f\n", grayThresh, count);

}
		    
     

