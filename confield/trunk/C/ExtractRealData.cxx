#include "cppheader.h"
#include "common.h"

// Extract slices from fmri data file. also return forward map 
// and backward map and region, etc.
void ExtractRealData(float* & im1,
		     float*  & im2,
		     unsigned int*  & fwdMap1,
		     unsigned int*  & fwdMap2,
		     unsigned short*  & bwdMap1,
		     unsigned short*  & bwdMap2,
		     unsigned int& regionSize1,
		     unsigned int& regionSize2,
		     unsigned short* size,
		     ConfigFile config)

{
//     char grayMaskFile[100];
     //   char fmriFile[100];
     string datadir;
     config.readInto(datadir, "datadir");
     string grayMaskFile(datadir);
     string fmriFile(datadir);
     grayMaskFile.append("/mask.nii");
     fmriFile.append("/fmri.nii");

     // reader.
     maskReaderType::Pointer maskReader = maskReaderType::New();
     fmriReaderType::Pointer fmriReader = fmriReaderType::New();
     fmriReader->SetFileName(fmriFile);
     maskReader->SetFileName(grayMaskFile);
     fmriReader->Update();
     maskReader->Update();

     ImageType::Pointer fmriPtr  =  fmriReader->GetOutput();
     MaskImType::Pointer maskPtr  =  maskReader->GetOutput();

     ImageType::SizeType itkSize =
	  fmriReader->GetOutput()->GetLargestPossibleRegion().GetSize();
     size[0] = itkSize[0];
     size[1] = itkSize[1];
     size[2] = itkSize[2];
     size[3] = itkSize[3];

     unsigned short s1 = 0, s2 = 0;
     unsigned short v1_xmin, v1_xmax, v1_ymin, v1_ymax, v1_zmin, v1_zmax;

     unsigned short sliceField, volField;
     config.readInto(sliceField, "sliceField");
     config.readInto(volField, "volField");
     if (sliceField){
	  config.readInto(s1, "s1");
	  config.readInto(s2, "s2");
	  if (_DBG >= 2) {
	       printf("ExtractRealData: s1=%d, s2=%d.\n", s1, s2);
	  }

	  regionSize1 = GetRegionSize(maskReader, GRAYTHRESH, 0, size[0]-1, 0, size[1]-1,
				      s1, s1);
	  regionSize2 = GetRegionSize(maskReader, GRAYTHRESH, 0, size[0]-1, 0, size[1]-1,
				      s2, s2);
     }
     else if (volField){
	  config.readInto(v1_xmin, "v1_xmin");
	  config.readInto(v1_xmax, "v1_xmax");
	  config.readInto(v1_ymin, "v1_ymin");
	  config.readInto(v1_ymax, "v1_ymax");
	  config.readInto(v1_zmin, "v1_zmin");
	  config.readInto(v1_zmax, "v1_zmax");

	  regionSize1 = GetRegionSize(maskReader, GRAYTHRESH, v1_xmin, v1_xmax,
				      v1_ymin, v1_ymax, v1_zmin, v1_zmax);
	  regionSize2 = GetRegionSize(maskReader, GRAYTHRESH, v1_xmin, v1_xmax,
				      v1_ymin, v1_ymax, v1_zmin, v1_zmax);
     }
     else{
	  fprintf(stderr, "either sliceField or volField have to be set. Check parameter file.");
	  exit(1);
     }

     /* Allocate time series in a general matrix. Each row is a time course. NOTE: we assume ROW MAJOR. */
     im1 = (PixelType*) malloc(regionSize1 * size[3] * sizeof(PixelType));
     im2 = (PixelType*) malloc(regionSize2 * size[3] * sizeof(PixelType));

     // Forward mapping matrix. We use calloc to init it to all zero. So in later
     // functoin, only those gray matter are assigned to value. others are still
     // zero. So this forward matrix can be used as a mask for gray matter. 
     fwdMap1 = (unsigned int*) 
	  calloc(size[0]*size[1]*size[2], sizeof(unsigned int));
     fwdMap2 = (unsigned int*) 
	  calloc(size[0]*size[1]*size[2], sizeof(unsigned int));
     bwdMap1 = (unsigned short*) 
	  malloc(regionSize1 * 3 *  sizeof(unsigned short));
     bwdMap2 = (unsigned short*) 
	  malloc(regionSize2 * 3 * sizeof(unsigned short));

     // Extract gray matter voxels into im1 and im2.
     if (sliceField){
	  ExtractGray(fmriReader, maskReader, im1, fwdMap1, bwdMap1, GRAYTHRESH,
		      0, size[0]-1, 0, size[1]-1, s1, s1);
	  ExtractGray(fmriReader, maskReader, im2, fwdMap2, bwdMap2, GRAYTHRESH,
		      0, size[0]-1, 0, size[1]-1, s2, s2);
     }
     else if (volField){
	  printf("ExtractRealData.cxx: mark.\n");
	  printf("v1_xmin=%d, v1_xmax=%d, v1_ymin=%d, v1_ymax=%d, v1_zmin=%d, v1_zmax=%d, v1_xmin=%d, v1_xmax=%d, v1_ymin=%d, v1_ymax=%d, v1_zmin=%d, v1_zmax=%d.\n", v1_xmin, v1_xmax, v1_ymin, v1_ymax, v1_zmin, v1_zmax, v1_xmin, v1_xmax, v1_ymin, v1_ymax, v1_zmin, v1_zmax);
	  ExtractGray(fmriReader, maskReader, im1, fwdMap1, bwdMap1, GRAYTHRESH,
		      v1_xmin, v1_xmax, v1_ymin, v1_ymax, v1_zmin, v1_zmax);
	  ExtractGray(fmriReader, maskReader, im2, fwdMap2, bwdMap2, GRAYTHRESH,
		      v1_xmin, v1_xmax, v1_ymin, v1_ymax, v1_zmin, v1_zmax);
     }
     else{
     }

     // debug code.
     if (_DBG >= 3){
	  int i, j;
	  for (i = 0; i < 10; i++){
	       printf("%4.2f ", im1[i]);
	  }
	  printf("\n");
	  printf(" ExtractRealData. size: %d, %d, %d %d\n", size[0], size[1], size[2], size[3]);
	  printf("extracted time series:\n");
	  for (i = 0; i < 5; i++){
	       for (j = 0; j < 10; j++){
		    printf("[%d][%d]=%4.2f ", i,j, im1[i*size[3] +j]);
	       }
	       printf("\n");
	  }
     }
     // detrend.
     GPUDetrend(im1, regionSize1, size[3]);
     GPUDetrend(im2, regionSize2, size[3]);

     if (_DBG >= 3){
	  printf("ExtractRealData.cxx: detrended data\n");
	  int i, j;
	  for (i = 0; i < 5; i++){
	       for (j = 0; j < 10; j++){
		    printf("[%d][%d]=%4.2f ", i, j, im1[i*size[3] +j]);
	       }
	       printf("\n");
	  }
     }


}


void ExtractGray(fmriReaderType::Pointer fmriReader,
		maskReaderType::Pointer maskReader,
		PixelType* im1,
		unsigned int* fwdMap,
		unsigned short* bwdMap,
		float grayThresh,
		 int xmin, int xmax,
		 int ymin, int ymax,
		 int zmin, int zmax)


{
     unsigned int x, y, z, t;

     ImageType::IndexType voxIdx;
     MaskImType::IndexType maskIdx;
     ImageType::SizeType size =
          fmriReader->GetOutput()->GetLargestPossibleRegion().GetSize();

     ImageType::Pointer fmriPtr  =  fmriReader->GetOutput();
     MaskImType::Pointer maskPtr  =  maskReader->GetOutput();

     unsigned int nIdx = 0; // linear index.

     for (x = xmin; x <= xmax; x++){
	  for (y = ymin; y <= ymax; y++){
	       for (z = zmin; z <= zmax; z++){
		    maskIdx[0] = x;
		    voxIdx[0] = x;
		    maskIdx[1] = y;
		    voxIdx[1] = y;
		    maskIdx[2] = z;
		    voxIdx[2] = z;
		    if (maskPtr->GetPixel(maskIdx) > grayThresh){
			 for (t = 0; t < size[3]; t++){
			      voxIdx[3] = t;
			      im1[nIdx*size[3] + t] = fmriPtr->GetPixel(voxIdx);
			 }

			 // Forward map. assume row major. a volume.
			 fwdMap[x*size[1]*size[2] + y*size[2] + z] = nIdx;
			 // backward map matrix. Nx3.
			 bwdMap[nIdx*3 + 0] = x;
			 bwdMap[nIdx*3 + 1] = y;
			 bwdMap[nIdx*3 + 2] = z;

			 nIdx++;
		    }
	       }
	  }
     }
     printf("ExtractGray: extracted gray matter voxels: %d\n", nIdx);
}

unsigned int GetRegionSize(maskReaderType::Pointer maskReader,
			   float grayThresh,
			   int xmin, int xmax,
			   int ymin, int ymax,
			   int zmin, int zmax)
{
     unsigned int x, y, z;
     MaskImType::IndexType maskIdx;
     MaskImType::Pointer maskPtr  =  maskReader->GetOutput();

     unsigned int regionSize = 0;
     for (x = xmin; x <= xmax ;  x++){
	  for (y = ymin; y <= ymax; y++){
	       for (z = zmin; z <= zmax; z++){
		    maskIdx[0] = x;
		    maskIdx[1] = y;
		    maskIdx[2] = z;
		    if (maskPtr->GetPixel(maskIdx) > grayThresh){
			 regionSize++;
		    }
	       }
	  }
     }
     if (_DBG >= 1){
	  printf("GetRegionSize: region size %d.\n", regionSize);
     }
     return regionSize;
}
