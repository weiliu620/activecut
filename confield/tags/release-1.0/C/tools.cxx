#include "common.h"
#include "cppheader.h"

int SaveResults(const float* post,
		unsigned short N1,
		unsigned short N2,
		const unsigned int* fwdMap1,
		const unsigned int* fwdMap2,
		ConfigFile config)
		
{
/*
     char grayMaskFile[100];
     char postFile[100];
     strcpy(grayMaskFile, datadir);
     strcpy(postFile, datadir);
     strcat(grayMaskFile, "/mask.nii");
     strcat(postFile, "/post.nii");
*/

     // seed voxel.
     unsigned short seed[3] = {0, 0, 0};

     string datadir;
     config.readInto(datadir, "datadir");
     string grayMaskFile(datadir);
     string postFile(datadir);
     grayMaskFile.append("/mask.nii");
     postFile.append("/post.nii");

     // reader.
     maskReaderType::Pointer maskReader = maskReaderType::New();
     WriterType::Pointer postWriter = WriterType::New();
     MaskImType::Pointer imPointer = maskReader->GetOutput();

     postWriter->SetInput(imPointer);

     maskReader->SetFileName(grayMaskFile);
     postWriter->SetFileName(postFile);
     // update reader to get image data.
     maskReader->Update();
     int i = 0, j = 0, k = 0;
     unsigned short n1;
     unsigned short n2;

     MaskImType::IndexType pixelIndex;
     
     MaskImType::SizeType size =  imPointer->GetLargestPossibleRegion().GetSize();

     printf("SaveResults, size = [%d][%d][%d].\n", size[0], size[1], size[2]);
     imPointer->FillBuffer(0);
     postWriter->Update();

     config.readInto(seed[0], "seed0");
     config.readInto(seed[1], "seed1");
     config.readInto(seed[2], "seed2");
     n1 = fwdMap1[seed[0] * (size[1]*size[2])
		  + seed[1] * size[2] + seed[2]];

     unsigned short s1 = 0, s2 = 0;
     unsigned short v1_xmin, v1_xmax, v1_ymin, v1_ymax, v1_zmin, v1_zmax;
     unsigned short v2_xmin, v2_xmax, v2_ymin, v2_ymax, v2_zmin, v2_zmax;

     unsigned short sliceField, volField;
     config.readInto(sliceField, "sliceField");
     config.readInto(volField, "volField");
     if (sliceField){
	  config.readInto(s1, "s1");
	  config.readInto(s2, "s2");

	  for (i = 0; i < size[0]; i++){
	       for (j = 0; j < size[1]; j++){	       
		    pixelIndex[0] = i;
		    pixelIndex[1] = j;
		    pixelIndex[2] = s2;
		    n2 = fwdMap2[i * size[1]*size[2] + j * size[2] + s2];
		    if (n2 != 0){
			 imPointer->SetPixel(pixelIndex, post[n1*N2 + n2]);
		    }
	       }
	  }
	  postWriter->Update();

     }
     else if (volField){
     
	  config.readInto(v1_xmin, "v1_xmin");
	  config.readInto(v1_xmax, "v1_xmax");
	  config.readInto(v1_ymin, "v1_ymin");
	  config.readInto(v1_ymax, "v1_ymax");
	  config.readInto(v1_zmin, "v1_zmin");
	  config.readInto(v1_zmax, "v1_zmax");

	  config.readInto(v2_xmin, "v2_xmin");
	  config.readInto(v2_xmax, "v2_xmax");
	  config.readInto(v2_ymin, "v2_ymin");
	  config.readInto(v2_ymax, "v2_ymax");
	  config.readInto(v2_zmin, "v2_zmin");
	  config.readInto(v2_zmax, "v2_zmax");


	  for (i = v2_xmin; i <= v2_xmax; i++){
	       for (j = v2_ymin; j <= v2_ymax; j++){
		    for (k = v2_zmin; k <= v2_zmax; k++){
			 pixelIndex[0] = i;
			 pixelIndex[1] = j;
			 pixelIndex[2] = k;
			 n2 = fwdMap2[i * size[1]*size[2] + j * size[2] + k];
			 if (n2 != 0){
			 imPointer->SetPixel(pixelIndex, post[n1*N2 + n2]);
			 }
		    }
	       }
	  }
	  postWriter->Update();
     }
     else{
	  fprintf(stderr, "either sliceField or volField have to be set. Check parameter file.");
	  exit(1);
     }
}


void MLEstep(float* R,
	     paraType par,
	     unsigned int N1,
	     unsigned int N2,
	     unsigned int Nk1,
	     unsigned int Nk2,
	     float* post2)

{
     unsigned int n;
     float gc1 = 1/sqrt(2*PI*par.sigma21);
     float gc2 = 1/sqrt(2*PI*par.sigma22);
     
     float lhood1 = 0;
     float lhood2 = 0;
     float prsum = 0;
     float lhoodsum = 0;
     float postSum = 0;
     float tpost1 = 0;
     float tpost2 = 0;
     for (n = 0; n < N1*N2; n++){
	  lhood1 = gc1 * exp(- (R[n]-par.mu1)*(R[n]-par.mu1)/(2*par.sigma21));
	  lhood2 = gc2 * exp(- (R[n]-par.mu2)*(R[n]-par.mu2)/(2*par.sigma22));
	  lhoodsum = lhood1 + lhood2;
	  lhood1 = lhood1/lhoodsum;
	  lhood2 = lhood2/lhoodsum;
	  
	  prsum = Nk1 + Nk2;
	  tpost1 = (Nk1/prsum) * lhood1;
	  tpost2 = (Nk2/prsum) * lhood2;
	  postSum = tpost1 + tpost2;
	  tpost1 = tpost1/postSum;
	  tpost2 = tpost2/postSum;

	  post2[n] = tpost2;
     }
}

int ThreshCorr(const float*R, 
	       float* post,
	       unsigned int N1,
	       unsigned int N2,
	       const unsigned int* fwdMap1,
	       const unsigned int* fwdMap2,
	       unsigned short* size,
	       ConfigFile config)

{
     unsigned int n;
     unsigned int n1, n2;
     int i, j, k;
     // seed voxel.
     unsigned short seed[3] = {0, 0, 0};

     printf("size is [%i][%i][%i]\n", size[0], size[1], size[2]);
     printf("N1 = %d, N2 = %d\n", N1, N2);

     for ( n = 0; n < N1*N2; n++){
	  if (R[n] > 0){
	       post[n] = R[n];
	  }
	  else{
	       post[n] = 0;
	  }
     }
}
