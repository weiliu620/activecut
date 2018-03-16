#include "cppheader.h"
#include "common.h"

int main(int argc, char* argv[])
{

     unsigned int i, j, t, n;
     // Two 1-D array, storing reshaped voxel data.
     PixelType* im1 = NULL;
     PixelType* im2 = NULL;
     // slice id to extract.
     float snr;

     // mapping matrix.
     unsigned int* fwdMap1 = NULL;
     unsigned int* fwdMap2 = NULL;
     unsigned short* bwdMap1 = NULL;
     unsigned short* bwdMap2 = NULL;
     float* R = NULL;
     float* C = NULL;
     float* post2 = NULL;
     // number of voxels in two regions.
     unsigned int regionSize1;
     unsigned int regionSize2;
     // 4d array for x, y, z, t dim size.
     unsigned short size[4] = {0, 0, 0, 0};

     ConfigFile config(PARAFILE);
     int toy = 0;
     config.readInto(toy, "toy");
     if(toy){
	  if (_DBG >= 1){
	       printf("Generating synthetic data...\n");
	  }
	  config.readInto(snr, "snr");
	  ExtractToyData(im1, im2, fwdMap1, fwdMap2, bwdMap1, bwdMap2,
			 regionSize1, regionSize2, size, snr);
     }
     else{
	  // compute gray matter voxel number.
	  if (config.read("graymattervox", 0)) {
	       graymatter(config);
	  }
	  
	  // im1, im2, fwdMap1...get allocated here. Need free at end of main().
	  ExtractRealData(im1, im2, fwdMap1, fwdMap2, bwdMap1, bwdMap2,
			regionSize1, regionSize2, size, config);
     }
     // debug output.
     unsigned int tsize = 5;
     if(_DBG >= 3){
	  printf("main. size: %d, %d, %d %d\n", size[0], size[1], size[2], size[3]);	 
	  if (im1 != NULL){
	       printf("\nim1:\n");
	       for (i = 0; i < tsize; i++){
		    for (j = 0; j < 10; j++){
			 printf("%4.2f ", im1[i*size[3] +j]);
		    }
		    printf("\n");
	       }
	  }
	  else{
	       fprintf(stderr, "im1 is null!\n");
	  }

	  if (im2 != NULL){
	       printf("main. im2:\n");
	       for (i = 0; i < tsize; i++){
		    for (j = 0; j < 10; j++){
			 printf("%4.2f ", im2[i*size[3] +j]);
		    }
		    printf("\n");
	       }
	  }
	  else{
	       fprintf(stderr, "im2 is null!\n");
	  }

	  printf("fwdMap1 and bwdMap1:\n");
	  int xx, yy, zz;
	  for (i = 0; i < tsize; i++){
	       printf("%d, %d, %d.   ",  bwdMap1[i*3 + 0], bwdMap1[i*3+1], bwdMap1[i*3+2]);
	       xx = bwdMap1[i*3+0];
	       yy = bwdMap1[i*3+1];
	       zz = bwdMap1[i*3+2];
	       printf("fwdMap1: %d\n", fwdMap1[xx*size[1]*size[2]+yy*size[2] + zz]);
	  }
	  printf("fwdMap2 and bwdMap2:\n");
	  for (i = 0; i < tsize; i++){
	       printf("%d, %d, %d.   ",  bwdMap2[i*3 + 0], bwdMap2[i*3+1], bwdMap2[i*3+2]);
	       xx = bwdMap2[i*3+0];
	       yy = bwdMap2[i*3+1];
	       zz = bwdMap2[i*3+2];
	       printf("fwdMap1: %d\n", fwdMap2[xx*size[1]*size[2]+yy*size[2] + zz]);
	  }

     }

     // allocate correlation matrix. Row major.
     R = (float*) malloc(regionSize1 * regionSize2 * sizeof(float));
     C = (float*) malloc(regionSize1 * regionSize2 * sizeof(float));
     post2 = (float*) calloc(regionSize1*regionSize2, sizeof(float));
     // Compute correlation. 
     GPUCorr(R, im1, im2, size[3], regionSize1, regionSize2);
     if (_DBG >= 4){
	  printf("R:\n");
	  for (i = RMIN; i <= RMAX; i++){
	       for (j = CMIN; j < CMAX; j++){
		    printf("[%d][%d]=%1.2f ", i, j, R[i*regionSize2 + j]);
	       }
	       printf("\n");
	  }
	  FieldViewer(R, regionSize1, regionSize2, fwdMap1, fwdMap2, "R");
     }

     if (THCORRMAP == 1){
	  ThreshCorr(R, post2, regionSize1, regionSize2, fwdMap1, fwdMap2, size, config);
	  SaveResults(post2, regionSize1, regionSize2, fwdMap1, fwdMap2, config);
	  return (0);
     }

     // Compute Posterior.
     unsigned short volSize[3];
     volSize[0] = size[0];
     volSize[1] = size[1];
     volSize[2] = size[2];
     printf("Original fmri volume size: %dx%dx%d\n", size[0], size[1], size[2]);
     EMPost(C, R, post2, fwdMap1, fwdMap2, bwdMap1, bwdMap2,
	    regionSize1, regionSize2, volSize);
     printf("main. EMPost done.\n");
     
     // save results.
     SaveResults(post2, regionSize1, regionSize2, fwdMap1, fwdMap2, config);
     
     // clean up.
     free(im1);
     free(im2);
     free(fwdMap1);
     free(fwdMap2);
     free(bwdMap1);
     free(bwdMap2);
     free(R);
     free(C);
     free(post2);

}
