#include "cppheader.h"
#include "common.h"

void ExtractToyData(float* & im1,
		    float*  & im2,
		    unsigned int*  & fwdMap1,
		    unsigned int*  & fwdMap2,
		    unsigned short*  & bwdMap1,
		    unsigned short*  & bwdMap2,
		    unsigned int& regionSize1,
		    unsigned int& regionSize2,
		    unsigned short* size,
		    float snr)
{
     const float F = 0.2; // freq of sine wave.
     const float mu = 800;
     const float fs = 1;
     const unsigned short D = 300; // # of time points.
     const unsigned short N = 100;
     const unsigned short A = 20; // amplitude.
     const float ts = 1/fs;
     const float sigma2 = A/snr;
     unsigned short* label = (unsigned short*)calloc(N, sizeof(unsigned short));
     unsigned int M = 0;
     int i, j, k;
     // allocate memory for im1, im2, etc. will be freed at upper level.
     im1 = (float*) malloc(N*D*sizeof(float));
     im2 = (float*) malloc(N*D*sizeof(float));
     // For this toy data, we also make a fwdmap volume with size 1 x 1 x N.
     size[0] = 1;
     size[1] = 1;
     size[2] = 100;
     size[3] = D;
     regionSize1 = N;
     regionSize2 = N;
     fwdMap1 = (unsigned int*) 
	  calloc(size[0]*size[1]*size[2], sizeof(unsigned int));
     fwdMap2 = (unsigned int*) 
	  calloc(size[0]*size[1]*size[2], sizeof(unsigned int));
     bwdMap1 = (unsigned short*) 
	  malloc(regionSize1 * 3 *  sizeof(unsigned short));
     bwdMap2 = (unsigned short*) 
	  malloc(regionSize2 * 3 * sizeof(unsigned short));

     // boost random number generator
//     base_generator_type generator(time(NULL));
     base_generator_type generator(2010);
     gen_type normalGenerator(generator, distribution_type(0, 1));
     
     for (i = 5; i < 15; i++){
	  label[i] = 1;
     }
     for (i = 70; i < 80; i++){
	  label[i] = 1;
     }
     for (i = 0; i < N; i++){
	  if (label[i]) M++;
     }

     if (_DBG >= 3){
	  printf("now generating im1 and im2.\n");
     } 
     for (i = 0; i < N; i++){
	  for (j = 0; j < D; j++){
	       if (label[i]){
		    im1[i*D+j] = mu + 
			 sigma2 * normalGenerator() + A * sin(2*PI*F*j*ts);
		    im2[i*D+j] = mu + 
			 sigma2 * normalGenerator() + A * sin(2*PI*F*j*ts);
	       }
	       else{
		    im1[i*D+j] = mu + sigma2 * normalGenerator();
		    im2[i*D+j] = mu + sigma2 * normalGenerator();
	       }
	  }
     }
      
     // make fwdMap matrix. identity mapping from (i,j,k) --> n = k
     if (_DBG >= 2){
	  printf("ExtractToyData, now make fwdmap.\n");
     }
     for (k = 0; k < size[2]; k++){
	  fwdMap1[0*size[1]*size[2] + 0*size[2] + k] = k;
	  fwdMap2[0*size[1]*size[2] + 0*size[2] + k] = k;
     }
     // bamke bwdMap matrix.
     for (k = 0; k < size[2]; k++){
	  // backward map matrix. Nx3.
	  bwdMap1[k*3 + 0] = 0;
	  bwdMap1[k*3 + 1] = 0;
	  bwdMap1[k*3 + 2] = k;

	  bwdMap2[k*3 + 0] = 0;
	  bwdMap2[k*3 + 1] = 0;
	  bwdMap2[k*3 + 2] = k;
     }
	  
     free(label);
}

