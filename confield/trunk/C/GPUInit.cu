#include "common.h"

void checkCUDAError(const char *msg);
__global__ void InitKernel(float* C, 
			   const float* im,
			   paraType par,
			   unsigned int M,  // row
			   unsigned int N, // col
			   unsigned short TL); // time series length.

__device__ float GetCorr(const float* im, 
			 unsigned int N, 
			 unsigned int M,  // time series length
			 int idx1, 
			 int idx2);

void GPUInit(float* C, // initial connectivity matrix.
	     const float* im,
	     paraType par,
	     unsigned int regionSize1, 
	     unsigned int regionSize2, 
	     unsigned short TL) // time series length.
{
     unsigned int BN = regionSize1 * (regionSize1+1)/2;
    // pointer to device memory.
     float* gpuC;
     float* gpuIm;
     /* create input and output array on GPU. */
     cudaMalloc((void**) &gpuC, sizeof(float)*regionSize1 * (regionSize1+1)/2);
     checkCUDAError("allocate gpuC");
     cudaMalloc((void**) &gpuIm, sizeof(float)*regionSize1*TL);
     checkCUDAError("CudaSampling, allocate gpuR");


     /* host to device memory. */
     cudaMemcpy(gpuC, C, sizeof(float) * regionSize1 * (regionSize1+1)/2, cudaMemcpyHostToDevice);
     checkCUDAError("Copy C to device");
     cudaMemcpy(gpuIm, im, sizeof(float)*regionSize1*TL, cudaMemcpyHostToDevice);
     checkCUDAError("CudaSampling,  copy im to gpuIm");     

     /* run the kernel function. */
     dim3 dimBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y);
     dim3 dimGrid(ceil((float)regionSize1 / BLOCK_SIZE_X) , 
		  ceil((float)regionSize2 / BLOCK_SIZE_Y));

     InitKernel<<<dimGrid, dimBlock>>>(gpuC, gpuIm, par, regionSize1, regionSize2, TL);
     
     /* Send results back to cpu memeory */
     cudaMemcpy(C, gpuC, sizeof(float)*regionSize1 * (regionSize1+1)/2, cudaMemcpyDeviceToHost);

     /* clean up. */
     cudaFree(gpuC);
     cudaFree(gpuIm);
}

/* Kernel function  */
__global__ void InitKernel(float* C, 
			   const float* im,
			   paraType par,
			   unsigned int M,  // row
			   unsigned int N, // col
			   unsigned short TL) // time series length.
{     

     uint idx1 = blockIdx.x*blockDim.x + threadIdx.x;
     uint idx2 = blockIdx.y*blockDim.y + threadIdx.y;

     float gc1 = 1/(sqrt(2*PI*par.sigma21));
     float gc2 = 1/(sqrt(2*PI*par.sigma22));

     // thread fall outside of matrix C, or mask is zero.
     if (idx1 >= M | idx2 >= N | idx1 > idx2) {return;}

     // compute correlation on the air.
     float yi = GetCorr(im, M, TL, idx1, idx2);

     float lh1 = gc1 * exp(-pow(yi-par.mu1,2)/(2*par.sigma21));
     float lh2 = gc2 * exp(-pow(yi-par.mu2,2)/(2*par.sigma22));     

     if (lh2 > lh1){
	  GTI(C, idx1, idx2, N) = 1;
     }
     else{
	  GTI(C, idx1, idx2, N) = -1;
     }
}


// compute correlation.
__device__ float GetCorr(const float* im, 
			 unsigned int N, 
			 unsigned int M,  // time series length
			 int idx1, 
			 int idx2)
{
     int m = 0;
     float mean1 = 0;
     float mean2 = 0;
     float std1 = 0;
     float std2 = 0;
     int lidx1, lidx2;
     float r = 0; // temp variable for sample correlation.
     if (idx1 == idx2){
	  return DIAGCORR;
     }

     // mean of vectors.
     for (m = 0; m < M; m++){
	  lidx1 = idx1*M + m;
	  lidx2 = idx2*M + m;
	  mean1 = mean1 + *(im + lidx1)/M;
	  mean2 = mean2 + *(im + lidx2)/M;
     }
     /* Standard deviation. */
     for (m = 0; m < M; m++){
	  lidx1 = idx1 * M + m;
	  lidx2 = idx2 * M + m;
	  std1 = std1 + pow((*(im+lidx1) - mean1), 2)/(M-1);
	  std2 = std2 + pow((*(im+lidx2) - mean2), 2)/(M-1);
     }
     std1 = sqrt(std1);
     std2 = sqrt(std2);
     /* Sample Correlation. */
     if (std1 == 0 | std2 == 0){
	  return 0;
     }
     else{
	  for (m = 0; m < M; m++){
	       lidx1 = idx1 * M + m;
	       lidx2 = idx2 * M + m;
	       r = r + (*(im + lidx1) - mean1) * (*(im + lidx2)-mean2)/((M-1)*std1*std2);
	  }
	  return r;
     }
}

