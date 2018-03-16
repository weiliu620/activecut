#include "common.h"

__global__ void CorrKernel(const float* gpuIm1, const float* gpuIm2, float* R, int N1, int N2, int M);
__device__ float GetCorr(const float* im, 
			 unsigned int N, 
			 unsigned int M,  // time series length
			 int idx1, 
			 int idx2);

void checkCUDAError(const char *msg);

void GPUCorr(float* R, const float* im1, const float* im2, int M, int N1, int N2)
{
    // pointer to device memory.
     float* gpuR;
     float* gpuIm1;
     float* gpuIm2;
     
     printf("GPUCorr: N1(voxel num) = %d, N2 = %d, M (time series length)= %d\n.", 
	       N1, N2, M);
          
     /* create input and output array on GPU. */
     cudaMalloc((void**) &gpuR, sizeof(float)*N1*(N1+1)/2);
     checkCUDAError("Allocate device gpuR");
     cudaMalloc((void**) &gpuIm1, sizeof(float)*N1*M);
     checkCUDAError("Allocate device gpuIm1");
     cudaMalloc((void**) &gpuIm2, sizeof(float)*N2*M);
     checkCUDAError("Allocate device, gpuIm2.");

     /* host to device memory. */
     cudaMemcpy(gpuR, R, sizeof(float)*N1*(N1+1)/2, cudaMemcpyHostToDevice);
     cudaMemcpy(gpuIm1, im1, sizeof(float)*N1*M, cudaMemcpyHostToDevice);
     cudaMemcpy(gpuIm2, im2, sizeof(float)*N2*M, cudaMemcpyHostToDevice);
     checkCUDAError("memory copy from Host to Device");

     /* run the kernel function. */
     int blockSize = 16;
     int gridDimx = N1/blockSize + (N1%blockSize == 0?0:1);
     int gridDimy = N2/blockSize + (N2%blockSize == 0?0:1);

     printf("GPUCorr, blockSize: %dx%d, gridDim = %dx%d\n", blockSize, 
	    blockSize, gridDimx, gridDimy);
     dim3 dimBlock(blockSize, blockSize);
     dim3 dimGrid(gridDimx, gridDimy);
     CorrKernel<<<dimGrid, dimBlock>>>(gpuIm1, gpuIm2, gpuR, N1, N2, M);
     checkCUDAError(" GPUCorr, main, call Kernel.");

     /* Send results back to cpu memeory */
     cudaMemcpy(R, gpuR, sizeof(float)*N1*(N1+1)/2, cudaMemcpyDeviceToHost);
     checkCUDAError("GPUCorr, main, memcopy from Device to Host.");
     
     /* clean up. */
     cudaFree(gpuR);
     cudaFree(gpuIm1);
     cudaFree(gpuIm2);
}

/* Kernel function  */
__global__ void CorrKernel(const float* gpuIm1, const float* gpuIm2, float* R, int N1, int N2, int M)
{
     int m = 0;
     float mean1 = 0;
     float mean2 = 0;
     float std1 = 0;
     float std2 = 0;
     int lidx1, lidx2;
     float r = 0; // temp variable for sample correlation.
     int idx1 = blockIdx.x * blockDim.x + threadIdx.x;
     int idx2 = blockIdx.y * blockDim.y + threadIdx.y;
     if ((idx1 >= N1) | (idx2 >= N2) | (idx1 > idx2)){ return; }
     if (idx1 == idx2){
	  GTI(R, idx1, idx2, N1) = DIAGCORR;
	  return;
     }

     // mean of vectors.
     for (m = 0; m < M; m++){
	  lidx1 = idx1*M + m;
	  lidx2 = idx2*M + m;
	  mean1 = mean1 + *(gpuIm1 + lidx1)/M;
	  mean2 = mean2 + *(gpuIm2 + lidx2)/M;
     }
     /* Standard deviation. */
     for (m = 0; m < M; m++){
	  lidx1 = idx1 * M + m;
	  lidx2 = idx2 * M + m;
	  std1 = std1 + pow((*(gpuIm1+lidx1) - mean1), 2)/(M-1);
	  std2 = std2 + pow((*(gpuIm2+lidx2) - mean2), 2)/(M-1);
     }
     std1 = sqrt(std1);
     std2 = sqrt(std2);
     /* Sample Correlation. */
     if (std1 == 0 | std2 == 0){
	  GTI(R, idx1, idx2, N1) = 0;
     }
     else{
	  for (m = 0; m < M; m++){
	       lidx1 = idx1 * M + m;
	       lidx2 = idx2 * M + m;
	       r = r + (*(gpuIm1+lidx1) - mean1) * (*(gpuIm2 + lidx2)-mean2)/((M-1)*std1*std2);
	  }
	  GTI(R, idx1, idx2, N1) = r;
     }
}
     
void checkCUDAError(const char *msg)
{
     cudaError_t err = cudaGetLastError();
     if( cudaSuccess != err) 
     {
	  fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
	  //exit(-1);
     }                         
}
