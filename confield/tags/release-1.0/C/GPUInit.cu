#include "common.h"

void checkCUDAError(const char *msg);
__global__ void InitKernel(float* C, 
			   const float* R, 
			   paraType par,
			   unsigned int M,  // row
			   unsigned int N); // col

void GPUInit(float* C, // initial connectivity matrix.
	     const float* R, // correlation (or other affinity) matrix, from data.
	     paraType par,
	     unsigned int regionSize1, 
	     unsigned int regionSize2)
{
    // pointer to device memory.
     float* gpuC;
     float* gpuR;

     /* create input and output array on GPU. */
     cudaMalloc((void**) &gpuC, sizeof(float)*regionSize1*regionSize2);
     checkCUDAError("allocate gpuC");
     cudaMalloc((void**) &gpuR, sizeof(float)*regionSize1*regionSize2);
     checkCUDAError("Allocate gpuR");

     /* host to device memory. */
     cudaMemcpy(gpuR, R, sizeof(float)*regionSize1*regionSize2, cudaMemcpyHostToDevice);
     checkCUDAError("copy R to device");

     cudaMemcpy(gpuC, C, sizeof(float)*regionSize1*regionSize2, cudaMemcpyHostToDevice);
     checkCUDAError("Copy C to device");


     /* run the kernel function. */
     dim3 dimBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y);
     dim3 dimGrid(ceil((float)regionSize1 / BLOCK_SIZE_X) , 
		  ceil((float)regionSize2 / BLOCK_SIZE_Y));

     InitKernel<<<dimGrid, dimBlock>>>(gpuC, gpuR, par, regionSize1, regionSize2);
     
     /* Send results back to cpu memeory */
     cudaMemcpy(C, gpuC, sizeof(float)*regionSize1*regionSize2, cudaMemcpyDeviceToHost);

     /* clean up. */
     cudaFree(gpuC);
     cudaFree(gpuR);

}

/* Kernel function  */
__global__ void InitKernel(float* C, 
			   const float* R, 
			   paraType par,
			   unsigned int M,  // row
			   unsigned int N) // col
{     

     uint idx1 = blockIdx.x*blockDim.x + threadIdx.x;
     uint idx2 = blockIdx.y*blockDim.y + threadIdx.y;

     float gc1 = 1/(sqrt(2*PI*par.sigma21));
     float gc2 = 1/(sqrt(2*PI*par.sigma22));

     // thread fall outside of matrix C, or mask is zero.
     if (idx1 >= M | idx2 >= N) {return;}

     float lh1 = gc1 * exp(-pow(R[idx1*N + idx2]-par.mu1,2)/(2*par.sigma21));
     float lh2 = gc2 * exp(-pow(R[idx1*N + idx2]-par.mu2,2)/(2*par.sigma22));     

     if (lh2 > lh1){
	  C[idx1*N + idx2] =  1;
     }
     else{
	  C[idx1*N + idx2] =  -1;
     }

}
