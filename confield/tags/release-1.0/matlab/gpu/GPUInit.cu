#include "mex.h"
#include "cuda.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <thrust/version.h>
// Forward Declaration.
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
#define PI 3.141592654
#define MYTEST 0
void GPUmain(float* C, // initial connectivity matrix.
	     const float* R, // correlation (or other affinity) matrix, from data.
	     float mu1, 
	     float mu2,
	     float sigma21,
	     float sigma22,
	     uint M, 
	     uint N,
	     const char* CMask);

__global__ void InitKernel(float* C, 
			   const float* R, 
			   float mu1, 
			   float mu2, 
			   float sigma21, 
			   float sigma22,
			   uint M,
			   uint N,
			   const char* CMask);
void checkCUDAError(const char *msg);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
/* Check for proper number of input and output arguments */
     if (nrhs != 6) {
	  mexErrMsgTxt("Arguments should be 5.\n");
     }
     if (nlhs > 1){
	  mexErrMsgTxt("Too many output arguments.");
     }

     float* R = (float*)mxGetPr(prhs[0]); 
     float mu1 = mxGetScalar(prhs[1]);
     float mu2 = mxGetScalar(prhs[2]);
     float sigma21 = mxGetScalar(prhs[3]);
     float sigma22 = mxGetScalar(prhs[4]);
     char* CMask = (char*)mxGetPr(prhs[5]);

     mwSize M = mxGetM(prhs[0]);
     mwSize N = mxGetN(prhs[0]);

     // Create output data.
     plhs[0] = mxCreateNumericMatrix(M, N, mxSINGLE_CLASS, mxREAL);
     float* C = (float*) mxGetPr(plhs[0]);

     GPUmain(C, R, mu1, mu2, sigma21, sigma22, (uint)M, (uint)N, CMask);
}

void GPUmain(float* C, // initial connectivity matrix.
	     const float* R, // correlation (or other affinity) matrix, from data.
	     float mu1, 
	     float mu2,
	     float sigma21,
	     float sigma22,
	     uint M, 
	     uint N,
	     const char* CMask)
{
    // pointer to device memory.
     float* gpuC;
     float* gpuR;
     char* gpuCMask;

     /* create input and output array on GPU. */
     cudaMalloc((void**) &gpuC, sizeof(float)*M*N);
     checkCUDAError("allocate gpuC");
     cudaMalloc((void**) &gpuR, sizeof(float)*M*N);
     checkCUDAError("Allocate gpuR");
     cudaMalloc((void**) &gpuCMask, sizeof(char)*M*N);
     checkCUDAError("Allocate device memory gpuCMask.");

     /* host to device memory. */
     cudaMemcpy(gpuR, R, sizeof(float)*M*N, cudaMemcpyHostToDevice);
     checkCUDAError("copy R to device");

     cudaMemcpy(gpuC, C, sizeof(float)*M*N, cudaMemcpyHostToDevice);
     checkCUDAError("Copy C to device");
     /* test */
     mexPrintf("uint size: %d\n", sizeof(uint));

     cudaMemcpy(gpuCMask, CMask, sizeof(char)*M*N, cudaMemcpyHostToDevice);
     checkCUDAError("Copy CMask to device");

     /* run the kernel function. */
     dim3 dimBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y);
     dim3 dimGrid(ceil((float)M / BLOCK_SIZE_X) , ceil((float)N / BLOCK_SIZE_Y));

     InitKernel<<<dimGrid, dimBlock>>>(gpuC, gpuR, mu1, mu2, sigma21, sigma22, M, N, CMask);
     
     /* Send results back to cpu memeory */
     cudaMemcpy(C, gpuC, sizeof(float)*M*N, cudaMemcpyDeviceToHost);

     /* clean up. */
     cudaFree(gpuC);
     cudaFree(gpuR);
     cudaFree(gpuCMask);
}

/* Kernel function  */
__global__ void InitKernel(float* C, 
			   const float* R, 
			   float mu1, 
			   float mu2, 
			   float sigma21, 
			   float sigma22,
			   uint M,
			   uint N,
			   const char* CMask)
{     

     uint idx1 = blockIdx.x*blockDim.x + threadIdx.x;
     uint idx2 = blockIdx.y*blockDim.y + threadIdx.y;

     float gc1 = 1/(sqrt(2*PI*sigma21));
     float gc2 = 1/(sqrt(2*PI*sigma22));

     // thread fall outside of matrix C, or mask is zero.
     if (idx1 >= M | idx2 >= N | *(CMask + idx2*M+idx1) == 0) {return;}

     float lh1 = gc1 * exp(-pow(R[idx2*M + idx1]-mu1,2)/(2*sigma21));
     float lh2 = gc2 * exp(-pow(R[idx1*M + idx1]-mu2,2)/(2*sigma22));     

     if (lh2 > lh1){
	  C[idx2*M + idx1] =  1;
     }
     else{
	  C[idx2*M + idx1] =  -1;
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
