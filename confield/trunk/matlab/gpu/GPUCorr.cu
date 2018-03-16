#include "mex.h"
#include "cuda.h"
#include <stdio.h>
// Forward Declaration.
void GPUmain(float* R, const float* im1, const float* im2, int M, int N1, int N2);
__global__ void CorrKernel(const float* gpuIm1, const float* gpuIm2, float* R, int N1, int N2, int M);
void checkCUDAError(const char *msg);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
/* Check for proper number of input and output arguments */
     if (nrhs != 2) {
	  mexErrMsgTxt("Two input argument required: matrix");
     }
     if (nlhs > 1){
	  mexErrMsgTxt("Too many output arguments.");
     }
     /* input and output data pointer. */
     float* im1 = (float*)mxGetPr(prhs[0]);
     float* im2 = (float*)mxGetPr(prhs[1]); 

 
     mwSize M = mxGetM(prhs[0]);
     mwSize N1 = mxGetN(prhs[0]);
     mwSize N2 = mxGetN(prhs[1]);
     // Create output data.
     plhs[0] = mxCreateNumericMatrix(N1, N2, mxSINGLE_CLASS, mxREAL);
     /* correlation (or other affinity) matrix, from data. */
     float* R = (float*)mxGetPr(plhs[0]); 

     GPUmain(R, im1, im2, M, N1, N2);

}

void GPUmain(float* R, const float* im1, const float* im2, int M, int N1, int N2)
{
    // pointer to device memory.
     float* gpuR;
     float* gpuIm1;
     float* gpuIm2;
     
     mexPrintf("in GPUmain: N1 = %d, N2 = %d, M = %d\n.", N1, N2, M);
          
     /* create input and output array on GPU. */
     cudaMalloc((void**) &gpuR, sizeof(float)*N1*N2);
     checkCUDAError("Allocate device gpuR");
     cudaMalloc((void**) &gpuIm1, sizeof(float)*M*N1);
     checkCUDAError("Allocate device gpuIm1");
     cudaMalloc((void**) &gpuIm2, sizeof(float)*M*N2);
     checkCUDAError("Allocate device, gpuIm2.");

     /* host to device memory. */
     cudaMemcpy(gpuR, R, sizeof(float)*N1*N2, cudaMemcpyHostToDevice);
     cudaMemcpy(gpuIm1, im1, sizeof(float)*M*N1, cudaMemcpyHostToDevice);
     cudaMemcpy(gpuIm2, im2, sizeof(float)*M*N2, cudaMemcpyHostToDevice);
     checkCUDAError("memory copy from Host to Device");

     /* run the kernel function. */
     int blockSize = 16;
     int gridDimx = N1/blockSize + (N1%blockSize == 0?0:1);
     int gridDimy = N2/blockSize + (N2%blockSize == 0?0:1);

     mexPrintf("blockSize: %dx%d, gridDim = %dx%d\n", blockSize, blockSize, gridDimx, gridDimy);
     dim3 dimBlock(blockSize, blockSize);
     dim3 dimGrid(gridDimx, gridDimy);
     CorrKernel<<<dimGrid, dimBlock>>>(gpuIm1, gpuIm2, gpuR, N1, N2, M);
     checkCUDAError(" GPUCorr, main, call Kernel.");

     /* Send results back to cpu memeory */
     cudaMemcpy(R, gpuR, sizeof(float)*N1*N2, cudaMemcpyDeviceToHost);
     checkCUDAError("GPUCorr, main, memcopy from Device to Host.");
/*
     for (int n = M; n < M+M; n++)
     {
	  mexPrintf("im1[%d] = %f, im2[%d] = %f\n",n,  *(im1+n), n, *(im2+n));
     }    
*/

     
     mexPrintf("R[1][1] = %f\n", *(R + N1*1+1));

     /* clean up. */
     cudaFree(gpuR);
     cudaFree(gpuIm1);
     cudaFree(gpuIm2);
}

/* Kernel function  */
__global__ void CorrKernel(const float* gpuIm1, const float* gpuIm2, float* R, int N1, int N2, int M)
{
/*
     int idx = blockIdx.x*blockDim.x+threadIdx.x;
     int idx1 = idx%N1;
     int idx2 = (idx - idx1)/N1;
*/
     int idx1 = blockIdx.x * blockDim.x + threadIdx.x;
     int idx2 = blockIdx.y * blockDim.y + threadIdx.y;
     if ((idx1 > N1) | (idx2 > N2))
     {
	  return;
     }
     else
     {
	  int m = 0;
	  float mean1 = 0;
	  float mean2 = 0;
	  float std1 = 0;
	  float std2 = 0;
	  int lidx1, lidx2;
	  float r = 0; // temp variable for sample correlation.

	  // mean of first vector.
	  for (m = 0; m < M; m++)
	  {
	       lidx1 = idx1*M + m;
	       lidx2 = idx2*M + m;
	       mean1 = mean1 + *(gpuIm1 + lidx1)/M;
	       mean2 = mean2 + *(gpuIm2 + lidx2)/M;
	  }

	  /* Standard deviation. */
	  for (m = 0; m < M; m++)
	  {
	       lidx1 = idx1 * M + m;
	       lidx2 = idx2 * M + m;
	       std1 = std1 + pow((*(gpuIm1+lidx1) - mean1), 2)/(M-1);
	       std2 = std2 + pow((*(gpuIm2+lidx2) - mean2), 2)/(M-1);
	  }
	  std1 = sqrt(std1);
	  std2 = sqrt(std2);
	  /* Sample Correlation. */
	  if (std1 == 0 | std2 == 0)
	  {
	       *(R + idx2*N1 + idx1) = 0;
	  }
	  else
	  {
	       for (m = 0; m < M; m++)
	       {
		    lidx1 = idx1 * M + m;
		    lidx2 = idx2 * M + m;
		    r = r + (*(gpuIm1+lidx1) - mean1) * (*(gpuIm2 + lidx2)-mean2)/((M-1)*std1*std2);
	       }
	       *(R+idx2*N1 + idx1) = r;
	  }

	  // *(R + idx2*N1 + idx1) = mean1;
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
