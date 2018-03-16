//#undef _GLIBCXX_ATOMIC_BUILTINS
#include "mex.h"
#include "cuda.h"
#include <math.h>

__constant__ float constData;
__global__ void square_elements(float* in, float* out, int N, float* randMat);
void base(const float* inputPtr, float* outputPtr, mwSize nRow, mwSize nCol);
__device__ float  myrand(int xz);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
/* Check for proper number of input and output arguments */
     if (nrhs != 1) {
	  mexErrMsgTxt("One input argument required.");
     }
     if (nlhs > 1){
	  mexErrMsgTxt("Too many output arguments.");
     }
/* Check data type of input argument */
     if (!(mxIsSingle(prhs[0]))) {
	  mexErrMsgTxt("Input array must be of type single.");
     }

     // from Matlab.
     float* inputPtr;
     float* outputPtr;
 
     mwSize M = mxGetM(prhs[0]);
     mwSize N = mxGetN(prhs[0]);
  

     // Create output data.
     plhs[0] = mxCreateNumericMatrix(M, N, mxSINGLE_CLASS, mxREAL);
     inputPtr = (float*) mxGetPr(prhs[0]);
     outputPtr = (float*) mxGetPr(plhs[0]);
     // enter cuda function.
     base(inputPtr, outputPtr, M, N);
}

void base(const float* inputPtr, float* outputPtr, mwSize M, mwSize N)
{
   /* generate random data by calling Matlab rand function. */
     mxArray* rprhs[3];
     mxArray* rplhs[1];
     rprhs[0] = mxCreateDoubleScalar(M);
     rprhs[1] = mxCreateDoubleScalar(N);
     rprhs[2] = mxCreateString("single");
     mexPrintf("test, M = %d, N = %d\n", M, N);
     mexCallMATLAB(1, rplhs, 3, rprhs, "rand");
     /* get pointer to random matrix. */
     float* randMat = (float*) mxGetPr(rplhs[0]);
     mwSize rM= mxGetM(rplhs[0]);
     mwSize rN= mxGetN(rplhs[0]);
     mexPrintf("output random matrix, M = %d, N = %d\n", rM, rN);
     /* test random matrix. */
     for (int n = 0; n < 10; n++){
	  mexPrintf("randMat[%d] = %f\n", n, *(randMat + n));
     }


    // for GPU
     float* inputGPUPtr;
     float* outputGPUPtr;

     /* create input and output array on GPU. */
     cudaMalloc((void**) &inputGPUPtr, sizeof(float)*M*N);
     cudaMalloc((void**) &outputGPUPtr, sizeof(float)*M*N);

/* The input array is single precision, it can be sent directly to the
   card */
     cudaMemcpy(inputGPUPtr, inputPtr, sizeof(float)*M*N, cudaMemcpyHostToDevice);

     /* test on constant memory. */

     float mydata[256];
     for (int n = 0; n < 256; n++)
     {
	  mydata[n] = 10;
     }

     cudaMemcpyToSymbol(constData, mydata, sizeof(float)); 
     
     /* run the kernel function. */
     int blockSize = 512;
     int nBlocks = (M*N)/blockSize + ((M*N)%blockSize == 0?0:1);
     mexPrintf("blockSize: %d, nBlocks = %d\n", blockSize, nBlocks);
     dim3 dimBlock(blockSize);
     dim3 dimGrid(nBlocks);
     square_elements<<<dimGrid, dimBlock>>>(inputGPUPtr, outputGPUPtr, M*N, randMat);
     
     /* Send results back to cpu memeory */
     cudaMemcpy(outputPtr, outputGPUPtr, sizeof(float)*M*N, cudaMemcpyDeviceToHost);

     /* clean up. */
     cudaFree(inputGPUPtr);
     cudaFree(outputGPUPtr);

     /* Scratch pad.  */
     unsigned long long bigint;
     mexPrintf("biging size: %d, %d\n", sizeof(bigint), sizeof(unsigned long long));



}


/* Kernel to square elements of the array on the GPU */
__global__ void square_elements(float* in, float* out, int N, float* randMat)
{
     int idx = blockIdx.x*blockDim.x+threadIdx.x;
     //if ( idx < N) out[idx]=in[idx]*in[idx];
     if ( idx < N) out[idx]= 1;
     out[idx] = randMat[idx];


}

/*
__device__ float  myrand(int xz)
{
     int N = 100;
     unsigned long long int a = 1664525;
     unsigned long long int c = 1013904223;
     unsigned long long int m = 4294967296; // 2^32;
     unsigned long long int x = xz;
     unsigned long long T = (unsigned)time
     for (int n = 0; n < N; n++)
     {
	  x = (a*xz + c)%m;
	  xz = x;
     }
	  
     return float(x)/float(m);
     
}
*/
