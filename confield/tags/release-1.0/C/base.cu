#include "cuda.h"
#include <math.h>
#include <stdio.h>
__global__ void square_elements(float* in, float* out, int M, int N);
int main()
{
     uint M = 5000;
     uint N = 5000;
     float* inputPtr = (float*)malloc(sizeof(float)*M*N);
     float* outputPtr = (float*)malloc(sizeof(float)*M*N);
     
     int m, n;
     for (m = 0; m < M; m++){
	  for (n = 0; n < N; n++){
	       *(inputPtr + m*N + n) = m+n;
	  }
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
     
     /* run the kernel function. */
     int blockSize = 256;
     int nBlocks = (M*N)/blockSize + ((M*N)%blockSize == 0?0:1);
     printf("blockSize: %d, nBlocks = %d\n", blockSize, nBlocks);
     dim3 dimBlock(blockSize);
     dim3 dimGrid(nBlocks);
     square_elements<<<dimGrid, dimBlock>>>(inputGPUPtr, outputGPUPtr, M, N);
     
     /* Send results back to cpu memeory */
     cudaMemcpy(outputPtr, outputGPUPtr, sizeof(float)*M*N, cudaMemcpyDeviceToHost);

     /* clean up. */
     cudaFree(inputGPUPtr);
     cudaFree(outputGPUPtr);
     free(inputPtr);
     free(outputPtr);

     /* Scratch pad.  */
}


/* Kernel to square elements of the array on the GPU */
__global__ void square_elements(float* in, float* out, int M, int N)
{
     int idx = blockIdx.x*blockDim.x+threadIdx.x;
     if ( idx < N) out[idx]= in[idx] * in[idx];

}
