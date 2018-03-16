#include "common.h"
__global__ void DetrendKernel(float* im, unsigned int N, unsigned int D);
int GPUDetrend(float* im, unsigned int N, unsigned int D)
{
     if (_DBG >= 3){
	  printf("GPUdetrend, N = %d, D = %d\n", N, D);
     }

     float* gpuIm;
     /* create input and output array on GPU. */
     cudaMalloc((void**) &gpuIm, sizeof(float)*N*D);
     checkCUDAError("GPUdetrend, allocate gpuIm.");
     
     /* host to device memory. */
     cudaMemcpy(gpuIm, im, sizeof(float)*N*D, cudaMemcpyHostToDevice);
     checkCUDAError("GPUdetrend, memecpy gpuIm");
     
     int gridDimx = N/BLOCK_SIZE_X + (N%BLOCK_SIZE_X == 0?0:1);

     dim3 dimBlock(BLOCK_SIZE_X, 1);
     dim3 dimGrid(gridDimx, 1);
     if(_DBG >= 2){
	  printf("GPUDetrend, block size: %dx%d\n", BLOCK_SIZE_X, 1);
	  printf("GPUDetrend, gridsize: %dx%d\n", gridDimx, 1);

     }
     // call kernel for detrending.
     DetrendKernel<<<dimGrid, dimBlock>>>(gpuIm, N, D);
     // send data back to cpu memory.
     cudaMemcpy(im, gpuIm, sizeof(float)*N*D, cudaMemcpyDeviceToHost);
     cudaFree(gpuIm);
     
}

__global__ void DetrendKernel(float* im, unsigned int N, unsigned int D)
{
     int d = 0;
     uint n = blockIdx.x*blockDim.x + threadIdx.x;
     float x_bar = ((float)D-1)/2;
     float y_bar = 0;
     float SSx =0;
     float SSxy = 0;
     float m = 0;
     float b = 0; // coefficients for the linear trend.
     if (n >= N) { return;}
     // least square estimation.
     for (d = 0; d < D; d++){
	  y_bar = y_bar + im[n*D+d]/D;
     }

     for (d = 0; d < D; d++){
	  SSx = SSx + (d - x_bar)*(d - x_bar);
	  SSxy = SSxy + (d - x_bar)*(im[n*D+d] - y_bar);
     }
     
     m = SSxy/SSx;
     b = y_bar - m * x_bar;
     
     // detrending.
     for (d = 0; d < D; d++){
	  im[n*D+d] = im[n*D+d] - (m*d + b);
     }
          
}
