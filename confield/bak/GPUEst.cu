#include "common.h"

__global__ void EstKernel(const float* C, 
			  const unsigned int* gpuFwdMap1,
			  const unsigned int* gpuFwdMap2,
			  const unsigned short* gpuBwdMap1,
			  const unsigned short* gpuBwdMap2,
			  paraType par,
			  unsigned int N1,
			  unsigned int N2,
			  dim3 size,
			  float* gpuSumMat);

void GPUEst(const float* C,
	    const unsigned int* fwdMap1,
	    const unsigned int* fwdMap2,
	    const unsigned short* bwdMap1,
	    const unsigned short* bwdMap2,
	    paraType* par,
	    unsigned int N1,
	    unsigned int N2,
	    const unsigned short* size)
{
     double Jacob = 0, Jacob_alpha = 0;
     double DELL = 0, DELL_alpha = 0;
     double step = 1e6;
     double step_alpha = 1e6;
     unsigned int n = 0, n1 = 0, n2 = 0;
     double nsum = 0;
     unsigned int mapSize = size[0]*size[1]*size[2];
     // define simSize as argument of kernel fun.
     dim3 dimSize(size[0], size[1], size[2]);
     // create nsum on host. to save results from gpu.
     float* sumMat = (float*) malloc(N1*N2*sizeof(float));

     // pointer to device memory.
     float* gpuC;
     unsigned int* gpuFwdMap1;
     unsigned int* gpuFwdMap2;
     unsigned short* gpuBwdMap1;
     unsigned short* gpuBwdMap2;
     float* gpuSumMat;

     /* create input and output array on GPU. */
     cudaMalloc((void**) &gpuC, sizeof(float)*N1*N2);
     checkCUDAError("GPUEst, allocate gpuC.");

     cudaMalloc((void**) &gpuFwdMap1, sizeof(unsigned int)*mapSize);
     checkCUDAError("GPUEst, allocate fwdMap1");     
     cudaMalloc((void**) &gpuFwdMap2, sizeof(unsigned int)*mapSize);
     checkCUDAError("GPUEst, allocate fwdMap2");     
     cudaMalloc((void**) &gpuBwdMap1, sizeof(unsigned short)*N1*3);
     checkCUDAError("GPUEst, allocate bwdMap1");     
     cudaMalloc((void**) &gpuBwdMap2, sizeof(unsigned short)*N2*3);
     checkCUDAError("GPUEst, allocate bwdMap2");     
     cudaMalloc((void**) &gpuSumMat, sizeof(float)*N1*N2);
     checkCUDAError("GPUEst, allocate gpuSumMat.");

     /* host to device memory. */
     cudaMemcpy(gpuC, C, sizeof(float)*N1*N2, cudaMemcpyHostToDevice);
     checkCUDAError("GPUEst, memcpy gpuC");

     cudaMemcpy(gpuFwdMap1, fwdMap1, sizeof(unsigned int)*mapSize, 
		cudaMemcpyHostToDevice);
     checkCUDAError("GPUEst, memcpy gpuFwdMap1");
     cudaMemcpy(gpuFwdMap2, fwdMap2, sizeof(unsigned int)*mapSize, 
		cudaMemcpyHostToDevice);
     checkCUDAError("GPUEst, memcpy gpuFwdMap2");

     cudaMemcpy(gpuBwdMap1, bwdMap1, sizeof(unsigned short)*N1*3, 
		cudaMemcpyHostToDevice);
     checkCUDAError("GPUEst, memcpy gpuBwdMap1");
     cudaMemcpy(gpuBwdMap2, bwdMap2, sizeof(unsigned short)*N2*3, 
		cudaMemcpyHostToDevice);
     checkCUDAError("GPUEst, memcpy gpuBwdMap2");
     cudaMemcpy(gpuSumMat, sumMat,  sizeof(float)*N1*N2, 
		cudaMemcpyHostToDevice);
     checkCUDAError("GPUEst, memcpy gpuSumMat");

     /* run the kernel function. */
     int gridDimx = N1/BLOCK_SIZE_X + (N1%BLOCK_SIZE_X == 0?0:1);
     int gridDimy = N2/BLOCK_SIZE_Y + (N2%BLOCK_SIZE_Y == 0?0:1);

     dim3 dimBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y);
     dim3 dimGrid(gridDimx, gridDimy);
     if (_DBG >= 1){
	  printf("GPUEst, block size: %dx%d. gridsize: %dx%d\n", BLOCK_SIZE_X, BLOCK_SIZE_Y, gridDimx, gridDimy);
     }
     // call kernel, only to compute sum x_j for each x_i.
     EstKernel<<<dimGrid, dimBlock>>>(gpuC, gpuFwdMap1, gpuFwdMap2,
				      gpuBwdMap1, gpuBwdMap2, *par, N1, N2,
				      dimSize, gpuSumMat);
     cudaMemcpy(sumMat, gpuSumMat, sizeof(float)*N1*N2, cudaMemcpyDeviceToHost);

     // estimate alpha and beta by Newton's method.
     // estimate alpha
     float alpha_temp = par->alpha;
     while((abs(step_alpha) > CUEPS) && (abs(alpha_temp * 2) < MAXEXP)){
	  Jacob_alpha = 0;
	  DELL_alpha = 0;
	  for (n1 = 0; n1 < N1; n1++){
	       for (n2 = 0; n2 < N2; n2++) {
		    if (n1 != n2) {
			 nsum = sumMat[n1*N2+n2];
			 // 2nd derivative.
			 Jacob_alpha = Jacob_alpha + 
			      tanh(alpha_temp + par->beta*nsum)*tanh(alpha_temp + par->beta*nsum) - 1;
			 DELL_alpha = DELL_alpha + C[n1*N2+n2] - tanh(alpha_temp + par->beta*nsum);
		    }
	       }
	  }
	  if (Jacob_alpha != 0) {
	       step_alpha = -DELL_alpha/Jacob_alpha;
	       alpha_temp = alpha_temp + step_alpha;
	       printf("GPUEst.cu: alpha = %f\n", par->alpha);
	  }
	  else {
	       fprintf(stderr, "GPUEst.cu: DELL_alpha zero. alpha does not change.");
	  }
     }
     par->alpha = alpha_temp;

     // estimate beta.
     while((abs(step) > CUEPS) && (par->beta < CUMAXBETA)){
	  Jacob = 0;
	  DELL = 0;
	  for (n1 = 0; n1 < N1; n1++){
	       for (n2 = 0; n2 < N2; n2++) {
		    if (n1 != n2) {
			 nsum = sumMat[n1*N2+n2];
			 // 2nd derivative.
			 Jacob = Jacob + 
			      nsum * nsum * (tanh(par->beta*nsum)*tanh(par->alpha + par->beta*nsum) - 1);
			 DELL = DELL + nsum * (C[n1*N2+n2] - tanh(par->alpha + par->beta*nsum));
		    }
	       }
	  }
	  if (Jacob != 0){
	       step = -DELL/Jacob;
	       par->beta = par->beta + step;
	       if (par->beta >= CUMAXBETA){
		    printf("GPUEst.cu: beta too large. Set Beta to largest value.\n");
		    par->beta = CUMAXBETA;
	       }
	  }
	  else{
	       fprintf(stderr, "GPUEst.cu: DELL zero. Netwon's method deviced by zero.\n");
	       printf("GPUEst.cu: jacobian is zero. Probably newBeta is too large, \n so we just set beta to largest value.\n");
	       par->beta = CUMAXBETA;
	  }
	  if (_DBG >= 1){
	       printf("GPUEst.cu: beta = %f\n", par->beta);
	  }

     }
     if (_DBG >= 1){
	  printf("GPUEst.cu: beta = %f\n", par->beta);
     }



     /* clean up. */
     cudaFree(gpuC);
     cudaFree(gpuFwdMap1);
     cudaFree(gpuFwdMap2);
     cudaFree(gpuBwdMap1);
     cudaFree(gpuBwdMap2);
     cudaFree(gpuSumMat);
     free(sumMat);

}

/* Kernel function  */
__global__ void EstKernel(const float* C, 
			  const unsigned int* gpuFwdMap1,
			  const unsigned int* gpuFwdMap2,
			  const unsigned short* gpuBwdMap1,
			  const unsigned short* gpuBwdMap2,
			  paraType par,
			  unsigned int N1,
			  unsigned int N2,
			  dim3 size,
			  float* gpuSumMat)

{     
     unsigned short i;
     unsigned short j;
     unsigned short k;
     int di, dj, dk;
     int ti, tj, tk;
     uint n;
     uint n1 = blockIdx.x*blockDim.x + threadIdx.x;
     uint n2 = blockIdx.y*blockDim.y + threadIdx.y;

     float gc1 = 1/(sqrt(2*PI*par.sigma21));
     float gc2 = 1/(sqrt(2*PI*par.sigma22));
     float sum = 0;

     // thread fall outside of matrix C, or mask is zero.
     if (n1 >= N1 | n2 >= N2) {return;}

     // image one's neighbors.
     i = gpuBwdMap1[n1*3 + 0];
     j = gpuBwdMap1[n1*3 + 1];
     k = gpuBwdMap1[n1*3 + 2];

     for (di = -1; di <= 1; di ++){
	  for (dj = -1; dj <= 1; dj ++){
	       for (dk = -1; dk <= 1; dk ++){
		    ti = i + di;
		    tj = j + dj;
		    tk = k + dk;
		    if ((ti >= 0 && ti < size.x
			 && tj >= 0 && tj < size.y
			 && tk >= 0 && tk < size.z)
			 && (gpuFwdMap1[ti * size.y*size.z +  tj * size.z + tk] > 0)){
		    n = gpuFwdMap1[ti * size.y * size.z +  tj * size.z + tk];
		    sum = sum + C[n*N2 + n2];
		    }
	       }
	  }
     }
     sum = sum - C[n1*N2 + n2];
     // image 2's neighbors.
     i = gpuBwdMap2[n2*3 + 0];
     j = gpuBwdMap2[n2*3 + 1];
     k = gpuBwdMap2[n2*3 + 2];
     for (di = -1; di <= 1; di ++){
	  for (dj = -1; dj <= 1; dj ++){
	       for (dk = -1; dk <= 1; dk ++){
		    ti = i + di;
		    tj = j + dj;
		    tk = k + dk;
		    if ((ti >= 0 && ti < size.x
			 && tj >= 0 && tj < size.y
			 && tk >= 0 && tk < size.z)
			 && (gpuFwdMap2[ti * size.y*size.z +  tj * size.z + tk] > 0)){
		    n = gpuFwdMap2[ti * size.y * size.z +  tj * size.z + tk];
		    sum = sum + C[n1*N2 + n];
		    }
	       }
	  }
     }
     sum = sum - C[n1*N2 + n2];
     
     gpuSumMat[n1*N2+n2] = sum;
#if __DEVICE_EMULATION__
     if (n1 == 0 && n2 == 0){

     }
#endif


}
