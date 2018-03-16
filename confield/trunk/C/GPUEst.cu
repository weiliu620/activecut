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
     double Jacob = 0;
     double DELL = 0;
     double step = 1e6;
     unsigned int n = 0, n1 = 0, n2 = 0;
     double nsum = 0;
     unsigned int mapSize = size[0]*size[1]*size[2];
     // connectivty matrix size.
     unsigned int BN = N1 * (N1+1)/2;
     // define simSize as argument of kernel fun.
     dim3 dimSize(size[0], size[1], size[2]);
     // create nsum on host. to save results from gpu.
     float* sumMat = (float*) malloc(BN*sizeof(float));

     // pointer to device memory.
     float* gpuC;
     unsigned int* gpuFwdMap1;
     unsigned int* gpuFwdMap2;
     unsigned short* gpuBwdMap1;
     unsigned short* gpuBwdMap2;
     float* gpuSumMat;

     /* create input and output array on GPU. */
     cudaMalloc((void**) &gpuC, sizeof(float)*BN);
     checkCUDAError("GPUEst, allocate gpuC.");

     cudaMalloc((void**) &gpuFwdMap1, sizeof(unsigned int)*mapSize);
     checkCUDAError("GPUEst, allocate fwdMap1");     
     cudaMalloc((void**) &gpuFwdMap2, sizeof(unsigned int)*mapSize);
     checkCUDAError("GPUEst, allocate fwdMap2");     
     cudaMalloc((void**) &gpuBwdMap1, sizeof(unsigned short)*N1*3);
     checkCUDAError("GPUEst, allocate bwdMap1");     
     cudaMalloc((void**) &gpuBwdMap2, sizeof(unsigned short)*N2*3);
     checkCUDAError("GPUEst, allocate bwdMap2");     
     cudaMalloc((void**) &gpuSumMat, sizeof(float)*BN);
     checkCUDAError("GPUEst, allocate gpuSumMat.");

     /* host to device memory. */
     cudaMemcpy(gpuC, C, sizeof(float)*BN, cudaMemcpyHostToDevice);
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
     cudaMemcpy(gpuSumMat, sumMat,  sizeof(float)*BN, 
		cudaMemcpyHostToDevice);
     checkCUDAError("GPUEst, memcpy gpuSumMat");

     /* run the kernel function. */
     int gridDimx = N1/BLOCK_SIZE_X + (N1%BLOCK_SIZE_X == 0?0:1);
     int gridDimy = N2/BLOCK_SIZE_Y + (N2%BLOCK_SIZE_Y == 0?0:1);

     dim3 dimBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y);
     dim3 dimGrid(gridDimx, gridDimy);
     if (_DBG >= 1){
	  printf("GPUEst, block size: %dx%d\n", BLOCK_SIZE_X, BLOCK_SIZE_Y);
	  printf("GPUEst, gridsize: %dx%d\n", gridDimx, gridDimy);
     }

     // call kernel, only to compute sum x_j for each x_i.
     EstKernel<<<dimGrid, dimBlock>>>(gpuC, gpuFwdMap1, gpuFwdMap2,
				      gpuBwdMap1, gpuBwdMap2, *par, N1, N2,
				      dimSize, gpuSumMat);
     cudaMemcpy(sumMat, gpuSumMat, sizeof(float)*BN, cudaMemcpyDeviceToHost);

     // estimate beta by Newton's method.
     while((step > CUEPS) && (par->beta < CUMAXBETA)){
	  Jacob = 0;
	  DELL = 0;

	  // iterate all data points.
	  for (n1 = 0; n1 < N1; n1++) {
	       for (n2 = n1+1; n2 < N2; n2++) {
		    nsum = GTI(sumMat, n1, n2, N1);
		    // 2nd derivative.
		    Jacob = Jacob + 
			 nsum * nsum * (tanh(par->beta*nsum)*tanh(par->beta*nsum) - 1);
		    DELL = DELL + nsum * (GTI(C, n1, n2, N2) - tanh(par->beta*nsum));
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
     if (n1 >= N1 | n2 >= N2 | n1 > n2) {return;}

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
			 n = gpuFwdMap1[ti * size.y*size.z +  tj * size.z + tk];
			 if (n <= n2) {
			      sum = sum + GTI(C, n, n2, N2);
			 }
			 else {
			      // lower traingle of matrix C. use upper triangle value.
			      sum = sum + GTI(C, n2, n, N2);
			 }
		    }
	       }
	  }
     }
     sum = sum - GTI(C, n1, n2, N2);

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
			&& (gpuFwdMap2[ti * size.y * size.z +  tj * size.z + tk] > 0)){
			 n = gpuFwdMap2[ti * size.y * size.z +  tj * size.z + tk];
			 if (n1 <= n) {
			      sum = sum + GTI(C, n1, n, N2);
			 }
			 else {
			      sum = sum + GTI(C, n, n1, N2);
			 }
		    }
	       }
	  }
     }
     sum = sum - GTI(C, n1, n2, N2);

     GTI(gpuSumMat, n1, n2, N2) = sum;

}
