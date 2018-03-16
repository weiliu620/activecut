#include "common.h"
__global__ void AnnealingKernel(float* C, 
				const float* R, 
				const unsigned int* gpuFwdMap1,
				const unsigned int* gpuFwdMap2,
				const unsigned short* gpuBwdMap1,
				const unsigned short* gpuBwdMap2,
				paraType par,
				unsigned int N1,
				unsigned int N2,
				float* gpuRandMat,
				dim3 size,
				uint job,
				uint checkerboard);

void GPUAnnealing(float* C, // initial connectivity matrix.
		 float* R,
		 const unsigned int* fwdMap1,
		 const unsigned int* fwdMap2,
		 const unsigned short* bwdMap1,
		 const unsigned short* bwdMap2,
		 paraType par,
		 unsigned int N1,
		 unsigned int N2,
		 const unsigned short* size,
		 float* post2)


{
     if (_DBG >= 2){
	  printf("GPUAnnealing row and col of C: N1 = %d, N2 = %d\n", N1, N2);
	  printf("GPUAnnealing, size = [%d][%d][%d]\n", size[0], size[1], size[2]);
     }
     uint iter; // Gibbs Sampling iteration time.
     uint i, j;
     unsigned int mapSize = size[0]*size[1]*size[2];
     // define simSize as argument of kernel fun.
     dim3 dimSize(size[0], size[1], size[2]);
     // random number matrix.
     float* randMat = (float*) calloc(N1*N2, sizeof(float));
     
     // pointer to device memory.
     float* gpuC;
     float* gpuR;
     unsigned int* gpuFwdMap1;
     unsigned int* gpuFwdMap2;
     unsigned short* gpuBwdMap1;
     unsigned short* gpuBwdMap2;
     float* gpuRandMat;
     
     /* create input and output array on GPU. */
     cudaMalloc((void**) &gpuC, sizeof(float)*N1*N2);
     checkCUDAError("CudaSampling, allocate gpuC.");

     cudaMalloc((void**) &gpuR, sizeof(float)*N1*N2);
     checkCUDAError("CudaSampling, allocate gpuR");

     cudaMalloc((void**) &gpuFwdMap1, sizeof(unsigned int)*mapSize);
     checkCUDAError("CudaSampling, allocate fwdMap1");     
     cudaMalloc((void**) &gpuFwdMap2, sizeof(unsigned int)*mapSize);
     checkCUDAError("CudaSampling, allocate fwdMap2");     
     cudaMalloc((void**) &gpuBwdMap1, sizeof(unsigned short)*N1*3);
     checkCUDAError("CudaSampling, allocate bwdMap1");     
     cudaMalloc((void**) &gpuBwdMap2, sizeof(unsigned short)*N2*3);
     checkCUDAError("CudaSampling, allocate bwdMap2");     

     cudaMalloc((void**) &gpuRandMat, sizeof(float)*N1*N2);
     checkCUDAError("CudaSampling, allocate gpuRandMat");

     /* host to device memory. */
     cudaMemcpy(gpuR, R, sizeof(float)*N1*N2, cudaMemcpyHostToDevice);
     checkCUDAError("CudaSampling, memcpy R");
     cudaMemcpy(gpuC, C, sizeof(float)*N1*N2, cudaMemcpyHostToDevice);
     checkCUDAError("CudaSampling, memcpy gpuC");

     cudaMemcpy(gpuFwdMap1, fwdMap1, sizeof(unsigned int)*mapSize, 
		cudaMemcpyHostToDevice);
     checkCUDAError("CudaSampling, memcpy gpuFwdMap1");
     cudaMemcpy(gpuFwdMap2, fwdMap2, sizeof(unsigned int)*mapSize, 
		cudaMemcpyHostToDevice);
     checkCUDAError("CudaSampling, memcpy gpuFwdMap2");

     cudaMemcpy(gpuBwdMap1, bwdMap1, sizeof(unsigned short)*N1*3, 
		cudaMemcpyHostToDevice);
     checkCUDAError("CudaSampling, memcpy gpuBwdMap1");
     cudaMemcpy(gpuBwdMap2, bwdMap2, sizeof(unsigned short)*N2*3, 
		cudaMemcpyHostToDevice);
     checkCUDAError("CudaSampling, memcpy gpuBwdMap2");


     /* run the kernel function. */
     int gridDimx = N1/BLOCK_SIZE_X + (N1%BLOCK_SIZE_X == 0?0:1);
     int gridDimy = N2/BLOCK_SIZE_Y + (N2%BLOCK_SIZE_Y == 0?0:1);

     dim3 dimBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y);
     dim3 dimGrid(gridDimx, gridDimy);
     if(_DBG >= 1){
	  printf("GPUAnnealing, block size: %dx%d\n", BLOCK_SIZE_X, BLOCK_SIZE_Y);
	  printf("GPUAnnealing, gridsize: %dx%d\n", gridDimx, gridDimy);

     }
     for (iter = 0; iter <= MAXANNITER; iter++){
	  iter++;
	  par.T = log(2) * 1 /log(iter + 1);
	  if (_DBG >= 1 & (iter)%1 == 0){
	       printf("GPUAnnealing, temperature = %f, begin iteration %d...\n", par.T, iter);
	  }
/*
	  // Generate random number.
//	  srand48(time(0));
	  srand48(iter);
	  for (i = 0; i < N1*N2; i++){
	       randMat[i] = drand48();
	  }

	  // send random number to gpu.
	  cudaMemcpy(gpuRandMat, randMat, sizeof(float)*N1*N2, cudaMemcpyHostToDevice);
	  checkCUDAError("GPUAnnealingling, memcpy gpuRandMat");
	  // call kernel for even voxels.
	  AnnealingKernel<<<dimGrid, dimBlock>>>(gpuC, gpuR, gpuFwdMap1, gpuFwdMap2,
						gpuBwdMap1, gpuBwdMap2, par, N1, N2,
						gpuRandMat, dimSize, 0, 0);
	  cudaThreadSynchronize();
	  // call kernel for odd voxels.
	  AnnealingKernel<<<dimGrid, dimBlock>>>(gpuC, gpuR, gpuFwdMap1, gpuFwdMap2,
						gpuBwdMap1, gpuBwdMap2, par, N1, N2,
						gpuRandMat, dimSize, 0, 1);

	  // Send results back to cpu memeory.
	  cudaMemcpy(C, gpuC, sizeof(float)*N1*N2, cudaMemcpyDeviceToHost);
*/
	  // mean field.
	  if (_DBG >= 1 & (iter+1)%1 == 0){
	       printf("GPUAnnealing, mean filed iteration: %d...\n", iter+1);
	  }
	  AnnealingKernel<<<dimGrid, dimBlock>>>(gpuC, gpuR, gpuFwdMap1, gpuFwdMap2,
						gpuBwdMap1, gpuBwdMap2, par, N1, N2,
						gpuRandMat, dimSize, 1, 0);

	  cudaThreadSynchronize();
	  AnnealingKernel<<<dimGrid, dimBlock>>>(gpuC, gpuR, gpuFwdMap1, gpuFwdMap2,
						gpuBwdMap1, gpuBwdMap2, par, N1, N2,
						gpuRandMat, dimSize, 1, 1);
	  /* Send results back to cpu memeory */
	  cudaMemcpy(C, gpuC, sizeof(float)*N1*N2, cudaMemcpyDeviceToHost);

     }

	  /* clean up. */
	  cudaFree(gpuC);
	  cudaFree(gpuR);
	  cudaFree(gpuFwdMap1);
	  cudaFree(gpuFwdMap2);
	  cudaFree(gpuBwdMap1);
	  cudaFree(gpuBwdMap2);
	  cudaFree(gpuRandMat);
	  free(randMat);
}

/* Kernel function  */
__global__ void AnnealingKernel(float* C, 
			       const float* R, 
			       const unsigned int* gpuFwdMap1,
			       const unsigned int* gpuFwdMap2,
			       const unsigned short* gpuBwdMap1,
			       const unsigned short* gpuBwdMap2,
			       paraType par,
			       unsigned int N1,
			       unsigned int N2,
			       float* gpuRandMat,
			       dim3 size,
			       uint job,
			       uint checkerboard)
{     
     unsigned short i;
     unsigned short j;
     unsigned short k;
     uint n;
     float postExp = 0;
     uint n1 = blockIdx.x*blockDim.x + threadIdx.x;
     uint n2 = blockIdx.y*blockDim.y + threadIdx.y;

     float gc1 = 1/(sqrt(2*PI*par.sigma21));
     float gc2 = 1/(sqrt(2*PI*par.sigma22));
     float sum = 0, denom = 0, pu1 = 0, pu2 = 0, postsum = 0;
     float lu1 = 0, lu2 = 0, post1 = 0, post2 = 0;

     // define checkerboard and update image in two steps.
     unsigned short checksum = gpuBwdMap1[n1*3 + 0]
	  + gpuBwdMap1[n1*3 + 1]
	  + gpuBwdMap1[n1*3 + 2]
	  + gpuBwdMap2[n2*3 + 0]
	  + gpuBwdMap2[n2*3 + 1]
	  + gpuBwdMap2[n2*3 + 2];

     if (checksum%2 != checkerboard) {return;}
     // thread fall outside of matrix C, or mask is zero.
     if (n1 >= N1 | n2 >= N2) {return;}

     // image one's neighbors.
     i = gpuBwdMap1[n1*3 + 0];
     j = gpuBwdMap1[n1*3 + 1];
     k = gpuBwdMap1[n1*3 + 2];

     if ((i > 0) && (gpuFwdMap1[(i-1) * size.y*size.z +  j * size.z + k] > 0)){
	  n = gpuFwdMap1[(i-1) * size.y * size.z +  j * size.z + k];
	  sum = sum + C[n*N2 + n2];}

     if ((j > 0) && (gpuFwdMap1[i * size.y*size.z +  (j-1) * size.z + k] > 0)){
	  n = gpuFwdMap1[i * size.y*size.z +  (j-1) * size.z + k];
	  sum = sum + C[n*N2 + n2];}

     if ((k > 0) && (gpuFwdMap1[i * size.y*size.z +  j * size.z + (k-1)] > 0)){
	  n = gpuFwdMap1[i * size.y*size.z +  j * size.z + (k-1)];
	  sum = sum + C[n*N2 + n2];}


     if ((i < size.x-1) && (gpuFwdMap1[(i+1) * size.y*size.z +  j * size.z + k] > 0)){
	  n = gpuFwdMap1[(i+1) * size.y*size.z +  j * size.z + k];
	  sum = sum + C[n*N2 + n2];}

     if ((j < size.y-1) && (gpuFwdMap1[i * size.y*size.z +  (j+1) * size.z + k] > 0)){
	  n = gpuFwdMap1[i * size.y*size.z +  (j+1) * size.z + k];
	  sum = sum + C[n*N2 + n2];}

     if ((k < size.z-1) && (gpuFwdMap1[i * size.y*size.z +  j * size.z + (k+1)] > 0)){
	  n = gpuFwdMap1[i * size.y*size.z +  j * size.z + (k+1)];
	  sum = sum + C[n*N2 + n2];}

     // image 2's neighbors.
     i = gpuBwdMap2[n2*3 + 0];
     j = gpuBwdMap2[n2*3 + 1];
     k = gpuBwdMap2[n2*3 + 2];

     if ((i > 0) && (gpuFwdMap2[(i-1) * size.y*size.z +  j * size.z + k] > 0)){
	  n = gpuFwdMap2[(i-1) * size.y * size.z +  j * size.z + k];
	  sum = sum + C[n1*N2 + n];}

     if ((j > 0) && (gpuFwdMap2[i * size.y*size.z +  (j-1) * size.z + k] > 0)){
	  n = gpuFwdMap2[i * size.y*size.z +  (j-1) * size.z + k];
	  sum = sum + C[n1*N2 + n];}


     if ((k > 0) && (gpuFwdMap2[i * size.y*size.z +  j * size.z + (k-1)] > 0)){
	  n = gpuFwdMap2[i * size.y*size.z +  j * size.z + (k-1)];
	  sum = sum + C[n1*N2 + n];}

     if ((i < size.x-1) && (gpuFwdMap2[(i+1) * size.y*size.z +  j * size.z + k] > 0)){
	  n = gpuFwdMap2[(i+1) * size.y*size.z +  j * size.z + k];
	  sum = sum + C[n1*N2 + n];}

     if ((j < size.y-1) && (gpuFwdMap2[i * size.y*size.z +  (j+1) * size.z + k] > 0)){
	  n = gpuFwdMap2[i * size.y*size.z +  (j+1) * size.z + k];
	  sum = sum + C[n1*N2 + n];}

     if ((k < size.z-1) && (gpuFwdMap2[i * size.y*size.z +  j * size.z + (k+1)] > 0)){
	  n = gpuFwdMap2[i * size.y*size.z +  j * size.z + (k+1)];
	  sum = sum + C[n1*N2 + n];}

     pu1 = par.beta * (-1)* sum;
     pu2 = par.beta * sum;
	  
     lu1 = pow(R[n1*N2+n2]-par.mu1, 2)/(2*par.sigma21)+ log(sqrt(par.sigma21));
     lu2 = pow(R[n1*N2+n2]-par.mu2, 2)/(2*par.sigma22) + log(sqrt(par.sigma22));
     // compute posterior and include annealing parameter T.
     // May need to consider precision issue!
     post1 = exp(-(pu1 + lu1)/par.T);
     post2 = exp(-(pu2 + lu2)/par.T);
     postsum = post1 + post2;
     post1 = post1/postsum;
     post2 = post2/postsum;
     
     if (job == 1){ // expectation. mean field.
	  C[n1*N2+n2] = post2 - post1;
     }
     else if (job == 0){ // sampling
	  if (gpuRandMat[n1*N2+n2] < post1){
	       C[n1*N2+n2] =  -1;
	  }

	  else{
	       C[n1*N2+n2] =  1;
	  }
     }
     else if (job == 2){ // compute posterior. p(c == -1)
	  gpuRandMat[n1*N2+n2] = post2;
     }
     else{

     }
}
