#include "common.h"

__global__ void SamplingKernel(float* C, 
			       const float* R, 
			       const unsigned int* fwdMap1,
			       const unsigned int* fwdMap2,
			       const unsigned short* bwdMap1,
			       const unsigned short* bwdMap2,
			       paraType par,
			       unsigned int N1,
			       unsigned int N2,
			       float* gpuRandMat,
			       dim3  size,
			       uint job,
			       uint checkerboard);

void GPUSampling(float* C, // initial connectivity matrix.
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
     if (_DBG >= 3){
	  printf("GPUSampling row and col of C: N1 = %d, N2 = %d\n", N1, N2);
	  printf("GPUSampling, size = [%d][%d][%d]\n", size[0], size[1], size[2]);
     }
     uint iter; // Gibbs Sampling iteration time.
     uint i, j;

     unsigned int mapSize = size[0]*size[1]*size[2];
     // define simSize as argument of kernel fun.
     dim3 dimSize(size[0], size[1], size[2]);
     // random number matrix.
     float* randMat = (float*) calloc(N1*N2, sizeof(float));
     // Read in parameter file.
     
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
     if(_DBG >= 3){
	  printf("GPUSampling, block size: %dx%d\n", BLOCK_SIZE_X, BLOCK_SIZE_Y);
	  printf("GPUsampling, gridsize: %dx%d\n", gridDimx, gridDimy);
     }

     for (iter = 0; iter < GIBBSMAXITER; iter++){
	  if (_DBG >= 1 & (iter+1)%10 == 0){
	       printf("GPUSampling, begin iteration %d...\n", iter+1);
	  }

	  // Generate random number.
//	  srand48(time(0));
	  srand48(iter);
	  for (i = 0; i < N1*N2; i++){
	       randMat[i] = drand48();
	  }

	  // send random number to gpu.
	  cudaMemcpy(gpuRandMat, randMat, sizeof(float)*N1*N2, cudaMemcpyHostToDevice);
	  checkCUDAError("CudaSampling, memcpy gpuRandMat");
	  // call kernel for even voxels.
	  SamplingKernel<<<dimGrid, dimBlock>>>(gpuC, gpuR, gpuFwdMap1, gpuFwdMap2,
						gpuBwdMap1, gpuBwdMap2, par, N1, N2,
						gpuRandMat, dimSize, 0, 0);
	  cudaThreadSynchronize();
	  // call kernel for odd voxels.
	  SamplingKernel<<<dimGrid, dimBlock>>>(gpuC, gpuR, gpuFwdMap1, gpuFwdMap2,
						gpuBwdMap1, gpuBwdMap2, par, N1, N2,
						gpuRandMat, dimSize, 0, 1);

	  // Send results back to cpu memeory.
	  cudaMemcpy(C, gpuC, sizeof(float)*N1*N2, cudaMemcpyDeviceToHost);
	  if (_DBG >= 4){
	       printf("GPUSampling, after samplingKernel iteration %d. C: \n", iter);
	       for (i = RMIN; i < RMAX; i++){
		    for (j = CMIN; j < CMAX; j++){
			 printf("[%d][%d]=%1.2f ", i, j, C[i*N2 + j]);
		    }
		    printf("\n");
	       }
	  }
     }

     if (_DBG >= 4){
	  printf("GPUSampling, all sampling iteration done.C: \n");
	  for (i = RMIN; i < RMAX; i++){
	       for (j = CMIN; j < CMAX; j++){
		    printf("[%d][%d]=%1.2f ", i, j, C[i*N2 + j]);
	       }
	       printf("\n");
	  }
	  
	  FieldViewer(C, N1, N2, fwdMap1, fwdMap2, "GPUSampling: Sampled C");

     }

     // Now compute mean field.
     for (iter = 0; iter <= MEANFIELDITER; iter++){
	  iter++;
	  if (_DBG >= 1 & iter%5 == 0){
	       printf("GPUSampling, mean filed iteration: %d...\n", iter);
	  }
	  SamplingKernel<<<dimGrid, dimBlock>>>(gpuC, gpuR, gpuFwdMap1, gpuFwdMap2,
						gpuBwdMap1, gpuBwdMap2, par, N1, N2,
						gpuRandMat, dimSize, 1, 0);

	  cudaThreadSynchronize();
	  SamplingKernel<<<dimGrid, dimBlock>>>(gpuC, gpuR, gpuFwdMap1, gpuFwdMap2,
						gpuBwdMap1, gpuBwdMap2, par, N1, N2,
						gpuRandMat, dimSize, 1, 1);

	  /* Send results back to cpu memeory */
	  cudaMemcpy(C, gpuC, sizeof(float)*N1*N2, cudaMemcpyDeviceToHost);

     }

     if (_DBG >= 3){
	  printf("Mean field all iteration done. C:\n");
	  for (i = RMIN; i < RMAX; i++){
	       for (j = CMIN; j < CMAX; j++){
		    printf("[%d][%d]=%1.2f ", i, j, C[i*N2 + j]);
	       }
	       printf("\n");
	  }

	  FieldViewer(C, N1, N2, fwdMap1, fwdMap2, "GPUSampling: Mean Field");

     }

     
     // Compute posterior probability p(c = -1), i.e. no connectivity, 
     // and save it in randMat. ( We dont bother use another argument, 
     // because the kernel function does not use randMat when computing 
     // posterior, so randMat is a good choice for other usage: saving 
     // the posterior array.
     SamplingKernel<<<dimGrid, dimBlock>>>(gpuC, gpuR, gpuFwdMap1, gpuFwdMap2,
					   gpuBwdMap1, gpuBwdMap2, par, N1, N2,
					   gpuRandMat, dimSize, 2, 0);     
     cudaThreadSynchronize();
     SamplingKernel<<<dimGrid, dimBlock>>>(gpuC, gpuR, gpuFwdMap1, gpuFwdMap2,
					   gpuBwdMap1, gpuBwdMap2, par, N1, N2,
					   gpuRandMat, dimSize, 2, 1);     
     // send randMat back. It contains posterior info.
     cudaMemcpy(post2, gpuRandMat, sizeof(float)*N1*N2, cudaMemcpyDeviceToHost);

     if (_DBG >= 4){
	  printf("After computing posterior. posterior p(c = 1):\n");
	  for (i = RMIN; i < RMAX; i++){
	       for (j = CMIN; j < CMAX; j++){
		    printf("[%d][%d]=%1.2f ", i, j, post2[i*N2 + j]);
	       }
	       printf("\n");
	  }

	  FieldViewer(post2, N1, N2, fwdMap1, fwdMap2, "GPUSampling: Sampled C");
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
__global__ void SamplingKernel(float* C, 
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
     int di, dj, dk;
     int ti, tj, tk;
     uint n;
     float postExp = 0;
     uint n1 = blockIdx.x*blockDim.x + threadIdx.x;
     uint n2 = blockIdx.y*blockDim.y + threadIdx.y;

     float gc1 = 1/(sqrt(2*PI*par.sigma21));
     float gc2 = 1/(sqrt(2*PI*par.sigma22));
     float sum = 0, denom = 0, pr1 = 0, pr2 = 0, postsum = 0;
     float lh1 = 0, lh2 = 0, post1 = 0, post2 = 0;

     // thread fall outside of matrix C, or mask is zero.
     if (n1 >= N1 | n2 >= N2) {return;}

     // define checkerboard and update image in two steps.
     unsigned short checksum = gpuBwdMap1[n1*3 + 0]
	  + gpuBwdMap1[n1*3 + 1]
	  + gpuBwdMap1[n1*3 + 2]
	  + gpuBwdMap2[n2*3 + 0]
	  + gpuBwdMap2[n2*3 + 1]
	  + gpuBwdMap2[n2*3 + 2];

     if (checksum%2 != checkerboard) {return;}


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


     denom = 2*cosh(par.beta*sum);

     // if exponential is too big, just mannually set pr1 and pr2.
     if (par.beta*sum > MAXEXP){
	  pr1 = 0;
	  pr2 = 1;
     }
     else if(par.beta*sum < -MAXEXP){
	  pr1 = 1;
	  pr2 = 0;
     }
     else{
	  pr1 = exp(par.beta * (-1)* sum)/denom;
	  pr2 = exp(par.beta * sum)/denom;
     }
     // if posterior is too big. just mannualy set post1 and post2.
     postExp = log(par.sigma21/par.sigma22)
	  + 2*par.beta*sum
	  - pow(R[n1*N2+n2]-par.mu1, 2)
	  + pow(R[n1*N2+n2]-par.mu2, 2);
     if (postExp > MAXEXP){
	  post1 = 0;
	  post2 = 1;
     }
     else if(postExp < - MAXEXP){
	  post1 = 1;
	  post2 = 0;
     }
     else{
	  
	  lh1 = gc1 * exp(-pow(R[n1*N2+n2]-par.mu1, 2)/(2*par.sigma21));
	  lh2 = gc2 * exp(-pow(R[n1*N2+n2]-par.mu2, 2)/(2*par.sigma22));
          
	  post1 = pr1 * lh1;
	  post2 = pr2 * lh2;
	  postsum = post1 + post2;
	  post1 = post1/postsum;
	  post2 = post2/postsum;
     }

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
     else{ };
 
}
