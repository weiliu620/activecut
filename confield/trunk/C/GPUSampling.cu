#include "common.h"

__global__ void SamplingKernel(float* C, 
			       const float* im,
			       const unsigned int* fwdMap1,
			       const unsigned int* fwdMap2,
			       const unsigned short* bwdMap1,
			       const unsigned short* bwdMap2,
			       paraType par,
			       unsigned int N1,
			       unsigned int N2,
			       unsigned short TL, // time series length.
			       float* gpuRandMat,
			       dim3  size,
			       uint job,
			       uint checkerboard);
__device__ float GetCorr(const float* im, 
			 unsigned int N, 
			 unsigned int M,  // time series length
			 int idx1, 
			 int idx2);

void GPUSampling(float* C, // initial connectivity matrix.
		 const float* im,
		 const unsigned int* fwdMap1,
		 const unsigned int* fwdMap2,
		 const unsigned short* bwdMap1,
		 const unsigned short* bwdMap2,
		 paraType par,
		 unsigned int N1,
		 unsigned int N2,
		 unsigned short TL, // time series length.
		 const unsigned short* size)


{
     if (_DBG >= 3){
	  printf("GPUSampling row and col of C: N1 = %d, N2 = %d\n", N1, N2);
	  printf("GPUSampling, size = [%d][%d][%d]\n", size[0], size[1], size[2]);
     }
     uint iter; // Gibbs Sampling iteration time.
     uint i, j;

     // connectivty matrix size.
     unsigned int BN = N1 * (N1+1)/2;
     unsigned int mapSize = size[0]*size[1]*size[2];
     // define simSize as argument of kernel fun.
     dim3 dimSize(size[0], size[1], size[2]);
     // random number matrix.
     float* randMat = (float*) calloc(BN, sizeof(float));
     
     // pointer to device memory.
     float* gpuC;
     float* gpuIm;
     unsigned int* gpuFwdMap1;
     unsigned int* gpuFwdMap2;
     unsigned short* gpuBwdMap1;
     unsigned short* gpuBwdMap2;
     float* gpuRandMat;
     

     /* create input and output array on GPU. */
     cudaMalloc((void**) &gpuC, sizeof(float)*BN);
     checkCUDAError("CudaSampling, allocate gpuC.");

     cudaMalloc((void**) &gpuIm, sizeof(float)*N1*TL);
     checkCUDAError("CudaSampling, allocate gpuR");

     cudaMalloc((void**) &gpuFwdMap1, sizeof(unsigned int)*mapSize);
     checkCUDAError("CudaSampling, allocate fwdMap1");     
     cudaMalloc((void**) &gpuFwdMap2, sizeof(unsigned int)*mapSize);
     checkCUDAError("CudaSampling, allocate fwdMap2");     
     cudaMalloc((void**) &gpuBwdMap1, sizeof(unsigned short)*N1*3);
     checkCUDAError("CudaSampling, allocate bwdMap1");     
     cudaMalloc((void**) &gpuBwdMap2, sizeof(unsigned short)*N2*3);
     checkCUDAError("CudaSampling, allocate bwdMap2");     



     /* host to device memory. */
     cudaMemcpy(gpuC, C, sizeof(float)*BN, cudaMemcpyHostToDevice);
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

     cudaMemcpy(gpuIm, im, sizeof(float)*N1*TL, cudaMemcpyHostToDevice);
     checkCUDAError("CudaSampling,  copy im to gpuIm");     

     /* run the kernel function. */
     int gridDimx = N1/BLOCK_SIZE_X + (N1%BLOCK_SIZE_X == 0?0:1);
     int gridDimy = N2/BLOCK_SIZE_Y + (N2%BLOCK_SIZE_Y == 0?0:1);

     dim3 dimBlock(BLOCK_SIZE_X, BLOCK_SIZE_Y);
     dim3 dimGrid(gridDimx, gridDimy);
     if(_DBG >= 3){
	  printf("GPUSampling, block size: %dx%d\n", BLOCK_SIZE_X, BLOCK_SIZE_Y);
	  printf("GPUsampling, gridsize: %dx%d\n", gridDimx, gridDimy);
     }

     if (DO_SAMPLING) {
	  // allocate gpu memory for random number.
	  cudaMalloc((void**) &gpuRandMat, sizeof(float)*BN);
	  checkCUDAError("CudaSampling, allocate gpuRandMat");

	  for (iter = 1; iter < GIBBSMAXITER; iter++){
	       if (_DBG >= 1 & iter%1 == 0){
		    printf("GPUSampling, begin iteration %d...\n", iter);
	       }
	       
	       // Generate random number.
	       //srand48(time(0));
	       srand48(iter);
	       for (i = 0; i < BN; i++){
		    randMat[i] = drand48();
	       }
	       
	       // send random number to gpu.
	       cudaMemcpy(gpuRandMat, randMat, sizeof(float)*BN, cudaMemcpyHostToDevice);
	       checkCUDAError("CudaSampling, memcpy gpuRandMat");
	       // call kernel for even voxels.
	       SamplingKernel<<<dimGrid, dimBlock>>>(gpuC, gpuIm, gpuFwdMap1, gpuFwdMap2,
						     gpuBwdMap1, gpuBwdMap2, par, N1, N2,
						     TL,
						     gpuRandMat, dimSize, 0, 0);
	       cudaThreadSynchronize();
	       // call kernel for odd voxels.
	       SamplingKernel<<<dimGrid, dimBlock>>>(gpuC, gpuIm, gpuFwdMap1, gpuFwdMap2,
						     gpuBwdMap1, gpuBwdMap2, par, N1, N2,
						     TL, 
						     gpuRandMat, dimSize, 0, 1);
	       
	       if (_DBG >= 4){
		    // Send results back to cpu memeory. Do we really need this step???
		    cudaMemcpy(C, gpuC, sizeof(float)*BN, cudaMemcpyDeviceToHost);
		    printf("GPUSampling, after samplingKernel iteration %d. C: \n", iter);
		    show_mat(C, N1, "matrix C");
	       }
	  }
	  
	  if (_DBG >= 4){
	       printf("GPUSampling, all sampling iteration done.C: \n");
	       show_mat(C, N1, "matrix C");
//	  FieldViewer(C, N1, N2, fwdMap1, fwdMap2, "GPUSampling: Sampled C");
	       
	  }
     } // end of DO_SAMPLING

     // Now compute mean field.
     for (iter = 1; iter <= MEANFIELDITER; iter++){
	  if (_DBG >= 1 & iter%1 == 0){
	       printf("GPUSampling, mean filed iteration: %d...\n", iter);
	  }
	  SamplingKernel<<<dimGrid, dimBlock>>>(gpuC, gpuIm, gpuFwdMap1, gpuFwdMap2,
						gpuBwdMap1, gpuBwdMap2, par, N1, N2,
						TL,
						gpuRandMat, dimSize, 1, 0);

	  cudaThreadSynchronize();
	  SamplingKernel<<<dimGrid, dimBlock>>>(gpuC, gpuIm, gpuFwdMap1, gpuFwdMap2,
						gpuBwdMap1, gpuBwdMap2, par, N1, N2,
						TL,
						gpuRandMat, dimSize, 1, 1);

	  if (_DBG >= 4){
	       /* Send results back to cpu memeory */
	       cudaMemcpy(C, gpuC, sizeof(float)*BN, cudaMemcpyDeviceToHost);
	       show_mat(C, N1, "expcted value C");
	       FieldViewer(C, N1, N2, fwdMap1, fwdMap2, "GPUSampling: Mean Field");
	  }

     }
     cudaMemcpy(C, gpuC, sizeof(float)*BN, cudaMemcpyDeviceToHost);
     
     /* clean up. */
     cudaFree(gpuC);
     cudaFree(gpuFwdMap1);
     cudaFree(gpuFwdMap2);
     cudaFree(gpuBwdMap1);
     cudaFree(gpuBwdMap2);
     cudaFree(gpuIm);
     free(randMat);

     if (DO_SAMPLING) {
	  cudaFree(gpuRandMat);
     }
}

/* Kernel function  */
__global__ void SamplingKernel(float* C, 
			       const float* im,
			       const unsigned int* gpuFwdMap1,
			       const unsigned int* gpuFwdMap2,
			       const unsigned short* gpuBwdMap1,
			       const unsigned short* gpuBwdMap2,
			       paraType par,
			       unsigned int N1,
			       unsigned int N2,
			       unsigned short TL, // time series length.
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

     uint n1 = blockIdx.x*blockDim.x + threadIdx.x;
     uint n2 = blockIdx.y*blockDim.y + threadIdx.y;

     float gc1 = 1/(sqrt(2*PI*par.sigma21));
     float gc2 = 1/(sqrt(2*PI*par.sigma22));
     float sum = 0, denom = 0, pr1 = 0, pr2 = 0, postsum = 0;
     float lh1 = 0, lh2 = 0, post1 = 0, post2 = 0;
     float upp = 0, upn = 0; // posterior energy for x = 1 (upp) and x = -1 (upn).

     // thread fall outside of matrix C, or mask is zero.
     if (n1 >= N1 | n2 >= N2 | n1 > n2) {return;}

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

     float yi = GetCorr(im, N1, TL, n1, n2);
     upn = par.beta*sum + (yi-par.mu1)*(yi-par.mu1)/(2*par.sigma21) + log(sqrt(par.sigma21));
     upp = - par.beta*sum + (yi-par.mu2)*(yi-par.mu2)/(2*par.sigma22) + log(sqrt(par.sigma22));
     if ((upn-upp) > MAXEXP) {
	  post1 = 0;
	  post2 = 1;
     }
     else if ((upn-upp) < -MAXEXP) {
	  post1 = 1;
	  post2 = 0;
     }
     else {
	  post1 = 1/(1 + exp(upn-upp));
	  post2 = 1 - post1;
     }


     if (job == 1){ // expectation. mean field.
	  GTI(C, n1, n2, N2) = post2 - post1;
     }
     else if (job == 0){ // sampling
	  if (GTI(gpuRandMat, n1, n2, N2) < post1) {
	       GTI(C, n1, n2, N2) = -1;
	  }

	  else{
	       GTI(C, n1, n2, N2) = 1;
	  }
     }
     else if (job == 2){ // compute posterior. p(c == -1)
	  GTI(gpuRandMat, n1, n2, N2) = post2;
//	  if (n1 == n2) 	  GTI(gpuRandMat, n1, n2, N2) = 1;
     }
     else{
	  // something must be wrong. Just exit.
     }
}


// compute correlation.
__device__ float GetCorr(const float* im, 
			 unsigned int N, 
			 unsigned int M,  // time series length
			 int idx1, 
			 int idx2)
{
     int m = 0;
     float mean1 = 0;
     float mean2 = 0;
     float std1 = 0;
     float std2 = 0;
     int lidx1, lidx2;
     float r = 0; // temp variable for sample correlation.
     if (idx1 == idx2){
	  return DIAGCORR;
     }

     // mean of vectors.
     for (m = 0; m < M; m++){
	  lidx1 = idx1*M + m;
	  lidx2 = idx2*M + m;
	  mean1 = mean1 + *(im + lidx1)/M;
	  mean2 = mean2 + *(im + lidx2)/M;
     }
     /* Standard deviation. */
     for (m = 0; m < M; m++){
	  lidx1 = idx1 * M + m;
	  lidx2 = idx2 * M + m;
	  std1 = std1 + pow((*(im+lidx1) - mean1), 2)/(M-1);
	  std2 = std2 + pow((*(im+lidx2) - mean2), 2)/(M-1);
     }
     std1 = sqrt(std1);
     std2 = sqrt(std2);
     /* Sample Correlation. */
     if (std1 == 0 | std2 == 0){
	  return 0;
     }
     else{
	  for (m = 0; m < M; m++){
	       lidx1 = idx1 * M + m;
	       lidx2 = idx2 * M + m;
	       r = r + (*(im + lidx1) - mean1) * (*(im + lidx2)-mean2)/((M-1)*std1*std2);
	  }
	  return r;
     }
}

