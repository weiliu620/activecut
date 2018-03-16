#include "cppheader.h"
#include "common.h"

void EMPost(float* C, 
	    float* R, 
	    float* post2,
	    const unsigned int* fwdMap1,
	    const unsigned int* fwdMap2,
	    const unsigned short* bwdMap1,
	    const unsigned short* bwdMap2,
	    unsigned int regionSize1, 
	    unsigned int regionSize2,
	    const unsigned short* size)
{
     ConfigFile config(PARAFILE);
     paraType para;
     paraType old_para;
     para.mu1 = 0;
     para.mu2 = 0.2;
     para.sigma21 = 0.05;
     para.sigma22 = 0.05;
     para.T = 1;
     para.beta = INITBETA;
     // one pixel clique coefficient alpha.
     para.alpha = config.read<float>("alpha", 0);

     unsigned annIter = 0;
     if (_DBG >= 1){
	  printf("parameters initialized:\n");
	  printf("mu1 = %5.2f\nmu2 = %5.2f\nsigma21 = %2.2f\nsigma22 = %2.2f\nbeta = %f\nalpha = %2.2f\n",
		 para.mu1, para.mu2, para.sigma21, para.sigma22, para.beta, para.alpha);
     }
     
     int i, j;
     unsigned int EMIter = 0;

     // Init connectivity matrix C.
     GPUInit(C, R, para, regionSize1, regionSize2);
     if (_DBG >= 3){
	  printf("EMPost.cxx. After GPUInit. C: \n");
	  for (i = RMIN; i <= RMAX; i++){
	       for (j = CMIN; j <= CMAX; j++){
		    printf("[%d][%d]=%1.1f ", i, j, C[i*regionSize2 + j]);
	       }
	       printf("\n");
	  }
     }
          
     
     float change = 1e10;
     unsigned int n = 0, n1 = 0, n2 = 0;

     // EM to estimate posterior and parameters.
     // init estimate of Nk1 and Nk2
     float  Nk1, Nk2;
     for (n1 = 0; n1 < regionSize1; n1++){
	  for (n2 = 0; n2 < regionSize2; n2++) {
	       if (n1 != n2) {
		    Nk2 = Nk2 + (C[n1+regionSize2+n2] > 0);
	       }
	  }
     }

     Nk1 = regionSize1 * regionSize2 - Nk2 - regionSize1;
     
     while (change > EPS){
	  EMIter++;
	  printf("-----------EM iteration %d.  E step.---------------------------------------\n", EMIter);

	  printf("Call GPU sampling now...\n");
	  // randMat used for random number when sampling, and used for 
	  // posterior p(c = 1) when computing posterior, i.e. 
	  // job = "post".
	  if (USEMRF == 1){
	  GPUSampling(C, R, fwdMap1, fwdMap2, bwdMap1, bwdMap2, 
		      para, regionSize1, regionSize2, size, post2);
	  }
	  else{

	       MLEstep(R, para, regionSize1, regionSize2, Nk1, Nk2, post2);
	  }


	  // M step, estimate parameters
	  printf("---------- EM iteration %d. M step. ---------------------------------------\n", EMIter);
	  // save parameters to old_para.
	  old_para = para;
	  Nk1 = 0;
	  Nk2 = 0;
	  for (n1 = 0; n1 < regionSize1; n1++){
	       for (n2 = 0; n2 < regionSize2; n2++) {
		    if (n1 != n2) {
			 Nk2 = Nk2 + post2[n1*regionSize2+n2];
		    }
	       }
	  }
	  Nk1 = regionSize1 * regionSize2 - Nk2 - regionSize1;
	  // also estimate mu1?
	  if (config.read("est_mu1", 1)) {
	       para.mu1 = 0;
	       for (n1 = 0; n1 < regionSize1; n1++){
		    for (n2 = 0; n2 < regionSize2; n2++) {
			 if (n1 != n2) {
			      para.mu1 = para.mu1 + (1- post2[n1*regionSize2+n2]) * R[n1*regionSize2+n2]/Nk1;
			 }
		    }
	       }	  
	  }

	  // estimate mu2
	  para.mu2 = 0;
	  for (n1 = 0; n1 < regionSize1; n1++){
	       for (n2 = 0; n2 < regionSize2; n2++) {
		    if (n1 != n2) {
			 para.mu2 = para.mu2 + (post2[n1*regionSize2+n2]) * R[n1*regionSize2+n2]/Nk2;
		    }
	       }
	  }
	  // estimate sigma21 and sigma22
	  para.sigma21 = 0;
	  para.sigma22 = 0;
	  for (n1 = 0; n1 < regionSize1; n1++){
	       for (n2 = 0; n2 < regionSize2; n2++) {
		    if (n1 != n2) {
			 para.sigma21 = para.sigma21 + 
		    (1 - post2[n1*regionSize2+n2]) * (R[n1*regionSize2+n2] - para.mu1)*(R[n1*regionSize2+n2] - para.mu1)/Nk1;
			 para.sigma22 = para.sigma22 + 
			      post2[n1*regionSize2+n2] * (R[n1*regionSize2+n2] - para.mu2)*(R[n1*regionSize2+n2] - para.mu2)/Nk2;
		    }
	       }
	  }
	  if (USEMRF == 1){
	  // Estimate alpha and beta. This need lots of work. We begin with beta
	  // in previous iteration as initital value for Nettonw's method.
	  GPUEst(C, fwdMap1, fwdMap2, bwdMap1, bwdMap2, &para, regionSize1, regionSize2, size);
	  }
	  change = fabs(para.mu1 - old_para.mu1) + 
	       fabs(para.mu2 - old_para.mu2) + 
	       fabs(para.sigma21 - old_para.sigma21) + 
	       fabs(para.sigma22 - old_para.sigma22) +
	       fabs(para.beta - old_para.beta) + 
	       fabs(para.alpha - old_para.alpha);
	  if (_DBG >= 2){
	       printf("EMPost, Nk1 = %f, Nk2 = %f, mu1 = %f, mu2 = %f,\n sigma21 = %f, sigma22 = %f\n",
		      Nk1, Nk2, para.mu1, para.mu2, para.sigma21, para.sigma22);
	  }

	  
	  if (_DBG >=4){
	       printf("EMPost, posterior p(c = 1):\n");
	       for (i = RMIN; i <= RMAX; i++){
		    for (j = CMIN; j <= CMAX; j++){
			 printf("[%d][%d]=%3.2f ", i, j, post2[i*regionSize2 + j]);
		    }
		    printf("\n");
	       }
	       FieldViewer(post2, regionSize1,regionSize2, fwdMap1, fwdMap2, "Posterior");
	  }

     }
     
     if (_DBG >=4){
	  printf("EMPost, posterior p(c = 1):\n");
	  for (i = RMIN; i <= RMAX; i++){
	       for (j = CMIN; j <= CMAX; j++){
		    printf("[%d][%d]=%3.2f ", i, j, post2[i*regionSize2 + j]);
	       }
	       printf("\n");
	  }

     }
     if (_DBG >= 3){
	  if (config.read("toy",0)){
	       ToyViewer(post2, regionSize1, regionSize2);
	  }
	  else {
	       FieldViewer(post2, regionSize1, regionSize2, 
			   fwdMap1, fwdMap2, "posterior map");
	  }
     }
     
/*	  
     // Annealing.
     GPUAnnealing(C, R, fwdMap1, fwdMap2, bwdMap1, bwdMap2, 
		 para, regionSize1, regionSize2, size, post2);
*/
}
