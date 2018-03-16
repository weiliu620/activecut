#include "cppheader.h"
#include "common.h"

void EMPost(float* C, 
	    const float* R,
	    const float* im,
	    const unsigned int* fwdMap1,
	    const unsigned int* fwdMap2,
	    const unsigned short* bwdMap1,
	    const unsigned short* bwdMap2,
	    unsigned int regionSize1, 
	    unsigned int regionSize2,
	    unsigned short TL, // time series length.
	    const unsigned short* size)
{
     ConfigFile config(PARAFILE);
     paraType para;
     paraType old_para;
     para.mu1 = 0;
     para.mu2 = 0.1;
     para.sigma21 = 0.05;
     para.sigma22 = 0.05;
     para.T = 1;
     para.beta = INITBETA;
     float post2 = 0;

     unsigned annIter = 0;
     int _DEBUG = config.read<int>("_DEBUG", 0);
     if (_DEBUG >= 1){
	  printf("parameters initialized:\n");
	  printf("mu1 = %5.2f\nmu2 = %5.2f\nsigma21 = %2.2f\nsigma22 = %2.2f\nbeta = %f\n",
		 para.mu1, para.mu2, para.sigma21, para.sigma22, para.beta);
     }
     
     int i, j;
     unsigned int EMIter = 0;
          
     // Init connectivity matrix C.
     GPUInit(C, im, para, regionSize1, regionSize2, TL);
     if (_DBG >= 3){
	  printf("EMPost.cxx. After GPUInit. C: \n");
	  show_mat(C, regionSize1, "matrix C");
	  FieldViewer(C, regionSize1, regionSize2, fwdMap1, fwdMap2, "Init'd C");
     }
     
     float change = 1e10;
     unsigned int n = 0, n1 = 0, n2 = 0;

     // EM to estimate posterior and parameters.
     // init estimate of Nk1 and Nk2
     float  Nk1, Nk2;
     for (n = 0; n < regionSize1 * (regionSize1+1)/2; n++){
	  Nk2 = Nk2 + (C[n] == 1);
     }
     Nk1 = regionSize1 * (regionSize1-1)/2 - Nk2;
     
     while (change > EPS){
	  EMIter++;
	  printf("-----------EM iteration %d.  E step.---------------------------------------\n", EMIter);
	  printf("Call GPU sampling now...\n");
	  // randMat used for random number when sampling, and used for 
	  // posterior p(c = 1) when computing posterior, i.e. 
	  // job = "post".
	  if (USEMRF == 1){
	  GPUSampling(C, im, fwdMap1, fwdMap2, bwdMap1, bwdMap2, 
		      para, regionSize1, regionSize2, TL, size);
	  }
	  else{

//	       MLEstep(R, para, regionSize1, regionSize2, Nk1, Nk2, post2);
	  }

	  if (_DBG >= 3){
	       printf("Mean field all iteration done. C:\n");
	       show_mat(C, regionSize1, "expcted value C");
	       FieldViewer(C, regionSize1, regionSize1, fwdMap1, fwdMap2, "GPUSampling: Mean Field");
	  }

	  // M step, estimate parameters
	  printf("---------- EM iteration %d. M step. ---------------------------------------\n", EMIter);
	  // save parameters to old_para.
	  old_para = para;
	  Nk1 = 0;
	  Nk2 = 0;
	  // estimate Nk. Do not count diag element, as they do not obey the likelihood model.
	  for (n1 = 0; n1 < regionSize1; n1++) {
	       for (n2 = n1+1; n2 < regionSize1; n2++) {
		    post2 = (1 + GTI(C, n1, n2, regionSize1))/2;
		    Nk2 = Nk2 + post2;
	       }
	  }
	  // all elements in upper triangular minus Nk2.
	  Nk1 = regionSize1 * (regionSize1-1)/2 - Nk2;
	  // also estimate mu1?
	  if (config.read("est_mu1", 1)) {
	       para.mu1 = 0;
	       for (n1 = 0; n1 < regionSize1; n1++) {
		    for (n2 = n1+1; n2 < regionSize1; n2++) {
			 post2 = (1 + GTI(C, n1, n2, regionSize1))/2;
			 para.mu1 = para.mu1 + (1- post2) * GTI(R, n1, n2, regionSize1)/Nk1;
		    }
	       }
	  }

	  // estimate mu2
	  para.mu2 = 0;
	  for (n1 = 0; n1 < regionSize1; n1++) {
	       for (n2 = n1+1; n2 < regionSize1; n2++) {
		    post2 = (1 + GTI(C, n1, n2, regionSize1))/2;		    
		    para.mu2 = para.mu2 + post2 * GTI(R, n1, n2, regionSize1)/Nk2;
	       }
	  }
	  // estimate sigma21 and sigma22
	  para.sigma21 = 0;
	  para.sigma22 = 0;

	  for (n1 = 0; n1 < regionSize1; n1++) {
	       for (n2 = n1+1; n2 < regionSize1; n2++) {
		    post2 = (1 + GTI(C, n1, n2, regionSize1))/2;		    
		    para.sigma21 = para.sigma21 + 
			 (1 - post2) * (GTI(R, n1, n2, regionSize1) - para.mu1) * (GTI(R, n1, n2, regionSize1) - para.mu1)/Nk1;
		    para.sigma22 = para.sigma22 + 
			 post2 * (GTI(R, n1, n2, regionSize1) - para.mu2) * (GTI(R, n1, n2, regionSize1) - para.mu2)/Nk2;
	       }
	  }

	  if (USEMRF == 1){
	  // Estimate beta. This need lots of work. We begin with beta
	  // in previous iteration as initital value for Nettonw's method.
	  GPUEst(C, fwdMap1, fwdMap2, bwdMap1, bwdMap2, &para, regionSize1, regionSize2, size);
	  }
	  change = fabs(para.mu1 - old_para.mu1) + 
	       fabs(para.mu2 - old_para.mu2) + 
	       fabs(para.sigma21 - old_para.sigma21) + 
	       fabs(para.sigma22 - old_para.sigma22) +
	       fabs(para.beta - old_para.beta);
	  if (_DBG >= 2){
	       printf("EMPost, Nk1 = %f, Nk2 = %f, mu1 = %f, mu2 = %f,\n sigma21 = %f, sigma22 = %f\n",
		      Nk1, Nk2, para.mu1, para.mu2, para.sigma21, para.sigma22);
	  }

     }
     
     if (_DBG >=3){
	  printf("EMPost, mean filed:\n");
	  show_mat(C, regionSize1, "mean field matrix:");
	  FieldViewer(C, regionSize1, regionSize2, 
			   fwdMap1, fwdMap2, "posterior map");
     }
}
