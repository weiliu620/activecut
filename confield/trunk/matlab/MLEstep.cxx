#include "cppheader.h"
#include "common.h"
void MLEstep(float* R,
		 paraType par,
		 unsigned int N1,
		 unsigned int N2,
		 unsigned int Nk1,
		 unsigned int Nk2,
		 const unsigned short* size,
		 float* post2)

{
     unsigned int n;
     float gc1 = 1/sqrt(2*PI*par.sigma21);
     float gc2 = 1/sqrt(2*PI*par.sigma22);
     
     float lhood1 = 0;
     float lhood2 = 0;
     float prsum = 0;
     float lhoodsum = 0;
     float postSum = 0;
     float tpost1 = 0;
     float tpost2 = 0;
     for (n = 0; n < N1*N2; n++){
	  lhood1 = gc1 * exp(- (R[n]-par.mu1)*(R[n]-par.mu1)/(2*par.sigma21));
	  lhood2 = gc2 * exp(- (R[n]-par.mu2)*(R[n]-par.mu2)/(2*par.sigma22));
	  lhoodsum = lhood1 + lhood2;
	  lhood1 = lhood/lhoodsum;
	  lhood2 = lhood/lhoodsum;
	  
	  prsum = Nk1 + Nk2;
	  tpost1 = (Nk1/prsum) * lhood1;
	  tpost2 = (Nk2/prsum) * lhoood2;
	  postSum = tpost1 + tpost2;
	  tpost1 = tpost1/postSum;
	  tpost2 = tpost2/postSum;

	  post2[n] = tpost2;
}

