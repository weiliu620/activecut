#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <cuda.h>

#ifndef PARATYPE
#define PARATYPE
typedef struct {
     float mu1;
     float mu2;
     float sigma21;
     float sigma22;
     float beta;
     float T;
} paraType;
#endif


#define _DBG 2
#define EPS 0.0001 // small number to see if parameters change.

//#define MAXBETA 15 // max beta.
#define MAXEMITER 50 // EM max iteration.
#define MAXANNITER 50 // EM max iteration.
#define FINALBETA 0.5
#define PI 3.14159265
#define INITTPRAUTRE 5
#define INITBETA 0.01

#define RMIN 0 // output matrix row min.
#define RMAX 10 // output matrix row max.

#define CMIN 0 // output matrix col min.
#define CMAX 10 // output matrix col max.

#define GIBBSMAXITER 10 // Gibbs max iteration.
#define MEANFIELDITER 10 // mean field max iteration.
#define CUEPS 0.0001
#define CUMAXBETA 10
#define MAXEXP 12 // max exponential for computing probability in GPUSampling.
// for debugging of ShowSlice.
#define SEED0 34
#define SEED1 18
#define SEED2 27
#define TARGETSLICE 27

// if we're testing toy example. decide how we visualize data.
#define _TOY 0

// use mrf or not.
#define USEMRF 1
#define THCORRMAP 0

// parameter file name
#define PARAFILE "parameters.inp"

// Forward Declaration.
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
#define MYTEST 1
//#define FORCE_ANNEALING 0

void checkCUDAError(const char *msg);
#define _EMU 0

int FieldViewer(const float* post, 
		unsigned short N1, 
		unsigned short N2, 
		const unsigned int* fwdMap1,
		const unsigned int* fwdmap2,
		char* title);

int ToyViewer(const float* post,
	      unsigned short N1,
	      unsigned short N2);


