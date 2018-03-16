
#include <iostream>
#include <string>

#include <itkAnalyzeImageIO.h>
#include <itkNiftiImageIO.h>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <boost/math/distributions/normal.hpp> // for normal_distribution.
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include "ConfigFile.h"

typedef float          PixelType;
typedef itk::Image< PixelType, 4 >   ImageType;
typedef itk::Image< PixelType, 3 >   MaskImType;
typedef itk::ImageFileReader< ImageType >  fmriReaderType;
typedef itk::ImageFileReader< MaskImType>  maskReaderType;
typedef itk::ImageFileWriter< MaskImType > WriterType;

typedef boost::minstd_rand base_generator_type;
typedef boost::normal_distribution<> distribution_type;
typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type;

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

unsigned int GetRegionSize(maskReaderType::Pointer maskReader,
			   float grayThresh,
			   int xmin, int xmax,
			   int ymin, int ymax,
			   int zmin, int zmax);

void ExtractGray(fmriReaderType::Pointer fmriReader,
		maskReaderType::Pointer maskReader,
		PixelType* im1,
		unsigned int* fwdMap,
		unsigned short* bwdMap,
		float grayThresh,
		 int xmin, int xmax,
		 int ymin, int ymax,
		 int zmin, int zmax);


void GPUCorr(float* R, const float* im1, const float* im2, int M, int N1, int N2);

void GPUInit(float* C, // initial connectivity matrix.
	     const float* R, // correlation (or other affinity) matrix, from data.
	     paraType par,
	     unsigned int regioinSize1, 
	     unsigned int regionSize2);
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
	    const unsigned short* size);

void GPUInit(float* C, // initial connectivity matrix.
	     const float* im,
	     paraType par,
	     unsigned int regionSize1, 
	     unsigned int regionSize2, 
	     unsigned short TL); // time series length.

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
		 const unsigned short* size);


void GPUEst(const float* C,
	    const unsigned int* fwdMap1,
	    const unsigned int* fwdMap2,
	    const unsigned short* bwdMap1,
	    const unsigned short* bwdMap2,
	    paraType* par,
	    unsigned int N1,
	    unsigned int N2,
	    const unsigned short* size);

void ExtractRealData(float* & im1,
		     float*  & im2,
		     unsigned int*  & fwdMap1,
		     unsigned int*  & fwdMap2,
		     unsigned short*  & bwdMap1,
		     unsigned short*  & bwdMap2,
		     unsigned int& regionSize1,
		     unsigned int& regionSize2,
		     unsigned short* size,
		     ConfigFile config);

void ExtractToyData(float* & im1,
		    float*  & im2,
		    unsigned int*  & fwdMap1,
		    unsigned int*  & fwdMap2,
		    unsigned short*  & bwdMap1,
		    unsigned short*  & bwdMap2,
		    unsigned int& regionSize1,
		    unsigned int& regionSize2,
		    unsigned short* size,
		    float snr);

int fltktest();

int SaveResults(const float* post,
		unsigned short N1,
		unsigned short N2,
		const unsigned int* fwdMap1,
		const unsigned int* fwdMap2,
		ConfigFile config);

int GPUDetrend(float* im, unsigned int N, unsigned int D);
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
		  float* post2);

int SaveResults(const unsigned short* seed,
		unsigned short targetSlice, // number of slice for saving.
		const float* post,
		unsigned short N1,
		unsigned short N2,
		const unsigned int* fwdMap1,
		const unsigned int* fwdMap2,
		const char datadir[]); // voume size. 
void MLEstep(float* R,
	     paraType par,
	     unsigned int N1,
	     unsigned int N2,
	     unsigned int Nk1,
	     unsigned int Nk2,
	     float* post2);

int ThreshCorr(const float*R, 
	       float* post,
	       unsigned int N1,
	       unsigned int N2,
	       const unsigned int* fwdMap1,
	       const unsigned int* fwdMap2,
	       unsigned short* size,
	       ConfigFile config);

