#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdio.h>
#include <cmath>
#include <cassert>
#include <new>
#include <iostream>
#include <iterator>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include <boost/math/distributions/normal.hpp> // for normal_distribution.
#include <boost/math/distributions/beta.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_on_sphere.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random.hpp>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/version.hpp>

/* #include <itkAnalyzeImageIO.h> */
#include <itkNiftiImageIO.h>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkExtractImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkMinimumMaximumImageCalculator.h"


namespace po = boost::program_options;
typedef boost::mt19937 twister_base_gen_type;

typedef itk::Image<float, 2> ImageType2D;
typedef itk::Image<short int, 3> ImageType3D;
typedef itk::Image<char, 3> ImageType3DChar;
typedef itk::Image<float, 3> ImageType3DFloat;

typedef itk::Image<char, 4> Image4DChar;
typedef itk::Image<char, 4> ImageType4DChar;
typedef itk::Image<float, 4> ImageType4DFloat;
typedef itk::Image<float, 4> Image4DFloat;
typedef itk::Image<char, 5> ImageType5DChar;

typedef itk::ImageFileReader< ImageType2D >  ReaderType;
typedef itk::ImageFileReader< ImageType3D >  ReaderType3D;
typedef itk::ImageFileReader< ImageType3DChar >  ReaderType3DChar;
typedef itk::ImageFileReader< ImageType3DFloat >  ReaderType3DFloat;
typedef itk::ImageFileReader< ImageType4DChar >  ReaderType4DChar;
typedef itk::ImageFileReader< ImageType4DFloat >  ReaderType4DFloat;

typedef itk::ImageFileWriter< ImageType2D >  WriterType;
typedef itk::ImageFileWriter< ImageType3D >  WriterType3D;
typedef itk::ImageFileWriter< ImageType3DFloat >  WriterType3DFloat;
typedef itk::ImageFileWriter< ImageType3DChar >  WriterType3DChar;
typedef itk::ImageFileWriter< ImageType4DChar >  WriterType4DChar;
typedef itk::ImageFileWriter< ImageType4DFloat >  WriterType4DFloat;

typedef itk::ImageRegionConstIterator< ImageType2D > ConstIteratorType2D;
typedef itk::ImageRegionConstIterator< ImageType3D > ConstIteratorType3D;

typedef itk::ImageRegionIterator< ImageType2D>       IteratorType2D;
typedef itk::ImageRegionIterator< ImageType3D>       IteratorType3D;
typedef itk::ImageRegionIterator< ImageType3DChar>       IteratorType3DChar;
typedef itk::ImageRegionIterator< ImageType3DFloat>       IteratorType3DFloat;
typedef itk::ImageRegionIteratorWithIndex< ImageType3DChar>       IteratorType3DCharIdx;
typedef itk::ImageRegionIterator< Image4DChar >      IteratorType4DChar;
typedef itk::ImageRegionIterator< Image4DFloat >      IteratorType4DFloat;

struct  CompType
{
     vnl_vector <float> mu;
     float meanNorm;
     float sigma;
     float kappa;
     int label;
     long numPoints;
};

struct  NeighborIndicatorType
{
     bool xplus;
     bool xminus;
     bool yplus;
     bool yminus;
     bool zplus;
     bool zminus;
};

int saveimage2d(ImageType2D::Pointer ip, std::string fname);

int saveimage3d(ImageType3D::Pointer ip, std::string fname);

int save3dchar(ImageType3DChar::Pointer ip, std::string fname);

int saveimage4dchar(ImageType4DChar::Pointer ip, std::string fname);

int saveExtractSlice(ImageType3D::Pointer inPtr,
		     unsigned sliceIdx, std::string fname);

int saveExtractVolume(ImageType4DChar::Pointer inPtr,
		      unsigned volumeIdx, std::string fname);
int generateVonMise(vnl_vector<float>& randomVector, float lambda);

int generateVonMiseWood(vnl_vector<float>& randomVector, float kappa);

int getRotationMat(vnl_matrix<float>& rotationMat,
		   vnl_vector<float>& srcVec, 
		   vnl_vector<float>& destVec);

template <class T>
int saveimage(T ip, std::string fname);

template <class T> int mytest(T var);

double singleSamplePriorLL(ImageType4DChar::Pointer samplePtr,
			   ImageType5DChar::Pointer neighborSumPtr,
			   unsigned short mcSampleIdx,
			   float beta,
			   unsigned short numClusters,
			   unsigned short numScan);

double singleSampleCondLL(ImageType4DFloat::Pointer imagePtr,
			  ImageType4DChar::Pointer samplePtr,
			  unsigned short mcSampleIdx,
			  std::vector<CompType> & cls);

int printVnlVector(vnl_vector<float> vec, unsigned numElements);
int printVnlMatrix(vnl_matrix<float> mat, unsigned numElements);
double logBesselI(float nu, double x);
#endif
