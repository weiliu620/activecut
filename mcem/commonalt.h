#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdio.h>
#include <cmath>
#include <new>
#include <iostream>
#include <iterator>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include <boost/math/distributions/normal.hpp> // for normal_distribution.
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
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

#include <itkAnalyzeImageIO.h>
#include <itkNiftiImageIO.h>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkExtractImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"


namespace po = boost::program_options;
typedef boost::mt19937 twister_base_gen_type;

typedef itk::Image<float, 2> ImageType2D;
typedef itk::Image<short int, 3> ImageType3D;
typedef itk::Image<unsigned char, 4> Image4DChar;
typedef itk::ImageFileReader< ImageType2D >  ReaderType;
typedef itk::ImageFileWriter< ImageType2D >  WriterType;


typedef itk::ImageFileReader< ImageType3D >  ReaderType3D;
typedef itk::ImageFileWriter< ImageType3D >  WriterType3D;
typedef itk::ImageRegionConstIterator< ImageType2D > ConstIteratorType2D;
typedef itk::ImageRegionConstIterator< ImageType3D > ConstIteratorType3D;

typedef itk::ImageRegionIterator< ImageType2D>       IteratorType2D;
typedef itk::ImageRegionIterator< ImageType3D>       IteratorType3D;
typedef itk::ImageRegionIterator< Image4DChar >      IteratorType4D;
struct  CompType
{
     float mu;
     float sigma;
     int label;
     long numPoints;
};

int saveimage2d(ImageType2D::Pointer ip, std::string fname);
int saveimage3d(ImageType3D::Pointer ip, std::string fname);
int saveExtractSlice(ImageType3D::Pointer inPtr,
		     unsigned sliceIdx, std::string fname);

double singleMCSampleLL(ImageType2D::Pointer imagePtr,
			ImageType3D::Pointer samplePtr,
			unsigned short mcSampleIdx,
			float beta,
			CompType * const cls,
			unsigned short numClusters);
#endif
