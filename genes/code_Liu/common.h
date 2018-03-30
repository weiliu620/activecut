#include <pthread.h>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cassert>
#include <new>
#include <iostream>
#include <iterator>
#include <string>
// #include <array>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/version.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/dynamic_bitset.hpp>

#include <unistd.h>

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
#include "itkImageSeriesReader.h"
#include "itkVectorImage.h"
#include "itkAddConstantToImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include <itkCastImageFilter.h>

namespace po = boost::program_options;

typedef itk::Image<float, 2> ImageType2D;
typedef itk::Image<short int, 3> ImageType3D;
typedef itk::Image<unsigned, 3> ImageType3DU;
typedef itk::Image<short int, 3> ImageType3DShort;
typedef itk::Image<char, 3> ImageType3DChar;
typedef itk::Image<float, 3> ImageType3DFloat;
typedef itk::Image<int, 4> ImageType4Int;
typedef itk::Image<short, 4> ImageType4DS;
typedef itk::Image<float, 4> ImageType4DF;

typedef itk::VectorImage<float, 3> ImageType3DVecF;
typedef itk::VectorImage<unsigned short, 3> ImageType3DVecUS;
typedef itk::VectorImage<bool, 3> ImageType3DVecB;
typedef itk::VectorImage<bool, 4> ImageType4DVecB;
typedef itk::Image<unsigned short, 4> ImageType4DUS;
typedef itk::Image< vnl_vector< float >, 3 > VnlVectorImageType;

typedef itk::Image<char, 4> Image4DChar;
typedef itk::Image<char, 4> ImageType4DChar;
typedef itk::Image<float, 4> ImageType4DFloat;
typedef itk::Image<float, 4> Image4DFloat;
typedef itk::Image<char, 5> ImageType5DChar;
typedef itk::Image<char, 6> ImageType6DChar;
typedef itk::Image<float, 5> ImageType5DFloat;
typedef itk::Image<std::vector< std::vector<double> >, 3 > ImageTypeMat;

typedef itk::ImageFileReader< ImageType2D >  ReaderType;
typedef itk::ImageFileReader< ImageType3D >  ReaderType3D;
typedef itk::ImageFileReader< ImageType3DShort >  ReaderType3DShort;
typedef itk::ImageFileReader< ImageType3DChar >  ReaderType3DChar;
typedef itk::ImageFileReader< ImageType3DFloat >  ReaderType3DFloat;
typedef itk::ImageFileReader< ImageType4DChar >  ReaderType4DChar;
typedef itk::ImageFileReader< ImageType4DFloat >  ReaderType4DFloat;
typedef itk::ImageFileReader< ImageType4DS >  ReaderType4DS;

typedef itk::ImageSeriesReader< ImageType4DFloat >  SeriesReaderType4DFloat;
typedef itk::ImageSeriesReader< ImageType5DFloat >  SeriesReaderType5DFloat;

typedef itk::ImageFileWriter< ImageType2D >  WriterType;
typedef itk::ImageFileWriter< ImageType3D >  WriterType3D;
typedef itk::ImageFileWriter< ImageType3DShort >  WriterType3DShort;
typedef itk::ImageFileWriter< ImageType3DFloat >  WriterType3DFloat;
typedef itk::ImageFileWriter< ImageType3DChar >  WriterType3DChar;
typedef itk::ImageFileWriter< ImageType4DChar >  WriterType4DChar;
typedef itk::ImageFileWriter< ImageType4DFloat >  WriterType4DFloat;
typedef itk::ImageFileWriter< ImageType3DVecF >  WriterType3DVecF;
typedef itk::ImageFileWriter< ImageType4DUS >  WriterType4DUS;
typedef itk::ImageFileWriter< ImageType4Int >  WriterType4Int;
typedef itk::ImageFileWriter< ImageType4DS >  WriterType4DS;


typedef itk::ImageRegionConstIterator< ImageType2D > ConstIteratorType2D;
typedef itk::ImageRegionConstIterator< ImageType3D > ConstIteratorType3D;
typedef itk::ImageRegionConstIterator< ImageType3DChar > ConstIteratorType3DChar;

typedef itk::ImageRegionIterator< ImageType2D>       IteratorType2D;
typedef itk::ImageRegionIterator< ImageType3D>       IteratorType3D;
typedef itk::ImageRegionIterator< ImageType3DChar>       IteratorType3DChar;
typedef itk::ImageRegionIterator< ImageType3DShort>       IteratorType3DShort;
typedef itk::ImageRegionIteratorWithIndex< ImageType3DChar>       IteratorType3DCharIdx;
typedef itk::ImageRegionIterator< ImageType3DFloat>       IteratorType3DFloat;
typedef itk::ImageRegionIteratorWithIndex< ImageType3DChar>       IteratorType3DCharIdx;
typedef itk::ImageRegionIterator< Image4DChar >      IteratorType4DChar;
typedef itk::ImageRegionIterator< Image4DFloat >      IteratorType4DFloat;
typedef itk::ImageRegionIterator< ImageType5DFloat >      IteratorType5DFloat;
typedef itk::ImageRegionIterator< ImageType3DVecF >       IteratorType3DVecF;
typedef itk::ImageRegionIterator< ImageType3DVecUS >      IteratorType3DVecUS;
typedef itk::ImageRegionIterator< VnlVectorImageType >      IteratorTypeVnlVector;
typedef itk::ImageRegionIterator< ImageType4DS >      IteratorType4DS;

typedef itk::AddConstantToImageFilter <ImageType3DChar, unsigned char, ImageType3DChar> AddConstantToImageFilterType;
