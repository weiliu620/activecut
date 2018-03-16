#ifndef __COMMON_H__
#define __COMMON_H__

#include <Eigen/Dense>
#include <Eigen/Core>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <new>
#include <iostream>
#include <iterator>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>

/* #include <boost/math/distributions/normal.hpp> // for normal_distribution. */
/* #include <boost/math/distributions/beta.hpp> */
/* #include <boost/math/special_functions/bessel.hpp> */
/* #include <boost/random/linear_congruential.hpp> */
/* #include <boost/random/uniform_int.hpp> */
/* #include <boost/random/uniform_real.hpp> */
/* #include <boost/random/normal_distribution.hpp> */

/* #include <boost/random/variate_generator.hpp> */
/* #include <boost/random/uniform_on_sphere.hpp> */
/* #include <boost/random/gamma_distribution.hpp> */
/* #include <boost/random.hpp> */
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

/* #include <itkNiftiImageIO.h> */
/* #include "itkImage.h" */
/* #include "itkImageFileReader.h" */
/* #include "itkImageFileWriter.h" */
/* #include "itkImageRegionIterator.h" */
/* #include "itkImageRegionConstIterator.h" */
/* #include "itkImageRegionIteratorWithIndex.h" */
/* #include "itkExtractImageFilter.h" */
/* #include "itkNeighborhoodIterator.h" */
/* #include "itkNeighborhoodIterator.h" */
/* #include "itkConstantBoundaryCondition.h" */
/* #include "itkImageSeriesReader.h" */
/* #include "itkVectorImage.h" */
/* #include "itkAddConstantToImageFilter.h" */
/* #include "itkMaskImageFilter.h" */
/* #include "itkMinimumMaximumImageCalculator.h" */

namespace po = boost::program_options;

struct Codescent_parameters {
     float finalLambda; 
     float alpha;
     float initLambda;
     float lambda;
     unsigned numSteps;
     float stopThr;
     unsigned numObs;
     unsigned numVar;
     float eps;
};



#endif
