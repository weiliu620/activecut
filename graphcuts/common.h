#include <stdio.h>
#include <cstddef>
#include <new>
#include <iostream>
#include <iterator>
#include <string>
#include <fstream>
#include <Eigen/Dense>

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

#include <itkNiftiImageIO.h>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "GCoptimization.h"

using Eigen::MatrixXd;
using Eigen::MatrixXf;

typedef boost::mt19937 twister_base_gen_type;

