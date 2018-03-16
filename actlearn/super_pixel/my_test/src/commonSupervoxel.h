//////////////////////////////
//
// for ITK filters
//
//////////////////////////////
// basic image operations
#include "itkImage.h"
// basic image files read and write
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
// don't know if necessary
#include "itkFixedArray.h"


//////////////////////////////
//
// standard c++
//
//////////////////////////////
#include <stdlib.h>     /* div, div_t */
#include <iostream>

#include <stddef.h>


//////////////////////////////
//
// for using boost
//
//////////////////////////////

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
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>

