
#include <stdio.h>
#include <cmath>
#include <new>
#include <iostream>
#include <iterator>
#include <string>
#include <fstream>

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

#include <Eigen/Core>

using namespace std;
namespace po = boost::program_options;
typedef boost::mt19937 twister_base_gen_type;
USING_PART_OF_NAMESPACE_EIGEN


typedef itk::Image<float, 2> ImageType2D;
typedef itk::Image<int, 3> ImageType3D;
typedef itk::ImageFileReader< ImageType2D >  ReaderType;
typedef itk::ImageFileWriter< ImageType2D >  WriterType;

typedef itk::ImageFileReader< ImageType3D >  ReaderType3D;
typedef itk::ImageFileWriter< ImageType3D >  WriterType3D;

typedef struct {
     float mu;
     float sigma;
     int label;
     long numPoints;
} CompType;

int SaveToy(MatrixXf image, string fname);
int saveimage2d(ImageType2D::Pointer ip, string fname);
int saveimage3d(ImageType3D::Pointer ip, string fname);


int est_ll(ImageType3D::Pointer samplePtr,
	   ImageType2D::Pointer imagePtr,
	   CompType* cls,
	   int numClusters);

int printcls(CompType* cls, int nLabels);

int est_ll(ImageType3D::Pointer samplePtr,
	   ImageType2D::Pointer imagePtr,
	   CompType* cls,
	   int numClusters);

int ll_over_beta(ImageType3D::Pointer samplePtr,
		 int numClusters);

int estep(ImageType2D::Pointer imagePtr,
	  ImageType3D::Pointer samplePtr,
	  CompType* cls,
	  int numClusters,
	  int burnin,
	  int numScan,
	  double beta);

double eval_ll(ImageType3D::Pointer samplePtr,
		  int numClusters,
		  double beta);

double eval_cll(ImageType3D::Pointer samplePtr,
		ImageType2D::Pointer imagePtr,
		CompType* cls);

double golden_section_search(double a, 
			     double b,
			     ImageType3D::Pointer samplePtr,
			     int numClusters);

double descent(ImageType3D::Pointer samplePtr,
	       int numClusters,
	       double beta,
	       double t);
