#include "commonalt.h"
#include "/home/sci/weiliu/packages/InsightToolkit-3.20.0/Utilities/vxl/core/vnl/algo/vnl_qr.h"
#include "/home/sci/weiliu/packages/InsightToolkit-3.20.0/Utilities/vxl/core/vnl/vnl_inverse.h"

using std::cout;
using std::endl;

twister_base_gen_type mygenerator(42u);

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     unsigned numSamples = 0;

     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "produce help message")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")
	  ("numsamples,n", po::value<unsigned >(&numSamples)->default_value(100), 
	   "num of smaples.");

     po::options_description cmdline_options;
     cmdline_options.add(generic);

     po::positional_options_description p;
     p.add("in", -1);

     po::variables_map vm;        
     po::store(po::command_line_parser(argc, argv).
	       options(cmdline_options).positional(p).run(), vm);
     po::notify(vm);    

     try {

	  if (vm.count("help")) {
	       std::cout << "Usage: generateimage [options]\n";
	       std::cout << cmdline_options << "\n";
	       return 0;
	  }
	  
     }

     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    


     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     const unsigned timeSeriesLength = 3;
     vnl_vector <float> timeSeries(timeSeriesLength, 0);

     // Compute the transformation matrix to rotate vector such that
     // it's aligned with one axis.

     vnl_matrix<float> hyperBases(timeSeriesLength, timeSeriesLength);
     vnl_matrix<float> naturalBases(timeSeriesLength, timeSeriesLength);
     vnl_matrix<float> transMat(timeSeriesLength, timeSeriesLength);

     // Construct original bases, which is for the space where the
     // hyperplane resides.

     hyperBases.fill(0);
     hyperBases.set_row(timeSeriesLength - 1, -1);
     hyperBases.fill_diagonal(1);
     hyperBases.set_column(timeSeriesLength - 1, 1);

     
     // Construct newBasis, which is the natural bases of destination
     // space.
     naturalBases.fill(0);
     naturalBases.fill_diagonal(1);

     vnl_qr<float> qrComp(hyperBases);
     vnl_matrix<float> orthHyperBases = qrComp.Q();

     // becuse naturalBases = orthHyperBases * transMat, we want to
     // get the matrix transMat. And coordinates_under_hyperBases =
     // transMat * coordinates_under_naturalBases.
     transMat = vnl_inverse(orthHyperBases) * naturalBases;

     
     unsigned sampleIdx = 0;
     for (sampleIdx = 0; sampleIdx < numSamples; sampleIdx ++) {

	  timeSeries[0] = uni();
	  timeSeries[1] = uni();
	  timeSeries[2] = uni();
	  printf("Orig 0 ");
	  printVnlVector(timeSeries, 3);

	  timeSeries = timeSeries - timeSeries.mean();
	  printf("demean 1 ");
	  printVnlVector(timeSeries, 3);

	  // root mean square of the time series, is actually the sample
	  // std deviation when sample is zero mean.

     
	  // Normalize by deviding by the magnitude, which is not equal to
	  // the rms (sample std deviation). So this normalization is not
	  // the one in the statistics sense. But because after this
	  // normalization, the time series will be on sphere, we use it.
	  timeSeries.normalize();
	  printf("onsphere 2 ");
	  printVnlVector(timeSeries, 3);

	  // all the data points are at a hyperplane x_1 + x_2 + ... = 1
	  // because we normalize it by subtract the mean. we need to
	  // project such that its last elements (or the first, depending
	  // on the implenmentation) is zero, and we discard it, so the
	  // remaning vector is in a p-1 space, and p-2 sphere.

	  // v2 = Matrix * v1, where v1 is the coordinates under natural
	  // bases, and v2 is coordinates under hyper bases.
	  timeSeries.pre_multiply(transMat);
	  printf("trans 3 ");
	  printVnlVector(timeSeries, 3);

	  // Now the coordinates should have one element always zero. 

     }

     return 0;
}
