#include "commonalt.h"
using std::cout;
using std::string;

int main (int argc, char* argv[])
{
     unsigned verbose = 0;
     float kappa;
     unsigned numSamples = 0;
     unsigned timeSeriesLength = 0;
     unsigned numClusters = 0;

     std::string observedimage, trueimage, maskimage;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "produce help message")
	  ("kappa,k", po::value<float>(&kappa)->default_value(50), "concentration parameter of von Mise. i.e. kappa.")	  
	  ("numsamples,n", po::value<unsigned>(&numSamples)->default_value(300), "Number of samples. ")	  
	  ("verbose,v", po::value<unsigned int>(&verbose)->default_value(0), 
	   "Verbose level.");

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
	       cout << cmdline_options << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  cout << e.what() << "\n";
	  return 1;
     }    


     vnl_vector<float> vMFSample, newSample;
     vMFSample.set_size(3);
     vnl_vector<float> mu, northPole;
     mu.set_size(3);

     northPole.set_size(3);
     northPole[0] = 0;
     northPole[1] = 0;
     northPole[2] = 1;
     vnl_matrix<float> rotationMat;
     rotationMat.set_size(3,3);


     
     for (int sampleIdx = 0; sampleIdx < numSamples; sampleIdx++) {
	  generateVonMiseWood(vMFSample, kappa);
	  printf("0 %f %f %f\n", vMFSample[0], vMFSample[1], vMFSample[2]);

	  // rotation.
	  mu[0] = 0;
	  mu[1] = 1;
	  mu[2] = 0;
	  getRotationMat(rotationMat, northPole, mu);
	  newSample = rotationMat * vMFSample;
	  printf("1 %f %f %f\n", newSample[0], newSample[1], newSample[2]);

	  // rotation.
	  mu[0] = 1;
	  mu[1] = 0;
	  mu[2] = 0;
	  getRotationMat(rotationMat, northPole, mu);
	  newSample = rotationMat * vMFSample;
	  printf("2 %f %f %f\n", newSample[0], newSample[1], newSample[2]);


	  // rotation.
	  mu[0] = 0;
	  mu[1] = -1;
	  mu[2] = 0;
	  getRotationMat(rotationMat, northPole, mu);
	  newSample = rotationMat * vMFSample;
	  printf("3 %f %f %f\n", newSample[0], newSample[1], newSample[2]);


	  // rotation.
	  mu[0] = -1;
	  mu[1] = 0;
	  mu[2] = 0;
	  getRotationMat(rotationMat, northPole, mu);
	  newSample = rotationMat * vMFSample;
	  printf("4 %f %f %f\n", newSample[0], newSample[1], newSample[2]);

     }
     
}
