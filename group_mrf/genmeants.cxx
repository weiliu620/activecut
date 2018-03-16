#include "commonalt.h"
using std::cout;
using std::cin;
using std::string;
twister_base_gen_type mygenerator(42u);

int GenerateMu(vnl_vector<float> & mu, float phi, float arvar);

int main (int argc, char* argv[])
{
     unsigned verbose = 0;
     float phi = 0, arvar = 0, tsvar = 0;
     unsigned seed = 0;
     unsigned timeSeriesLength = 0;
     unsigned numClusters = 0;
     std::string observedimage, trueimage, meantsFile;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Command line options.");
     mydesc.add_options()
	  ("help,h", "Mean time series are from a AR process with phi and arvar as the correlation coefficient and noise variance.")
	  ("phi,p", po::value<float>(&phi)->default_value(0.5), "AR process coefficient.")	  
	  ("arvar,s", po::value<float>(&arvar)->default_value(0.5), "AR process white noise variance.")	  
	  ("timepoints,d", po::value<unsigned>(&timeSeriesLength)->default_value(300), "Time series length.")	
	  ("numClusters,k", po::value<unsigned>(&numClusters)->default_value(8),
	       "Number of labels. Default is 5.")
	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")  
	  ("out,o", po::value<std::string>(&meantsFile)->default_value("meants.txt"),
	   "ascii file that contain mean time series.")
	  ("verbose,v", po::value<unsigned int>(&verbose)->default_value(0), 
	   "Verbose level.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: generateimage [options]\n";
	       cout << mydesc << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  cout << e.what() << "\n";
	  return 1;
     }    

     // init random generator.
     mygenerator.seed(static_cast<unsigned int>(seed));


     // Allocate memory for mean time series.
     std::vector< vnl_vector<float> > allMu(numClusters);
     for (unsigned clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  allMu.at(clsIdx).set_size(timeSeriesLength);
	  allMu.at(clsIdx).fill(0);
	  GenerateMu(allMu[clsIdx], phi, arvar);
     }


     // write back.
     std::ofstream outstream;
     outstream.open(meantsFile.c_str());
     
     for (unsigned clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  for (unsigned tIdx = 0; tIdx < timeSeriesLength; tIdx ++) {
	       outstream << allMu[clsIdx][tIdx] << " ";
	  }
	  outstream << std::endl;
     }
     
     outstream.close();

}

int GenerateMu(vnl_vector<float> & mu, float phi, float arvar)
{
  
     // Uniform real random generator.
     boost::normal_distribution<double> normal_dist(0,sqrt(arvar)); // normal distribution
     boost::variate_generator<twister_base_gen_type&, boost::normal_distribution<> > nor(mygenerator, normal_dist);
     unsigned tsLength = mu.size();
     
     mu[0] = nor();

     for (unsigned tsIdx = 1; tsIdx < tsLength; tsIdx ++) {
	  mu[tsIdx] = phi * mu[tsIdx - 1] + nor();
     }

     
     
     
}
