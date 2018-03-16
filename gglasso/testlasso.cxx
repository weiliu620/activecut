#include <common.h>
using std::cout;
using std::cin;

int codescent(Eigen::VectorXf & Y, // Nx1 response vector.
	      Eigen::MatrixXf & X, // NxP regressor matrix.
	      Eigen::MatrixXf & beta,
	      Codescent_parameters & cdparam,
	      std::string method,
	      unsigned verbose);

int main (int argc, char* argv[])
{

     std::string designmatrixfile, responsefile;
     float lambda = 0.5, lambda0 = 0.5, alpha = 0.5;
     unsigned int verbose = 0;
     unsigned int numObs = 0, numVar = 0, numSteps = 100;
     bool standardize = 1;
     std::string method;

     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Test lasso (and elastic net) regression.")
	  ("designmatrix,x", po::value<std::string>(&designmatrixfile)->default_value("designmatrix.txt"), "NxP regression matrix X.")
	  ("response,y", po::value<std::string>(&responsefile)->default_value("response.txt"), "Nx1 response variable Y.")
	  ("lambda,l", po::value<float>(&lambda)->default_value(0.5), "regularization parameter.")	  
	  ("initlambda,i", po::value<float>(&lambda0)->default_value(0.5), "Initial regularization parameter.")	  
	  ("numsweeps,s", po::value<unsigned>(&numSteps)->default_value(100), "Number of observations.")	  
	  ("method,m", po::value<std::string>(&method)->default_value("naive"), "optimization method. Can be 'naive', 'covupd'. ")	  
	  ("alpha,a", po::value<float>(&alpha)->default_value(0.5), "Tradeoff parameter between ridge (alpah=0) and lasso (alpha=1).")	  
	  ("numobs,N", po::value<unsigned>(&numObs)->default_value(100), "Number of observations.")	  
	  ("numpred,P", po::value<unsigned>(&numVar)->default_value(1000), "Number of predictors.")	  
	  ("standardize,t", po::value<bool>(&standardize)->default_value(1), "If standardize the X matrix such that each column has mean zero and variance one. Binary variable either 0 or 1.")	  
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

     Eigen::MatrixXf predictor(numObs, numVar);
     // read predictor samples from design matrix file.
     std::string varValue;
     std::ifstream predStream(designmatrixfile.c_str() );
     if (predStream.is_open() ) {
	  for (unsigned obsIdx = 0; obsIdx < numObs; obsIdx ++) {
	       for (unsigned varIdx = 0; varIdx < numVar; varIdx ++) {
		    predStream >> varValue;
		    predictor(obsIdx, varIdx) = boost::lexical_cast<float> (varValue);

	       }
	  }
     }
     else {
	  std::cout << "can not open file "  << designmatrixfile << std::endl;
     }

     // std::cout << "predictor:\n" << predictor << std::endl;

     // read response variables.
     Eigen::VectorXf response(numObs);
     std::ifstream responseStream(responsefile.c_str() );
     if (responseStream.is_open() ) {
	  for (unsigned obsIdx = 0; obsIdx < numObs; obsIdx ++) {
	       responseStream >> varValue;
	       response(obsIdx) = boost::lexical_cast<float> (varValue);
	  }
     }
     else {
	  std::cout << "can not open file "  << responsefile << std::endl;
     }

     // std::cout << "response:\n" << response << std::endl;

     // normalize such that mean is zero and variance is 1.
     if (standardize) {
	  for (unsigned varIdx = 0; varIdx < numVar; varIdx ++) {

	       // convert to array to do element wise sum. It's weird that "/" below
	       // does not need that.
	       predictor.col(varIdx) = predictor.col(varIdx).array() - predictor.col(varIdx).mean();
	       predictor.col(varIdx) = predictor.col(varIdx) / ((1/sqrt((float)numObs)) * (predictor.col(varIdx).norm()));
	  }
     }

     Eigen::MatrixXf beta(numVar, numSteps);

     Codescent_parameters cdparam;
     cdparam.finalLambda = lambda;
     cdparam.initLambda = lambda0;
     cdparam.alpha = alpha;
     cdparam.numSteps = numSteps;
     cdparam.stopThr = 0.00001;
     cdparam.numObs = numObs;
     cdparam.numVar = numVar;
     cdparam.eps = 1e-6;

     codescent(response,
	       predictor,
	       beta,
	       cdparam,
	       method,
	       verbose);

     // Write back.
     std::ofstream outstream;
     outstream.open("beta.txt");
     for (unsigned sweepIdx = 0; sweepIdx < numSteps; sweepIdx ++) {     
	  for (unsigned varIdx = 0; varIdx < numVar; varIdx ++) {
	       outstream << beta(varIdx, sweepIdx) << " ";
	  }
	  outstream << std::endl;
     }
     
     outstream.close();


}
