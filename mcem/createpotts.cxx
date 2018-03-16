#include "common.h"

// Init random generator as a global object.
twister_base_gen_type mygenerator(42u);

int main (int argc, char* argv[])
{
     int _DBG = 0;
     int N = 200; // image size NxN.
     float alpha = 0; 
     float beta = 0.7;
     float betaUp = 1, betaDown = 1, betaLeft = 1, betaRight = 1;
     float sigma = 1;
     float scale = 10;
     float offset = 100;

     int i = 0, j = 0;
     int scan = 0;
     int num_scan;
     int nlabels;
     float denergy = 0; // energy change = candidate energy - current energy.
     int cand;
     double p_acpt = 0; // probability to accept the candidate sample.
     string observedimage, trueimage;
     string config_file;
     unsigned seed = 0;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "produce help message")

	  ("imagesize,z", po::value<int>(&N)->default_value(128),
	   "Image size. Image is square.")

	  ("beta", po::value<float>(&beta)->default_value(2), 
	   "Set MRF parameter beta.")
	  ("betaup", po::value<float>(&betaUp)->default_value(2), 
	   "Set MRF parameter betaUp")
	  ("betadown", po::value<float>(&betaDown)->default_value(2), 
	   "Set MRF parameter beta.")
	  ("betaleft", po::value<float>(&betaLeft)->default_value(2), 
	   "Set MRF parameter betaLeft.")
	  ("betaright", po::value<float>(&betaRight)->default_value(2), 
	   "Set MRF parameter betaRight.")

	  ("sigma", po::value<float>(&sigma)->default_value(2), 
	   "Set noise stand deviation sigma. Default is 2.")

	  ("nlabels", po::value<int>(&nlabels)->default_value(4),
	       "Number of labels. Default is 5.")

	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	       "random generator seed.")

	  ("scan,n", po::value<int>(&num_scan)->default_value(500),
	   "Number of scan on MRF.")

	  ("observed,o", po::value<string>(), "output noise image file")
	  ("true,t", po::value<string>(), "output true image file")

	  ("alpha", po::value<float>(&alpha)->default_value(0),
	   "Set MRF parameter alpha.")

	  ("scale", po::value<float>(&scale)->default_value(10),
	   "likelihood function mu = scale * label + offset. Default 10.")
	  ("offset", po::value<float>(&offset)->default_value(100),
	   "likelihood function mu = scale * label + offset. Default 100.")

	  ("debuglevel", po::value<int>(&_DBG)->default_value(0), 
	   "set debug level. 0 for no debug output. 3 for most debug output.");

     // Declare a group of options that will be allowed only in config
     // file. 
     
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
	       cout << "Usage: generateimage [options]\n";
	       cout << cmdline_options << "\n";
	       return 0;
	  }

	  if (vm.count("observed")) {
	       observedimage = vm["observed"].as<string>();
	  }
	  else {
	       observedimage = "observedimage.nii";
	  }
	  
	  if (vm.count("true")) {
	       trueimage = vm["true"].as<string>();
	  }
	  else {
	       trueimage = "trueimage.nii";
	  }
     }
     catch(exception& e)
     {
	  cout << e.what() << "\n";
	  return 1;
     }    
     
     MatrixXf image(N, N);
//     mygenerator.seed(static_cast<unsigned int>(0));
     mygenerator.seed(static_cast<unsigned int>(seed));

     // Uniform integer generator.
     boost::uniform_int<> uni_int(1, nlabels); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);
     // normal distribution generator.
     boost::normal_distribution<> normal_dist(0, 1); // Normal distribution.
     boost::variate_generator<twister_base_gen_type&, boost::normal_distribution<> > nor(mygenerator, normal_dist);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     // Init image.
     for (i = 0; i < image.rows(); i++) {
	  for (j = 0; j < image.cols(); j++) {
	       image(i,j) = roll_die() - 1;
	  }
     }     

     SaveToy(image, "randomimage.nii");
     
     // sampling Markov Random Field.
     while (scan < num_scan) {
	  scan ++;
	  for (i = 0; i < image.rows(); i++) {
	       for (j = 0; j < image.cols(); j++) {
		    denergy = 0;
		    cand = roll_die() - 1;
		    
		    if (i > 0) {denergy = denergy + 
			      betaLeft * (int(cand != image(i-1,j)) - int(image(i,j) != image(i-1,j)) );}
		    if (i < image.rows()-1) {denergy = denergy +
			      betaRight * (int(cand != image(i+1,j)) - int(image(i,j) != image(i+1,j)) );}
		    if (j > 0) {denergy = denergy + 
			      betaDown * ( int(cand != image(i,j-1)) - int(image(i,j) != image(i,j-1)) );}
		    if (j < image.cols()-1) {denergy = denergy + 
			      betaUp * ( int(cand != image(i,j+1)) - int(image(i,j) != image(i,j+1)) );}
//		    denergy = beta * denergy;
		    // if energy change less than zero, just accept
		    // candidate. otherwise accept with exp(- energy
		    // change).
		    if (denergy <= 0) {
			 image(i,j) = cand;
		    }
		    else {
			 p_acpt = exp(-denergy);
			 if (uni() < p_acpt) {
			 image(i,j) = cand;
			 }
		    }
	       }
	  }
	  if (scan%50 == 0 ) {
	       printf("scan %d. \n", scan);
	       SaveToy(image, "trueimage.nii");
	  }
     }
     SaveToy(image, "trueimage.nii");

     // Emit observation variables.
     for (i = 0; i < image.rows(); i ++) {
	  for (j = 0; j < image.cols(); j++) {
	       image(i,j) = image(i,j) * scale + offset + nor() * sigma;
	       // in case sample is less than zero. Truncate it. This
	       // should be very rare, and have no effect on model
	       // parameter estimation.
	       image(i,j) = image(i,j)>0?image(i,j):0;
	  }
     }

     if (_DBG > 1) {
	  cout << "After adding noise:\n" << image << "\n";     
     }
     SaveToy(image, observedimage);

}

