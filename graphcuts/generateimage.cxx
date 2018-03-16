#include <common.h>
using namespace std;
int SaveToy(MatrixXf image, string fname);
int label2hidden(MatrixXf & image);
int hidden2label(MatrixXf & image);

namespace po = boost::program_options;

int main (int argc, char* argv[])
{
     int _DBG = 0;
     int N = 200; // image size NxN.
     float alpha = 0; 
     float beta = 0.7;
     float sigma = 1;

     int i = 0, j = 0;
     int scan = 0;
     int num_scan;
     float energy = 0;
     float p_keep = 0; // probability to keep current sampled value.
     string observedimage, trueimage;
     string config_file;

     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "produce help message")
	  ("version,v", "print version string")
	  ("config,c", po::value<string>(&config_file)->default_value("graphcuts.cfg"),
	   "name of a file of a configuration.")
	  ("imagesize,z", po::value<int>(&N)->default_value(200),
	   "Image size. Image is square.")
	  ("scan,n", po::value<int>(&num_scan)->default_value(500),
	   "Number of scan on MRF.")
	  ("observed,o", po::value<string>(), "output noise image file")
	  ("true,t", po::value<string>(), "output true image file");

     // Declare a group of options that will be allowed both on
     // command line and in config file
     po::options_description config("Configuration options (can be given at both commandline and config file");
     config.add_options()

	  ("alpha,a", po::value<float>(&alpha)->default_value(0),
	   "Set MRF parameter alpha.")
	  ("beta,b", po::value<float>(&beta)->default_value(0.7), 
	   "Set MRF parameter beta.")
	  ("sigma,s", po::value<float>()->default_value(0.5), 
	   "Set noise stand deviation sigma.")
	  ("debuglevel", po::value<int>(&_DBG)->default_value(0), 
	   "set debug level. 0 for no debug output. 3 for most debug output.");
     
     po::options_description cmdline_options;
     cmdline_options.add(generic).add(config);

     po::options_description config_file_options;
     config_file_options.add(config);

     po::positional_options_description p;
     p.add("in", -1);

     po::variables_map vm;        
     po::store(po::command_line_parser(argc, argv).
	       options(cmdline_options).positional(p).run(), vm);
     po::notify(vm);    

     try {
	  ifstream ifs(config_file.c_str());
	  if (!ifs)
	  {
	       cout << "can not open config file: " << config_file << "\n";
	       return 0;
	  }
	  else
	  {
	       store(parse_config_file(ifs, config_file_options), vm);
	       notify(vm);
	  }

	  if (vm.count("help")) {
	       cout << "Usage: generateimage [options]\n";
	       cout << cmdline_options << "\n";
	       return 0;
	  }
     
	  if (vm.count("version")) {
	       cout << " version 1.0.\n";
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

	  if (vm.count("sigma")) {
	       sigma = vm["sigma"].as<float>();

	  }
     }
     catch(exception& e)
     {
	  cout << e.what() << "\n";
	  return 1;
     }    
     
     // Init random generator.
     twister_base_gen_type mygenerator(42u);
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::normal_distribution<> normal_dist(0, 1); // Normal distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);
     boost::variate_generator<twister_base_gen_type&, boost::normal_distribution<> > nor(mygenerator, normal_dist);
     mygenerator.seed(static_cast<unsigned int>(std::time(0)));


     MatrixXf image(N, N);

     // Init image.

     for (i = 0; i < image.rows(); i++) {
	  for (j = 0; j < image.cols(); j++) {
	       image(i,j) = uni()>0.5?1:-1;
	  }
     }     

     // sampling Markov Random Field.
     while (scan < num_scan) {
	  scan ++;
	  for (i = 0; i < image.rows(); i++) {
	       for (j = 0; j < image.cols(); j++) {
		    energy = 0;
		    energy = energy - alpha * image(i,j);
		    if (i > 0) {energy = energy - beta * image(i,j) * image(i-1,j);}
		    if (i < image.rows()-1) {energy = energy - beta * image(i,j) * image(i+1,j);}
		    if (j > 0) {energy = energy - beta * image(i,j) * image(i,j-1);}
		    if (j < image.cols()-1) {energy = energy - beta * image(i,j) * image(i,j+1);}
		    p_keep = exp(-energy)/(2*cosh(-energy));
		    if (uni() > p_keep) {
			 // flip the value.
			 image(i,j) = (-1) * image(i,j);
		    }
	       }
	  }
	  if (scan%10 == 0 && _DBG > 0) {
	       printf("scan %d. Savd image.\n", scan);
	       SaveToy(image, "testimage");
	  }
     }

     // Convert hidden variable map to label map.
     hidden2label(image);
     SaveToy(image, trueimage);
     if (_DBG > 1) {
	  cout << "True label image: \n" << image << "\n";
     }
     // Convert label map to hidden variable map.
     label2hidden(image);
     
     // Emit observation variables.
     for (i = 0; i < image.rows(); i ++) {
	  for (j = 0; j < image.cols(); j++) {
	       image(i,j) = image(i,j) + nor() * sigma;
	  }
     }

     if (_DBG > 1) {
	  cout << "After adding noise:\n" << image << "\n";     
     }
     SaveToy(image, observedimage);

}

int SaveToy(MatrixXf image, string fname)

{
     typedef itk::Image<float, 2> ImageType2D;
     typedef itk::ImageFileWriter< ImageType2D >  WriterType;
     WriterType::Pointer writer = WriterType::New();

     ImageType2D::Pointer imP = ImageType2D::New();
     ImageType2D::IndexType start;               
     
     start[0] =   0;  // first index on X
     start[1] =   0;  // first index on Y

     ImageType2D::SizeType  size;
     size[0]  = image.rows();  // size along X
     size[1]  = image.cols();  // size along Y
     printf("in Savetoy, image.rows = %d, image.cols = %d\n", image.rows(), image.cols());

     ImageType2D::RegionType region;
     region.SetSize( size );
     region.SetIndex( start );

     imP->SetRegions( region );
     imP->Allocate();

     // The image buffer is initialized to a particular value
     imP->FillBuffer(0);

     ImageType2D::IndexType pixelIndex;

     for (pixelIndex[0] = 0; pixelIndex[0] < image.rows(); pixelIndex[0]++){
          for (pixelIndex[1] = 0; pixelIndex[1] < image.cols(); pixelIndex[1]++){
	       imP->SetPixel(pixelIndex, image(pixelIndex[0], pixelIndex[1]));
	  }
     }

     writer->SetInput(imP);
     writer->SetFileName(fname);
     try 
     { 
     	  writer->Update(); 
     } 
     catch( itk::ExceptionObject & err ) 
     { 
     	  std::cerr << "ExceptionObject caught !" << std::endl; 
     	  std::cerr << err << std::endl; 
     	  return EXIT_FAILURE;
     } 
     return 0;
}

