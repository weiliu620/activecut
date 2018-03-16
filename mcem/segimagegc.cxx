#include "commonalt.h"
#include <string.h>
#include <time.h>


#include <Eigen/Core>
#include "GCoptimization.h"

using namespace std;
namespace po = boost::program_options;
USING_PART_OF_NAMESPACE_EIGEN

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(cout, " ")); 
    return os;
}

int main(int argc, char* argv[])
{
     int _DBG = 0;
     float alpha = 0; 
     float beta = 0.7;
     float sigma[2] = {2, 2};
     float mu[2] = {100, 110};
     int label = 0;
     float energy_ll = 0; // likelihood energy. i.e. the exponential.
     int num_labels = 2;
     int i = 0, j = 0;

     string observedimage, recoveredimage;
     string config_file;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help", "produce help message")
	  ("version,v", "print version string")
	  ("config,c", po::value<string>(&config_file)->default_value("graphcuts.cfg"),
	   "name of a file of a configuration.")
	  ("observed", po::value<string>(), "output noise image file")
	  ("recovered", po::value<string>(), "output recovered image file");

     // Declare a group of options that will be allowed both on
     // command line and in config file
     po::options_description config("Configuration options (can be given at both commandline and config file");
     config.add_options()
	  ("alpha", po::value<float>(&alpha)->default_value(0),
	   "Set MRF parameter alpha.")
	  ("beta", po::value<float>(&beta)->default_value(0.7), 
	   "Set MRF parameter beta.")
	  ("sigma", po::value<float>()->default_value(0.5), 
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
	  if (!ifs) {
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
	  
	  if (vm.count("recovered")) {
	       recoveredimage = vm["recovered"].as<string>();
	  }
	  else {
	       recoveredimage = "recoveredimage.nii";
	  }

	  if (vm.count("sigma")) {
	       sigma[0] = vm["sigma"].as<float>();
	       sigma[1] = vm["sigma"].as<float>();
	  }
     }
     catch(exception& e)
     {
	  cout << e.what() << "\n";
	  return 1;
     }    

     typedef itk::Image<float, 2> ImageType2D;
     typedef itk::ImageFileReader< ImageType2D >  ReaderType;
     typedef itk::ImageFileWriter< ImageType2D >  WriterType;

     ReaderType::Pointer reader = ReaderType::New();
     WriterType::Pointer writer = WriterType::New();

     reader->SetFileName(observedimage);
     writer->SetFileName(recoveredimage);
     
     reader->Update();
     ImageType2D::Pointer imP = reader->GetOutput();
     writer->SetInput(imP);

     ImageType2D::IndexType pixelIndex;
     ImageType2D::SizeType size =
          reader->GetOutput()->GetLargestPossibleRegion().GetSize();

     MatrixXd image;
     image = MatrixXd::Zero(size[0], size[1]);

     for (pixelIndex[0] = 0; pixelIndex[0] < size[0]; pixelIndex[0]++){
          for (pixelIndex[1] = 0; pixelIndex[1] < size[1]; pixelIndex[1]++){
	       image(pixelIndex[0], pixelIndex[1]) = imP->GetPixel(pixelIndex);
	  }
     }

     GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(image.cols(), image.rows(), num_labels);
     // data costs.
     for (i = 0; i < image.rows(); i++) {
     	  for (j = 0; j < image.cols(); j++) {
     	       for (label = 0; label < num_labels; label++) {
     		    energy_ll =  (image(i,j)-mu[label])*(image(i,j)-mu[label])/(2*sigma[label]*sigma[label]) + log(sigma[label]);
     		    gc->setDataCost(i*image.cols()+j, label, energy_ll);
		    
	       }
	  }
     }
	  
     // next set up smoothness costs individually
     gc->setSmoothCost(0, 0, -alpha - beta); 
     gc->setSmoothCost(0, 1, -alpha + beta); 
     gc->setSmoothCost(1, 0, -alpha + beta); 
     gc->setSmoothCost(1, 1, -alpha - beta); 

     printf("Before optimization. Total energy: %f, data enerty: %f, smooth energy: %f, label energy: %f.\n",gc->compute_energy(), gc->giveDataEnergy(), gc->giveSmoothEnergy(), gc->giveLabelEnergy());

     // run expansion for 2 iterations. 
     gc->expansion(1);
     printf("After expansion 1, Total energy: %f, data enerty: %f, smooth energy: %f, label energy: %f.\n ",gc->compute_energy(), gc->giveDataEnergy(), gc->giveSmoothEnergy(), gc->giveLabelEnergy());
     
     // Write back.
     for (pixelIndex[0] = 0; pixelIndex[0] < size[0]; pixelIndex[0]++){
          for (pixelIndex[1] = 0; pixelIndex[1] < size[1]; pixelIndex[1]++){
	       
	       image(pixelIndex[0], pixelIndex[1]) = imP->GetPixel(pixelIndex);
	       imP->SetPixel(pixelIndex, gc->whatLabel(pixelIndex[0]*image.cols()+pixelIndex[1]));
	  }
     }
     writer->Update();
}

