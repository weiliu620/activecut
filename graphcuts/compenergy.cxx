#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <new>
#include <iostream>
#include <iterator>
#include <string>
#include <fstream>

#include <itkNiftiImageIO.h>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
//#include <boost/program_options/cmdline.hpp>
//#include <boost/program_options/errors.hpp>
//#include <boost/program_options/option.hpp>
//#include <boost/program_options/value_semantic.hpp>
//#include <boost/program_options/version.hpp>

#include <Eigen/Core>
#include "GCoptimization.h"

using namespace std;
namespace po = boost::program_options;
// USING_PART_OF_NAMESPACE_EIGEN

int main(int argc, char* argv[])
{
     int _DBG = 0;
     float alpha = 0; 
     float beta = 0.7;
     float sigma[2] = {0.1, 0.1};
     float mu[2] = {-1, 1};
     int label = 0;
     float energy_ll = 0; // likelihood energy. i.e. the exponential.
     int num_labels = 2;
     int i = 0, j = 0;
     string labeledimage, observedimage;
     string config_file;

     typedef itk::Image<float, 2> ImageType2D;
     typedef itk::ImageFileReader< ImageType2D >  ReaderType;
     ReaderType::Pointer obs_reader = ReaderType::New();
     ReaderType::Pointer label_reader = ReaderType::New();

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help", "Given labeledimage, compute energy of this labelling.")
	  ("version,v", "print version string")
	  ("config,c", po::value<string>(&config_file)->default_value("graphcuts.cfg"),
	   "name of a file of a configuration.")
	  ("labeledimage,l", po::value<string>(), "labeled image file as input. ")
	  ("observed", po::value<string>(), "Noise image file as input. Data energy is between labeledimage and observed image.");

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

	  if (vm.count("labeledimage")) {
	       labeledimage = vm["labeledimage"].as<string>();
	  }
	  else {
//	       labeledimage = "labeledimage.nii";
	       cout << "need labeledimage as a parameter. Program terminated.\n";
	  }
	  if (vm.count("observed")) {
	       observedimage = vm["observed"].as<string>();
	  }
	  else {
	       observedimage = "observedimage.nii";

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

     obs_reader->SetFileName(observedimage);
     label_reader->SetFileName(labeledimage);
     obs_reader->Update();
     label_reader->Update();
     ImageType2D::Pointer label_pointer = label_reader->GetOutput();
     ImageType2D::Pointer obs_pointer = obs_reader->GetOutput();

     ImageType2D::IndexType pixelIndex;
     ImageType2D::SizeType size =
          label_reader->GetOutput()->GetLargestPossibleRegion().GetSize();


     GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(size[0], size[1], num_labels);
     printf(" size[0] = %d, size[1] = %d\n", size[0], size[1]);
     // data costs.
     float tmu = 0, tsigma = 0, tobs = 0;
     for (pixelIndex[0] = 0; pixelIndex[0] < size[0]; pixelIndex[0]++) {
     	  for (pixelIndex[1] = 0; pixelIndex[1] < size[1]; pixelIndex[1]++) {
	       tmu = mu[(int)label_pointer->GetPixel(pixelIndex)];
	       tsigma = sigma[(int)label_pointer->GetPixel(pixelIndex)];
	       tobs = obs_pointer->GetPixel(pixelIndex);
	       energy_ll = (tobs - tmu)*(tobs - tmu)/(2*tsigma*tsigma)
		    + log(tsigma);
	       
	       gc->setDataCost(pixelIndex[0]*size[1]+pixelIndex[1], label_pointer->GetPixel(pixelIndex), energy_ll);
		    
	  }
     }
	  
     // next set up smoothness costs individually
     gc->setSmoothCost(0, 0, -alpha - beta); 
     gc->setSmoothCost(0, 1, -alpha + beta); 
     gc->setSmoothCost(1, 0, -alpha + beta); 
     gc->setSmoothCost(1, 1, -alpha - beta); 

     // Assign labels according to labeled image.
     for (pixelIndex[0] = 0; pixelIndex[0] < size[0]; pixelIndex[0]++) {
     	  for (pixelIndex[1] = 0; pixelIndex[1] < size[1]; pixelIndex[1]++) {
	       // cout << "pixelIndex[0]*size[1]+pixelIndex[1] = " << pixelIndex[0]*size[1]+pixelIndex[1] << ", (int)label_pointer->GetPixel(pixelIndex) = " << (int)label_pointer->GetPixel(pixelIndex) << "\n";
	       gc->setLabel(pixelIndex[0]*size[1]+pixelIndex[1], (int)label_pointer->GetPixel(pixelIndex));
	  }
     }

     printf("labeled image. Total energy: %f, data enerty: %f, smooth energy: %f, label energy: %f.\n",gc->compute_energy(), gc->giveDataEnergy(), gc->giveSmoothEnergy(), gc->giveLabelEnergy());

}

