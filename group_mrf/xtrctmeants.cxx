#include "commonalt.h"
#include "itkCastImageFilter.h"
twister_base_gen_type mygenerator(42u);
#include <fstream>

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;

     int numClusters = 5;
     std::string fmriFile, meantsFile, maskFile, trueFile, labelFile;

     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "Given fmri file and label file, extract mean time series of each component and save it to a file")
	  ("numClusters,k", po::value<int>(&numClusters)->default_value(4),
	       "Number of labels. Default is 5.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")
	  ("fmri,f", po::value<std::string>(&fmriFile)->default_value("observedimage.nii"),
	   "fmri image file")
	  ("labelmap, l", po::value<std::string>(&labelFile)->default_value("labelmap.nii"),
	   "label map image.")
	  ("out,o", po::value<std::string>(&meantsFile)->default_value("meants.txt"),
	   "ascii file that contain mean time series.");

     po::options_description cmdline_options;
     cmdline_options.add(generic);

     po::positional_options_description p;
     p.add("in", -1);

     po::variables_map vm;        
     po::store(po::command_line_parser(argc, argv).
	       options(cmdline_options).positional(p).run(), vm);
     po::notify(vm);    

     try {

	  if (vm.count("help") | argc == 1) {
	       std::cout << "Usage: segimagealt [options]\n";
	       std::cout << cmdline_options << "\n";
	       return 0;
	  }
     }

     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // read in fmri file.
     ReaderType4DFloat::Pointer fmriReader = ReaderType4DFloat::New();
     fmriReader->SetFileName(fmriFile);
     fmriReader->Update();
     ImageType4DFloat::Pointer fmriPtr = fmriReader->GetOutput();
     ImageType4DFloat::SizeType fmriSize = fmriPtr->GetLargestPossibleRegion().GetSize();
     ImageType4DFloat::IndexType fmriIdx;

     // read label file.
     ReaderType3DChar::Pointer labelReader = ReaderType3DChar::New();
     labelReader->SetFileName(labelFile);
     ImageType3DChar::Pointer labelPtr = labelReader->GetOutput();
     labelPtr->Update();
     IteratorType3DChar labelIt(labelPtr, labelPtr->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType labelIdx;


     // computer mean time series for all components.
     // init cls
     vnl_vector <float> timeSeries(fmriSize[3], 0);
     std::vector < CompType >  cls;
     cls.resize(numClusters);

     unsigned clsIdx = 0;
     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  cls[clsIdx].mu.set_size(fmriSize[3]);
     }

     // reset all mu to zero.
     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  cls[clsIdx].numPoints = 0;
	  cls[clsIdx].mu = 0;
     }

     for (labelIt.GoToBegin(); !labelIt.IsAtEnd(); ++ labelIt) {
	  if (labelIt.Get() > 0) {
	       labelIdx = labelIt.GetIndex();
	       fmriIdx[0] = labelIdx[0];
	       fmriIdx[1] = labelIdx[1];
	       fmriIdx[2] = labelIdx[2];
	  
	       // First save time series at this position into a vector.
	       for (fmriIdx[3] = 0; fmriIdx[3] < fmriSize[3]; fmriIdx[3] ++) {
		    timeSeries[fmriIdx[3]] = fmriPtr->GetPixel(fmriIdx);
	       }

	       clsIdx = labelIt.Get() - 1;
	       cls[clsIdx].mu += timeSeries;
	       cls[clsIdx].numPoints ++;
	  }
     }

     // compute mu_k and meanNorm.

     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  // meanNorm see page 1350 of "Clustering on the Unit
	  // Hypersphere using von Mises-Fisher distributions" by
	  // Banerjee.
	  cls[clsIdx].meanNorm = cls[clsIdx].mu.two_norm() / cls[clsIdx].numPoints;
	  cls[clsIdx].mu.normalize();
     }



     // write back.
     std::ofstream outstream;
     outstream.open(meantsFile.c_str());
     
     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  for (unsigned tIdx = 0; tIdx < fmriSize[3]; tIdx ++) {
	       outstream << cls[clsIdx].mu[tIdx] << " ";
	  }
	  outstream << std::endl;
     }
     
     outstream.close();

}
