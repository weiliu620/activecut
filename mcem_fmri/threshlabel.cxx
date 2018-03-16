#include "commonalt.h"
#include "itkCastImageFilter.h"
int ReadInitLabel(std::string initLabelFile, ImageType4DChar::Pointer samplePtr);

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;

     int numClusters = 5;
     unsigned short targetLabel = 0;
     float thresh = 0.3;
     bool doAllLabels = 0;
     std::string fmriFile, newLabelFile, maskFile, trueFile, labelFile;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "threshold on the dot product between fmri time series and its mean time series, for single label, or for all labels.")
	  ("numClusters,k", po::value<int>(&numClusters)->default_value(4),
	       "Number of labels. Default is 5.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")
	  ("fmri,f", po::value<std::string>(&fmriFile)->default_value("observedimage.nii"),
	   "fmri image file")
	  ("labelmap, l", po::value<std::string>(&labelFile)->default_value("labelmap.nii"),
	   "label map image.")
	  ("label, c", po::value<unsigned short>(&targetLabel),
	   "label map image. 1 based index.")
	  ("thresh, t", po::value<float>(&thresh)->default_value(0.3),
	   "maximal dot product between a time series and its mean")
	  ("out,o", po::value<std::string>(&newLabelFile)->default_value("newlabelfile.nii"),
	   "output labeled image file");

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
	  
	  if (vm.count("label")) {
	       doAllLabels = 0;
	  }
	  else {
	       doAllLabels = 1;
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



     // compute mean of the sample correlation.
     vnl_vector <float> mean(numClusters, 0);
     vnl_vector <float> numPoints(numClusters, 0);
     float thisDotProd = 0;
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
	       thisDotProd = inner_product(timeSeries, cls[clsIdx].mu);
	       
	       mean[clsIdx] += thisDotProd;
	       numPoints[clsIdx] ++;
	  } // labelIt.Get() > 0
	       
     }
     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  mean[clsIdx] = mean[clsIdx] / numPoints[clsIdx];
     }

     // compute std of the sample correlation.
     vnl_vector <float> std(numClusters, 0);
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
	       thisDotProd = inner_product(timeSeries, cls[clsIdx].mu);
	       std[clsIdx] += (thisDotProd - mean[clsIdx]) * (thisDotProd - mean[clsIdx]);
	       // if (clsIdx == 1) {
	       // 	    printf("idx[%i][%i][%i, ]thisDotProd = %f, std = %f\n", fmriIdx[0], fmriIdx[1], fmriIdx[2], thisDotProd, std[clsIdx]);
	       // }
	  } // labelIt.Get() > 0
     }

     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  std[clsIdx] = sqrt (std[clsIdx] / (numPoints[clsIdx] - 1) );
     }


     for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	  printf("cls[%i], mean = %f, std = %f.\n", clsIdx, mean[clsIdx], std[clsIdx]);
     }     

     // Threshold label map.
     bool keep = 0;
     unsigned numKeep = 0, numEra = 0;
     for (labelIt.GoToBegin(); !labelIt.IsAtEnd(); ++ labelIt) {
	  if (labelIt.Get() > 0) {
	       labelIdx = labelIt.GetIndex();
	       fmriIdx[0] = labelIdx[0];
	       fmriIdx[1] = labelIdx[1];
	       fmriIdx[2] = labelIdx[2];
	  
	  
	       if (doAllLabels |  labelIt.Get() == targetLabel ) {

		    // First save time series at this position into a vector.
		    for (fmriIdx[3] = 0; fmriIdx[3] < fmriSize[3]; fmriIdx[3] ++) {
			 timeSeries[fmriIdx[3]] = fmriPtr->GetPixel(fmriIdx);
		    }

		    clsIdx = labelIt.Get() - 1;
		    thisDotProd = inner_product(timeSeries, cls[clsIdx].mu);
		    keep = (thisDotProd - mean[clsIdx])/std[clsIdx] > thresh;

		    if (!keep) { labelIt.Set(0); numEra ++;}
		    else { numKeep ++; }
	       } 
	       else {
		    labelIt.Set(0); 
	       }		    
	  } // labelIt.Get() > 0
     }
     
     printf("numKeep = %i, numEra = %i\n", numKeep, numEra);

     // write back.

     // fslview bug...can not recognize char, have to use float.
     // WriterType3DChar::Pointer writer = WriterType3DChar::New();
     typedef itk::CastImageFilter<ImageType3DChar, ImageType3DFloat> CastFilterType;
     CastFilterType::Pointer castFilter = CastFilterType::New();
     castFilter->SetInput(labelPtr);
     
     WriterType3DFloat::Pointer writer = WriterType3DFloat::New();
     writer->SetInput(castFilter->GetOutput() );
     writer->SetFileName(newLabelFile);
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

     std::cout << "File " << newLabelFile << " saved.\n";



}
