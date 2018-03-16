#include <common.h>
#include <itkMaskedImageToHistogramFilter.h>

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     std::string mask_file, in_file;
     unsigned nbin = 100;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "histogram")
	  ("input,i", po::value<std::string>(&in_file)->default_value("image.nii.gz"), 
	   "input nifti file")
	  ("mask,m", po::value<std::string>(&mask_file)->default_value("mask.nii.gz"), 
	   "mask binary image.")
	  
	  ("nbin,n", po::value<unsigned>(&nbin)->default_value(100),
	   "number of bins.")

	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: rain [options]\n";
	       std::cout << mydesc << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    


     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(mask_file);
     ImageType3DChar::Pointer maskPtr = maskReader->GetOutput();
     maskPtr->Update();

     ReaderType3DFloat::Pointer imageReader = ReaderType3DFloat::New();
     imageReader->SetFileName(in_file);
     ImageType3DFloat::Pointer imagePtr = imageReader->GetOutput();
     imagePtr->Update();

     typedef itk::Statistics::MaskedImageToHistogramFilter< ImageType3DFloat, ImageType3DChar>   HistogramFilterType;
     typedef HistogramFilterType::HistogramMeasurementVectorType             HistogramMeasurementVectorType;
     typedef HistogramFilterType::HistogramSizeType                          HistogramSizeType;
     typedef HistogramFilterType::HistogramType                              HistogramType;

     HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();
     histogramFilter->SetInput(imagePtr);
     histogramFilter->SetMaskImage(maskPtr);
     histogramFilter->SetAutoMinimumMaximum(false);

     HistogramFilterType::HistogramMeasurementVectorType lowerBound( 1 );
     HistogramFilterType::HistogramMeasurementVectorType upperBound( 1 );
     lowerBound[0] = -1;
     upperBound[0] = 1;

     histogramFilter->SetHistogramBinMinimum( lowerBound );
     histogramFilter->SetHistogramBinMaximum( upperBound );
     
     const unsigned int MeasurementVectorSize = 1;
     HistogramSizeType histogramSize( MeasurementVectorSize );
     histogramSize[0] = nbin;  // number of bins for the Red   channel

     histogramFilter->SetHistogramSize(histogramSize);
     histogramFilter->SetMarginalScale(10); // Required (could this be set in the filter?)
     histogramFilter->Update();

     const HistogramType * histogram = histogramFilter->GetOutput();
     HistogramType::ConstIterator histogramIterator = histogram->Begin();

     while( histogramIterator  != histogram->End() )
     {
	  std::cout << "Index = " << histogramIterator.GetMeasurementVector() << " Frequency = " << histogramIterator.GetFrequency() << std::endl;
	  std::cout << "Index = " << histogramIterator.GetIndex() << "  Frequency = " << histogramIterator.GetFrequency() << std::endl;
	  ++histogramIterator ;
     }

}
