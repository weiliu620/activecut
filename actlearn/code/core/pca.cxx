#include <common.h>
#include <gmm.h>
#include <graphcuts.h>
#include <utility.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
     std::string input_image_file, out_file, mask_file;
     unsigned outdim = 0, verbose = 0;
     
     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Grab cut segmentation of 3D multi-channel traumatic brain injury (TBI) images.")

	  ("data,d", po::value<std::string>(&input_image_file),
	   "Input all-channel image file. A 4D gipl or nii or nii.gz file.")

	  ("mask,m", po::value<std::string>(&mask_file),
	   "binary mask file.")

	  ("outdim,p", po::value<unsigned>(&outdim)->default_value(3),
	   "output image dimension. must be smaller than input dim.")

	  ("out", po::value<std::string>(&out_file),
	   "output file.")

	  ("verbose,v", po::value<unsigned>(&verbose)->default_value(0),
	   "verbose level in [0, 3]. ")
	  ;

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: pca [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // read data images.
     ReaderType4DF::Pointer dataReader = ReaderType4DF::New();
     dataReader->SetFileName(input_image_file);
     dataReader->Update();
     ImageType4DF::Pointer dataPtr = dataReader->GetOutput();
     ImageType4DF::SizeType dataSize = dataPtr->GetLargestPossibleRegion().GetSize();
     unsigned n_channels = dataSize[3];
     ImageType4DF::IndexType dataIdx;

     // read mask images.
     ReaderType3DI::Pointer maskReader = ReaderType3DI::New();
     maskReader->SetFileName(mask_file);
     maskReader->Update();
     ImageType3DI::Pointer maskPtr = maskReader->GetOutput();
     IteratorType3DI maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());
     ImageType3DI::IndexType maskIdx;

     unsigned n_samples = 0;
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       n_samples ++;
	  } // mask > 0
     }

     if (verbose >= 1 )
	  printf("n_samples: %i\n", n_samples);

     vnl_matrix<double> data(n_samples, n_channels);
     vnl_matrix<double> out_data(n_samples, outdim);

     // move the data into 'data' matrix.
     unsigned n = 0; // sample id.
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
     	  if (maskIt.Get() > 0) {
     	       maskIdx = maskIt.GetIndex();
     	       dataIdx[0] = maskIdx[0];
     	       dataIdx[1] = maskIdx[1];
     	       dataIdx[2] = maskIdx[2];
     	       for (unsigned channel_id = 0; channel_id < n_channels; channel_id ++) {
     		    dataIdx[3] = channel_id;
     		    data(n, channel_id) = dataPtr->GetPixel(dataIdx);
     	       }
     	       n++;
     	  } // mask > 0
     }

     // compute mean of each channel.
     vnl_vector<double> mean(n_channels, 0);
     for (unsigned k = 0; k < n_channels; k++) {
     	  mean[k] = data.get_column(k).mean();
     }
     
     // substract mean
     for (unsigned k = 0; k < n_channels; k++) {
     	  data.set_column(k, data.get_column(k) - mean[k] );
     	  // assert (data.get_column(k).mean() == 0);
     }

     // compute covariance mat
     vnl_matrix<double> cov_mat = data.transpose() * data;
     vnl_symmetric_eigensystem <double> eig(cov_mat);

     out_data = data * eig.V.fliplr();

     
     // extract diagonal as vector, then flip so it's in decending order.
     vnl_vector<double> eigvalues = eig.D.diagonal();
     eigvalues.flip();
     printf("eigenvalues: ");
     print_vnl_vec(eigvalues, n_channels);
     
     vnl_vector<double> ratio = eigvalues / eigvalues.sum();
     printf("ratios of eigenvalues:  ");
     print_vnl_vec(ratio, n_channels);

     // add mean back.
     for (unsigned k = 0; k < n_channels; k++) {
     	  out_data.set_column(k, out_data.get_column(k) + mean[k] );
     }

     // create outimage.
     ImageType4DF::IndexType start;
     start.Fill(0);
     ImageType4DF::SizeType outimageSize = dataSize;
     outimageSize[3] = outdim;

     ImageType4DF::RegionType outimageRegion;
     outimageRegion.SetSize(outimageSize);
     outimageRegion.SetIndex(start);
     ImageType4DF::Pointer outimagePtr = ImageType4DF::New();
     outimagePtr->SetRegions(outimageRegion);
     outimagePtr->Allocate();
     outimagePtr->FillBuffer( 0 ); // init to zero.

     outimagePtr->SetOrigin( dataPtr->GetOrigin() );
     outimagePtr->SetSpacing(dataPtr->GetSpacing() );
     outimagePtr->SetDirection(dataPtr->GetDirection() );


     n = 0;
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
     	  if (maskIt.Get() > 0) {

     	       maskIdx = maskIt.GetIndex();
     	       dataIdx[0] = maskIdx[0];
     	       dataIdx[1] = maskIdx[1];
     	       dataIdx[2] = maskIdx[2];
     	       for (unsigned channel_id = 0; channel_id < outdim; channel_id ++) {
     		    dataIdx[3] = channel_id;
     		    outimagePtr->SetPixel(dataIdx, out_data(n, channel_id));
     	       }
     	       n++;
     	  } // mask > 0
     }


     WriterType4DF::Pointer writer = WriterType4DF::New();
	  
     writer->SetInput(outimagePtr);
     writer->SetFileName(out_file);
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

     std::cout << "save_gmm_labelmap(): File " << out_file << " saved.\n";

     return 0;
}
