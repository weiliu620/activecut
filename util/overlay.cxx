#include <common.h>
#include <itkLabelOverlayImageFilter.h>
#include <itkRGBPixel.h>

int main(int argc, char* argv[])
{
     std::string gray_file, label_file, out_file;
     unsigned short dimen = 0, slice_id = 0, verbose = 0;
     double opacity = 0.5;

     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "overlay label image on to gray scale image and draw color map. ")
	  ("input,i", po::value<std::string>(&gray_file)->default_value("gray.nii.gz"), 
	   "input gray scale nifti file")
	  ("label,l", po::value<std::string>(&label_file)->default_value("label.nii.gz"), 
	   "label image.")
 
	  ("out,o", po::value<std::string>(&out_file)->default_value("out.nii.gz"), 
	   "output 2d overlaied map.")

	  ("dimension,d", po::value<unsigned short>(&dimen)->default_value(2), 
	   "which dimensino to extract slice, 0 for x, 1 for y, and 2 for z.")

	  ("slice,s", po::value<unsigned short>(&slice_id)->default_value(30), 
	   "the id of the slice to be extracted. Zero based value")

	  ("opacity,p", po::value<double>(&opacity)->default_value(0.5), 
	   "opacity of the label map. 0 for all transparent, and 1 for all opaque (only show labels).")

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


     // read in gray intensity and label file. 
     ReaderType3DS::Pointer labelReader = ReaderType3DS::New();
     labelReader->SetFileName(label_file);
     ImageType3DS::Pointer labelPtr = labelReader->GetOutput();
     labelPtr->Update();

     ReaderType3DFloat::Pointer grayReader = ReaderType3DFloat::New();
     grayReader->SetFileName(gray_file);
     ImageType3DFloat::Pointer grayPtr = grayReader->GetOutput();
     grayPtr->Update();

     // extract slice from label image
     ImageType3DS::IndexType label_idx;
     label_idx.Fill(0);
     label_idx[dimen] = slice_id;
     ImageType3DS::SizeType label_size = labelPtr->GetLargestPossibleRegion().GetSize();
     label_size[dimen] = 0;
     typedef itk::ExtractImageFilter<ImageType3DS, ImageType2DS> ExtFilter;
     ExtFilter::Pointer ext_filter = ExtFilter::New();

     ImageType3DS::RegionType label_region(label_idx, label_size);
     ext_filter->SetExtractionRegion(label_region);
     ext_filter->SetInput(labelPtr);
     ext_filter->Update();
     ImageType2DS::Pointer labelSlicePtr = ext_filter->GetOutput();

     // extract slice from label image
     ImageType3DF::IndexType gray_idx;
     gray_idx.Fill(0);
     gray_idx[dimen] = slice_id;
     ImageType3DF::SizeType gray_size = grayPtr->GetLargestPossibleRegion().GetSize();
     gray_size[dimen] = 0;
     typedef itk::ExtractImageFilter<ImageType3DF, ImageType2DF> ExtFilter_gray;
     ExtFilter_gray::Pointer ext_filter_gray = ExtFilter_gray::New();

     ImageType3DF::RegionType gray_region(gray_idx, gray_size);
     ext_filter_gray->SetExtractionRegion(gray_region);
     ext_filter_gray->SetInput(grayPtr);
     ext_filter_gray->Update();
     ImageType2DF::Pointer graySlicePtr = ext_filter_gray->GetOutput();


     // do the overlay
     typedef itk::RGBPixel<unsigned char> RGBPixelType;
     typedef itk::Image<RGBPixelType> RGBImageType;
 
     typedef itk::LabelOverlayImageFilter<ImageType2DF, ImageType2DS, RGBImageType> 
     	  LabelOverlayImageFilterType;
     LabelOverlayImageFilterType::Pointer labelOverlayImageFilter = LabelOverlayImageFilterType::New();
     labelOverlayImageFilter->SetInput(graySlicePtr);
     labelOverlayImageFilter->SetLabelImage(labelSlicePtr);
     labelOverlayImageFilter->SetOpacity(opacity);
     labelOverlayImageFilter->Update();

     typedef  itk::ImageFileWriter< RGBImageType  > WriterType;
     WriterType::Pointer writer = WriterType::New();

     writer->SetFileName(out_file);
     writer->SetInput(labelOverlayImageFilter->GetOutput());
     writer->Update();
}
