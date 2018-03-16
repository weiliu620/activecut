#include <common.h>
#include <itkRGBPixel.h>
#include "itkScalarToRGBColormapImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

typedef itk::ExtractImageFilter<ImageType3DF, ImageType2DF> ExtFilter;
typedef itk::ExtractImageFilter<ImageType3DF, ImageType2DF> ExtFilter_gray;
typedef itk::RGBPixel<unsigned char> RGBPixelType;
typedef itk::Image<RGBPixelType> RGBImageType;
typedef  itk::ImageFileWriter< RGBImageType  > RGBWriterType;
typedef itk::ImageRegionIterator< RGBImageType>       IteratorTypeRGB;

struct Triplet {
unsigned r;
unsigned g;
unsigned b; 
};

int main(int argc, char* argv[])
{
     std::string gray_file, stat_file, out_file, colormap_file;
     unsigned short dimen = 0, slice_id = 0, verbose = 0;
     double opacity = 0.5, stat_min = 0, stat_max = 1;
     bool bgwhite = false, clip_high = false, clip_low = false;

     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "overlay statistics value onto gray scale image and draw a heat map. ")
	  ("background,b", po::value<std::string>(&gray_file)->default_value("gray.nii.gz"), 
	   "Background structural image (T1)")
	  ("stat,t", po::value<std::string>(&stat_file)->default_value("stat.nii.gz"), 
	   "statistical image. Need to be the same size with background T1 image.")
 
	  ("out,o", po::value<std::string>(&out_file)->default_value("out.png"), 
	   "output 2d overlaid image (usually png)")

	  ("colormap,c", po::value<std::string>(&colormap_file)->default_value("colormap.txt"), 
	   "Plain text file of Nx3, for the color map.")

	  ("dimension,d", po::value<unsigned short>(&dimen)->default_value(2), 
	   "which dimensino to extract slice, 0 for x, 1 for y, and 2 for z.")

	  ("slice,s", po::value<unsigned short>(&slice_id)->default_value(30), 
	   "the id of the slice to be extracted. Zero based value")

	  ("bgwhite,w", po::bool_switch(&bgwhite), 
	   "if use white background. Default is no and use black.")

	  ("cliphigh,h", po::bool_switch(&clip_high), 
	   "set true to make stat above statmax map to the color of statmax. Set false to not show it. Default is false. ")

	  ("cliplow,l", po::bool_switch(&clip_low), 
	   "set true to make stat below statmin map to the color of statmin. Set false to not show it. Default is false. ")

	  ("statmin,m", po::value<double>(&stat_min)->default_value(0), 
	   "the min value of stat map to show. Values below this value not shown. ")

	  ("statmax,x", po::value<double>(&stat_max)->default_value(1), 
	   "the max value of stat map to show. Values ablove this value not shown. ")

	  ("opacity,p", po::value<double>(&opacity)->default_value(0.5), 
	   "opacity of the stat map. 0 for all transparent, and 1 for all opaque (only show stats).")

	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: overlayheat [options]\n";
	       std::cout << mydesc << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // read colormap. 
     std::string rgb_string;
     std::ifstream instream(colormap_file.c_str() );

     Triplet triplet;
     std::vector<Triplet> colormap;
     while(std::getline(instream, rgb_string)) {
	  std::istringstream iss(rgb_string);
	  iss >> triplet.r;
	  iss >> triplet.g;
	  iss >> triplet.b;
	  colormap.push_back(triplet);
     }

     // read in stat file. 
     ReaderType3DF::Pointer statReader = ReaderType3DF::New();
     statReader->SetFileName(stat_file);
     ImageType3DF::Pointer statPtr = statReader->GetOutput();
     statPtr->Update();

     ReaderType3DFloat::Pointer grayReader = ReaderType3DFloat::New();
     grayReader->SetFileName(gray_file);
     ImageType3DFloat::Pointer grayPtr = grayReader->GetOutput();
     grayPtr->Update();

     // extract slice from stat image
     ImageType3DF::IndexType stat_idx;
     stat_idx.Fill(0);
     stat_idx[dimen] = slice_id;
     ImageType3DF::SizeType stat_size = statPtr->GetLargestPossibleRegion().GetSize();
     stat_size[dimen] = 0;

     ExtFilter::Pointer ext_filter = ExtFilter::New();

     ImageType3DF::RegionType stat_region(stat_idx, stat_size);
     ext_filter->SetExtractionRegion(stat_region);
     ext_filter->SetInput(statPtr);
     ext_filter->Update();
     ImageType2DF::Pointer statSlicePtr = ext_filter->GetOutput();

     // extract slice from gray image
     ImageType3DF::IndexType gray_idx;
     gray_idx.Fill(0);
     gray_idx[dimen] = slice_id;
     ImageType3DF::SizeType gray_size = grayPtr->GetLargestPossibleRegion().GetSize();
     gray_size[dimen] = 0;

     ExtFilter_gray::Pointer ext_filter_gray = ExtFilter_gray::New();

     ImageType3DF::RegionType gray_region(gray_idx, gray_size);
     ext_filter_gray->SetExtractionRegion(gray_region);
     ext_filter_gray->SetInput(grayPtr);
     ext_filter_gray->Update();
     ImageType2DF::Pointer grayUnnormalizedSlicePtr = ext_filter_gray->GetOutput();

     // normalize to 0 to 255. 
     typedef itk::RescaleIntensityImageFilter< ImageType2DF, ImageType2DF > RescaleFilterType;
     RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
     rescaleFilter->SetInput(grayUnnormalizedSlicePtr);
     rescaleFilter->SetOutputMinimum(0);
     rescaleFilter->SetOutputMaximum(255);
     rescaleFilter->Update();
     ImageType2DF::Pointer graySlicePtr = rescaleFilter->GetOutput();

      // create out image
     RGBImageType::IndexType outIdx;
     outIdx.Fill(0);
     RGBImageType::SizeType outSize = graySlicePtr->GetLargestPossibleRegion().GetSize();
     RGBImageType::RegionType outRegion;
     outRegion.SetSize(outSize);
     outRegion.SetIndex(outIdx);
     RGBImageType::Pointer outPtr = RGBImageType::New();
     outPtr->SetRegions(outRegion);
     outPtr->Allocate();
     outPtr->FillBuffer(itk::NumericTraits< RGBImageType::PixelType >::Zero);

     // since ITK can not do the overlap between gray image and a continuous
     // image, we have to manually convert the stat image to rgb image, and
     // overlay the rgb image onto the gray image.

     IteratorTypeRGB outIt(outPtr, outPtr->GetLargestPossibleRegion());
     IteratorType2DF grayIt(graySlicePtr, graySlicePtr->GetLargestPossibleRegion());
     IteratorType2DF statIt(statSlicePtr, statSlicePtr->GetLargestPossibleRegion());
     RGBImageType::PixelType rgb_pixel = itk::NumericTraits< RGBImageType::PixelType >::Zero;
     RGBImageType::PixelType out_pixel = itk::NumericTraits< RGBImageType::PixelType >::Zero;
     RGBImageType::PixelType rgb_bg = itk::NumericTraits< RGBImageType::PixelType >::Zero;     

     if (bgwhite) {
	  rgb_bg[0] = 255;
	  rgb_bg[1] = 255;
	  rgb_bg[2] = 255;
     }

     // get the size of the colormap.
     double colormap_size = colormap.size();

     unsigned col_id = 0;
     for (outIt.GoToBegin(), grayIt.GoToBegin(), statIt.GoToBegin(); 
	  !outIt.IsAtEnd(); ++ grayIt, ++ outIt, ++ statIt) {

	  // this is dangeous since inside the brain, intensity can be zero
	  // (rarely).
	  // outside brain.
	  if (grayIt.Get() == 0 ) {
	       out_pixel[0] = rgb_bg[0];
	       out_pixel[1] = rgb_bg[1];
	       out_pixel[2] = rgb_bg[2];
	  }

	  // inside the brain. 
	  else if ( statIt.Get() < stat_min) {
	       if (clip_low) {
		    // set to same value of statmin.
		    rgb_pixel[0] = colormap[0].r;
		    rgb_pixel[1] = colormap[0].g;
		    rgb_pixel[2] = colormap[0].b;

		    // blending
		    out_pixel[0] = grayIt.Get() * (1-opacity) + rgb_pixel[0] * opacity;
		    out_pixel[1] = grayIt.Get() * (1-opacity) + rgb_pixel[1] * opacity;
		    out_pixel[2] = grayIt.Get() * (1-opacity) + rgb_pixel[2] * opacity;
	       }
	       else {
	       // below stat_min, not shown.
		    out_pixel[0] = grayIt.Get();
		    out_pixel[1] = grayIt.Get();
		    out_pixel[2] = grayIt.Get();
	       }
	  }
	  
	  // above stat_max, show as stat_max.
	  else if (statIt.Get() > stat_max) {
	       if (clip_high) {
		    col_id = unsigned(colormap_size - 1);
		    rgb_pixel[0] = colormap[col_id].r;
		    rgb_pixel[1] = colormap[col_id].g;
		    rgb_pixel[2] = colormap[col_id].b;

		    // blending
		    out_pixel[0] = grayIt.Get() * (1-opacity) + rgb_pixel[0] * opacity;
		    out_pixel[1] = grayIt.Get() * (1-opacity) + rgb_pixel[1] * opacity;
		    out_pixel[2] = grayIt.Get() * (1-opacity) + rgb_pixel[2] * opacity;
	       }
	       else {
		    // not show
		    out_pixel[0] = grayIt.Get();
		    out_pixel[1] = grayIt.Get();
		    out_pixel[2] = grayIt.Get();
	       }
	  }


	  // index of the colors in the color map.
	  else {
	       col_id = unsigned ( (colormap_size-1) * (statIt.Get() - stat_min) / (stat_max - stat_min) );
	       rgb_pixel[0] = colormap[col_id].r;
	       rgb_pixel[1] = colormap[col_id].g;
	       rgb_pixel[2] = colormap[col_id].b;
	       // blending
	       out_pixel[0] = grayIt.Get() * (1-opacity) + rgb_pixel[0] * opacity;
	       out_pixel[1] = grayIt.Get() * (1-opacity) + rgb_pixel[1] * opacity;
	       out_pixel[2] = grayIt.Get() * (1-opacity) + rgb_pixel[2] * opacity;
	  }

	  outIt.Set(out_pixel);

     }

     // write to file
     typedef  itk::ImageFileWriter< RGBImageType  > WriterType;
     WriterType::Pointer writer = WriterType::New();
     writer->SetFileName(out_file);
     writer->SetInput(outPtr);
     writer->Update();
}
