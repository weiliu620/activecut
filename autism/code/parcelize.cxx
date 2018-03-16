#include <common.h>
#include <utility.h>
namespace po = boost::program_options;
struct CoType {
     double x;
     double y;
     double z;
};

// convert RAS coordinates to LPS.
ImageType4DF::PointType to_LPS(const ImageType4DF::PointType  in_point);
// LPS coordinates to RAS.
ImageType4DF::PointType to_RAS(const ImageType4DF::PointType  in_point);
// if lower bound is larger than higher bound, swap them.
int sort_bound(unsigned & lb, unsigned & hb);

int main(int argc, char* argv[])
{
     std::string input_image_file, out_file, mask_file, seeds_file, avgbold_file;
     unsigned verbose = 0;
     double radius = 0;
     
     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Given the seeds MNI cooridnates, parcelate the voxels in to regions..")

	  ("data,d", po::value<std::string>(&input_image_file),
	   "input fmri images.")

	  ("seeds,s", po::value<std::string>(&seeds_file),
	   "a text file with each row the MNI coordinates of the seeds.")

	  ("radius,r", po::value<double>(&radius)->default_value(3),
	   "Search radius for averaging the BOLD signal around the seed. Currently averagign is within a cube with 2r as the edge length.  Unit is mm.")

	  ("out,o", po::value<std::string>(&out_file)->default_value("features.txt"),
	   "output correlation txt file.")

	  ("avgbold,a", po::value<std::string>(&avgbold_file)->default_value("avgbold.txt"),
	   "output NxT average bold signal file. N is number of ROIs, and T is time series length. ")

	  ("verbose,v", po::value<unsigned>(&verbose)->default_value(0),
	   "verbose level in [0, 3]. ")
	  ;

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: parcelize [options]\n";
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
     ImageType4DF::IndexType dataIdx;
     dataIdx.Fill(0);

     // read seed MNI coordinates.
     
     std::vector<ImageType4DF::PointType> seeds;
     std::ifstream seeds_stream(seeds_file.c_str() );
     ImageType4DF::PointType mniIdx;
     mniIdx.Fill(0);
     std::string line;
     while (getline(seeds_stream, line)) {
	  std::istringstream ss(line);
	  ss >>  mniIdx[0] >>  mniIdx[1] >> mniIdx[2];
	  seeds.push_back(mniIdx);
     }

     // compute average time series for each MNI seed.
     vnl_matrix<double> bold(seeds.size(), dataSize[3], 0);
     vnl_vector<unsigned> n_pts(seeds.size(), 0);

     dataIdx.Fill(0);
     dataPtr->TransformIndexToPhysicalPoint(dataIdx, mniIdx);
     mniIdx = to_RAS(mniIdx);
     // std::cout << "mniIdx: " << mniIdx << " <- dataIdx: " << dataIdx << std::endl;

     unsigned il = 0, ih = 0, jl = 0, jh = 0, kl = 0, kh = 0; // lower/upper limit.
     ImageType4DF::IndexType nbr_idx;
     nbr_idx.Fill(0);
     double vox_radius = 0;
     for (unsigned n = 0; n < seeds.size(); n ++) { // seed index.
	  // get to know the cube boundary in voxel domain.
	  mniIdx = seeds[n];
	  mniIdx[0] -=  radius;
	  mniIdx[1] -=  radius;
	  mniIdx[2] -=  radius;
	  dataPtr->TransformPhysicalPointToIndex(to_LPS(mniIdx), dataIdx);
	  il = dataIdx[0], jl = dataIdx[1], kl = dataIdx[2];

	  mniIdx = seeds[n];
	  mniIdx[0] += radius;
	  mniIdx[1] += radius;
	  mniIdx[2] += radius;
	  dataPtr->TransformPhysicalPointToIndex(to_LPS(mniIdx), dataIdx);
	  ih = dataIdx[0], jh = dataIdx[1], kh = dataIdx[2];
	  sort_bound(il, ih);
	  sort_bound(jl, jh);
	  sort_bound(kl, kh);
	  // printf("lower: [%i %i %i]. Higher bound: [%i %i %i]\n", il, jl, kl, ih, jh, kh);

	  mniIdx = seeds[n];
	  dataPtr->TransformPhysicalPointToIndex(to_LPS(mniIdx), dataIdx);
	  // std::cout << "mniIdx: " << mniIdx << " -> dataIdx: " << dataIdx << std::endl;
	  
	  vox_radius = (ih - il) * 0.5; // a shortcut to get radius in voxel domain.

	  // make sure <= in case il=ih
	  for (nbr_idx[0] = il; nbr_idx[0] <= ih; nbr_idx[0] ++) {
	       for (nbr_idx[1] = jl; nbr_idx[1] <= jh; nbr_idx[1] ++) {
		    for (nbr_idx[2] = kl; nbr_idx[2] <= kh; nbr_idx[2] ++) {
			 // in the sphere?
			 // std::cout << "nbr_idx: " << nbr_idx << ".   dataIdx: " << dataIdx << std::endl;
			 // make sure <= in case vol_radius = 0;
			 if (pow(nbr_idx[0]-dataIdx[0], 2) + pow(nbr_idx[1]-dataIdx[1], 2) + pow(nbr_idx[2]-dataIdx[2], 2) <= pow(vox_radius, 2) ) {
			      n_pts[n] ++;
			      for ( nbr_idx[3] = 0; nbr_idx[3] < dataSize[3]; nbr_idx[3] ++) {
				   bold(n, nbr_idx[3]) += dataPtr->GetPixel(nbr_idx);
			      }
			 } // in sphere.
		    } // [2]
	       } // [1]
	  } // [0]
     } // sample n

     // get the mean from the sum.
     for (unsigned n = 0; n < seeds.size(); n ++) { // seed index.
	  bold.scale_row(n, 1/ double(n_pts[n]));
     }	  

     // begin compute correlation. First substract the mean.
     for (unsigned n = 0; n < seeds.size(); n ++) { // seed index.
	  bold.set_row(n, bold.get_row(n) - bold.get_row(n).mean () );
     }

     // divide by the variance.
     for (unsigned n = 0; n < seeds.size(); n ++) { // seed index.
	  double var = bold.get_row(n).squared_magnitude() / bold.cols();
	  if (var > 0 ) {
	       bold.scale_row(n, sqrt(1/var));
	  }
	  else {
	       // all elements must be all zero. Do nothing.
	  }

	  if (n==1) {
	       printf("res = %f\n", bold.get_row(n).squared_magnitude());
	       // vnl_matlab_print(vcl_cout, bold.get_row(n));
	  }
     }
     
     // compute correlation matrix.
     vnl_matrix<double> corr_mat = bold * bold.transpose() / bold.cols();

     // save to file
     std::ofstream outstr;
     outstr.open(out_file.c_str());
     corr_mat.print(outstr);
     outstr.close();

     outstr.open(avgbold_file.c_str());
     bold.print(outstr);
     outstr.close();

     return 0;
}

// ITK use RAS->LPS coordinates, but MNI use RAS. Need add minus sign
// in front of x and y of MNI coordinates to convert it to LPI.
ImageType4DF::PointType to_LPS(const ImageType4DF::PointType  in_point)
{
     ImageType4DF::PointType out_point = in_point;
     out_point[0] = - out_point[0];
     out_point[1] = - out_point[1];
     return out_point;
}

ImageType4DF::PointType to_RAS(const ImageType4DF::PointType  in_point)
{
     ImageType4DF::PointType out_point = in_point;
     out_point[0] = - out_point[0];
     out_point[1] = - out_point[1];
     return out_point;
}

int sort_bound(unsigned & lb, unsigned & hb)
{
     if (lb > hb) {
	  unsigned lbt = lb;
	  lb = hb;
	  hb = lbt;
     }
     return 0;
}
