#include <common.h>

std::vector<double> readRow(std::string row);
int draw_roi(ImageType3DFloat::Pointer maskPtr,
	     ImageType3DFloat::IndexType maskIdx,
	     double pixelValue,
	     double radius);

int main(int argc, char* argv[])
{
     std::string eigen_vectors_file, outprefix, coordinate_file, mask_file;
     unsigned short verbose;

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Given the eigen vector, map it to the image.")
	   ("eigenvectors,e", po::value<std::string>(&eigen_vectors_file)->default_value("eigenvectorfile.txt"), 
	    "eigen vector files. K by R matrix.")
	   ("outprefix,o", po::value<std::string>(&outprefix)->default_value("eigen_image_"), 
	    "eigen image nifti file prefix.")

	  ("coord,c", po::value<std::string>(&coordinate_file)->default_value("coordinates.txt"), 
	   "coordinates file. Kx3 file")
	  ("mask,m", po::value<std::string>(&mask_file)->default_value("mask.nii.gz"), 
	   "mask file as a input file.")

	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: geteigenimage [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }    


     ReaderType3DFloat::Pointer maskReader = ReaderType3DFloat::New();
     maskReader->SetFileName(mask_file);
     maskReader->Update();
     ImageType3DFloat::Pointer maskPtr = maskReader->GetOutput();

     IteratorType3DFloat maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());
     ImageType3DFloat::IndexType maskIdx;
     ImageType3DFloat::PointType MNIdx;
     
     // read coordinates.
     std::vector< std::vector<double> > coord;
     std::ifstream coord_stream(coordinate_file.c_str());
     std::vector<double> this_coord(3, 0);
     while(coord_stream >> this_coord[0] >> this_coord[1] >> this_coord[2]) {
     	  coord.push_back(this_coord);
     }

     // read eigenvecotrs.
     std::vector< std::vector<double> > eigen_vectors;
     std::ifstream eigen_vector_stream(eigen_vectors_file.c_str());
     std::string line;
     while (std::getline(eigen_vector_stream, line))
	  eigen_vectors.push_back(readRow(line));


// test transform.
     MNIdx[0] = 0, MNIdx[1] = 0, MNIdx[2] = 0;
     maskPtr->TransformPhysicalPointToIndex(MNIdx, maskIdx);
     printf("MNIdx: %f %f %f, maskIdx: %i %i %i\n", MNIdx[0], MNIdx[1], MNIdx[2], maskIdx[0], maskIdx[1], maskIdx[2]);

// convert MNI coordinates to voxel index.     
     unsigned num_eig_vectors = eigen_vectors[0].size();
     printf("num of vectors: %d\n", num_eig_vectors);
     for (unsigned eigIdx = 0; eigIdx < num_eig_vectors; eigIdx ++) {
	  maskPtr->FillBuffer ( 0);
	  for (unsigned i = 0; i < coord.size(); i ++) {
	       MNIdx[0] = coord[i][0];
	       MNIdx[1] = coord[i][1];
	       MNIdx[2] = coord[i][2];
	       maskPtr->TransformPhysicalPointToIndex(MNIdx, maskIdx);

	       // printf("MNIdx: %f %f %f, maskIdx: %i %i %i\n", MNIdx[0], MNIdx[1], MNIdx[2], maskIdx[0], maskIdx[1], maskIdx[2]);

	       // draw_roi(maskPtr, maskIdx, eigen_vectors[i][eigIdx], 5);
	       draw_roi(maskPtr, maskIdx, 1, 5);
	  }

	  WriterType3DFloat::Pointer writer = WriterType3DFloat::New();
	  writer->SetInput( maskPtr );
	  std::stringstream ss;
	  ss << std::setw(3) << std::setfill('0') << eigIdx;
	  std::string outfile = outprefix + ss.str() + ".nii.gz";
	  writer->SetFileName(outfile);

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
	  std::cout << "applymask(): File  " << outfile << " saved.\n";
     }
     return 0;
}


std::vector<double> readRow(std::string row) {
     std::vector<double> retval;
     std::istringstream is(row);
     double num;
     while (is >> num) {
	  // std::cout << num << " ";
	  retval.push_back(num);
     }
     // std::cout << "\n";
     return retval;
     
}

// draw a ball around the seed.
int draw_roi(ImageType3DFloat::Pointer maskPtr,
		  ImageType3DFloat::IndexType maskIdx,
		  double pixelValue,
		  double radius)
{
     
     radius = int(radius);
     ImageType3DFloat::IndexType thisIdx;
     for (unsigned i = maskIdx[0] - radius; i < maskIdx[0] + radius; i++ ) {
	  for (unsigned j = maskIdx[1] - radius; j < maskIdx[1] + radius; j++ ) {
	       for (unsigned k = maskIdx[2] - radius; k < maskIdx[2] + radius; k++ ) {
		    if ( (i-maskIdx[0])*(i-maskIdx[0]) + (j-maskIdx[1])*(j-maskIdx[1]) + (k-maskIdx[2])*(k-maskIdx[2]) < radius*radius ) {
			 thisIdx[0] = i;
			 thisIdx[1] = j;
			 thisIdx[2] = k;
			 maskPtr->SetPixel(thisIdx, pixelValue);
		    }
	       } // k
	  } // j
     } // i

     return 0;
}
