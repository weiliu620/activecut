#include <common.h>
#include <utility.h>
twister_base_gen_type mygenerator;

int main(int argc, char* argv[])
{
     unsigned seed = 0,  n_nbrs, n_repeat;
     std::string init_file, out_file;
     unsigned short verbose = 0;

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Given initial label map, apply smooth filter by majority voting. ")

	  ("seed,s", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")

	  ("nbr,n", po::value<unsigned int>(&n_nbrs)->default_value(6),
	   "Number of spatial neighbors. Choose 6, 18 or 26.")

	  ("repeat,r", po::value<unsigned int>(&n_repeat)->default_value(1),
	   "Number of times to apply smoothing filter.")

	   ("init,i", po::value<std::string>(&init_file)->default_value("init.nii.gz"), 
	    "initial label file.")

	   ("out,o", po::value<std::string>(&out_file)->default_value("outlabel.nii.gz"), 
	    "Output smoothed label file.")

	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");


     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: gentestlabel [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // init random generator.
     mygenerator.seed(static_cast<unsigned int>(seed));

     // read in initial label file. Voxel value in [1, K].  Also act as mask
     // since outside brain is zero.
     typedef itk::ImageFileReader< ImageType3DU >  ReaderType3DU;
     ReaderType3DU::Pointer labelReader = ReaderType3DU::New();
     labelReader->SetFileName(init_file);
     labelReader->Update();
     ImageType3DU::Pointer labelPtr = labelReader->GetOutput();

     // get the number of clusters.
     typedef itk::MinimumMaximumImageCalculator <ImageType3DU> ImageCalculatorFilterType;
 
     ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New ();
     imageCalculatorFilter->SetImage(labelReader->GetOutput());
     imageCalculatorFilter->Compute();
     unsigned n_clusters = imageCalculatorFilter->GetMaximum();

     // define neighborhood iterator
     typedef itk::ConstantBoundaryCondition<ImageType3DU>  BoundaryConditionType;
     typedef itk::NeighborhoodIterator< ImageType3DU, BoundaryConditionType > NeighborhoodIteratorType;
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType labelIt(radius, labelPtr, labelPtr->GetLargestPossibleRegion());
     unsigned int nei_set_array[] = {4, 10, 12, 14, 16, 22, // 6 neighborhood
				     1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25, // 18 neighborhood
				     0, 2, 6, 8, 18, 20, 24, 26}; // 26 neighborhood

     if (n_nbrs != 6 && n_nbrs != 18 && n_nbrs != 26) {
	  printf("labelmapsim(): number of neighbors must be 6, 18, or 26. Other values may give inacruate results!\n");
	  exit(1);
     }


     unsigned offset = 0, nei_label;
     
     // repeat the smoothing multiple times. 
     for (unsigned r = 0; r < n_repeat; r ++) {
	  // for (labelIt.GoToBegin(); !labelIt.IsAtEnd(); ++ labelIt) {
	  for (labelIt.GoToEnd(); !labelIt.IsAtBegin(); -- labelIt) {
	       // a label histogram to save the frequency of the labels in the
	       // neighboring voxels.
	       vnl_vector <unsigned> label_hist(n_clusters, 0);
	       // a list of most frequenent lables, in case some labels are equal
	       // frequency.
	       std::vector<unsigned> argmax_list(1, 0);

	       if ( labelIt.GetCenterPixel() > 0) {
		    for (unsigned neiIdx = 0; neiIdx < n_nbrs; neiIdx ++) {
			 offset = nei_set_array[neiIdx];
			 nei_label = labelIt.GetPixel(offset);
			 if (nei_label > 0) {// neighbor is also in mask.
			      label_hist[nei_label - 1] ++;
			 }
		    }

		    if (verbose >=1) {
			 std::cout << "label_hist: " << label_hist << "    ";
		    }

		    // done with collecting frequency. now look for the most frequent
		    // label.

		    // make sure the first run the histogram can not be equal to
		    // tmp_max_count. 
		    int tmp_max_count = -1;
		    for (unsigned k = 0; k < n_clusters; k ++) {
			 if (int(label_hist[k]) > tmp_max_count) {
			      tmp_max_count = label_hist[k];
			      argmax_list.resize(1);
			      argmax_list[0] = k;
			 }
			 else if (int(label_hist[k]) == tmp_max_count) {
			      argmax_list.push_back(k);
			 }
			 else { // smaller freq label, just ignore. 
			 }
		    }

		    if (verbose >=1) {
			 std::cout << " argmax_list: ";
			 for (std::vector<unsigned>::iterator myit = argmax_list.begin(); 
			      myit != argmax_list.end(); ++ myit) {
			      std::cout << *myit << ' ';
			 }
			 std::cout << "    ";
		    }


		    // now I have a list of most frequent labels. Choose one by
		    // uniform sampling.
		    unsigned list_size = argmax_list.size(); // how many equally good labels. 
		    boost::random::uniform_int_distribution<> uni_generator(0, list_size - 1);
		    unsigned random_pos = uni_generator(mygenerator);
		    if (verbose >=1) {
			 std::cout << "    best label: " << argmax_list[random_pos] << "\n";
		    }
		    // update center voxel value. 
		    labelIt.SetCenterPixel(argmax_list[random_pos] + 1);
	       } // labelIt > 0
	  } // for
     } // r 

     // save the label map. 

     typedef itk::ImageFileWriter< ImageType3DU >  WriterType3DU;
     WriterType3DU::Pointer writer = WriterType3DU::New();
     writer->SetInput( labelPtr );
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
     std::cout << "labelmapsim(): File  " << out_file << " saved.\n";

}


	       



