#include "common_mcem.h"
using std::cout;
using std::string;

// Init random generator as a global object.
twister_base_gen_type mygenerator(42u);
int main (int argc, char* argv[])
{
     unsigned verbose = 0;
     float alpha = 0.1;
     float beta = 0.7;
     int scanIdx = 0;
     int numScans;
     int numClusters = 1;
     unsigned seed = 0;
     int grpLabel = 0, currentLabel = 0, cand = 0;
     float denergy = 0;
     float p_acpt = 0;
     std::string truegrpmap, submap;


     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Given group map, create synthetic subject level's label map.")
	  ("alpha,a", po::value<float>(&alpha)->default_value(0.1), 
	   "alpha parameter, connection between group and subject level.")
	  ("beta,b", po::value<float>(&beta)->default_value(1), 
	   "Set MRF parameter beta at subject level.")
	  ("scan,s", po::value<int>(&numScans)->default_value(100),
	   "Number of scan on MRF.")
	  ("grp,g", po::value<std::string>(&truegrpmap)->default_value("truegrpmap.nii"), "true group map.")
	  ("out,o", po::value<std::string>(&submap)->default_value("submap.nii"), "output subject label map.")
	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("verbose,v", po::value<unsigned int>(&verbose)->default_value(0), 
	   "Verbose level.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: generateimage [options]\n";
	       cout << mydesc << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  cout << e.what() << "\n";
	  return 1;
     }    

     mygenerator.seed(static_cast<unsigned int>(seed));     

     // read in grp image
     ReaderType3DChar::Pointer grpReader = ReaderType3DChar::New();
     grpReader->SetFileName(truegrpmap);
     grpReader->Update();


     // we convert 1-based label to 0-based label here.
     AddConstantToImageFilterType::Pointer minusoneFilter = AddConstantToImageFilterType::New();
     minusoneFilter->SetInput(grpReader->GetOutput() );
     minusoneFilter->SetConstant(-1);
     minusoneFilter->Update();


     ImageType3DChar::Pointer grpPtr = minusoneFilter->GetOutput();
     ImageType3DChar::SizeType grpSize = grpPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DChar grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType grpIdx;     

     // compute number of components.
     typedef itk::MinimumMaximumImageCalculator <ImageType3DChar>
	  ImageCalculatorFilterType;

     ImageCalculatorFilterType::Pointer imageCalculatorFilter
	  = ImageCalculatorFilterType::New ();

     imageCalculatorFilter->SetImage(grpPtr);
     imageCalculatorFilter->ComputeMaximum();
     imageCalculatorFilter->ComputeMinimum();
     numClusters = imageCalculatorFilter->GetMaximum() - imageCalculatorFilter->GetMinimum();

     if (verbose >= 1) {
	  printf("numClusters = %i\n", numClusters);
     }

     // Allocate memory for subject label.
     ImageType3DChar::Pointer subPtr = ImageType3DChar::New();
     ImageType3DChar::IndexType start;
     start.Fill(0);
     ImageType3DChar::SizeType subSize = grpSize; 
     ImageType3DChar::RegionType subRegion;
     subRegion.SetSize(subSize);
     subRegion.SetIndex(start);
     subPtr->SetRegions(subRegion);
     subPtr->Allocate();

     // Define neighborhood iterator
     typedef itk::NeighborhoodIterator< ImageType3DChar> NeighborhoodIteratorType;
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);

     NeighborhoodIteratorType subIt( radius, subPtr, 
				     subPtr->GetLargestPossibleRegion() );

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1}};

     // Uniform integer generator.
     boost::uniform_int<> uni_int(0, numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     // init sub map.
     for (grpIt.GoToBegin(), subIt.GoToBegin(); !grpIt.IsAtEnd(); ++ grpIt, ++ subIt)
     {
	  grpLabel = grpIt.Get();
	  if (grpLabel >= 0) {
	       subIt.SetCenterPixel( roll_die() );
	  }
	  else {
	       subIt.SetCenterPixel( -1 );
	  }
     }

     // real sampling work begin here.
     for (scanIdx = 0; scanIdx < numScans; scanIdx ++) {

	  for (grpIt.GoToBegin(), subIt.GoToBegin(); !grpIt.IsAtEnd(); ++ grpIt, ++ subIt) {
	       grpLabel = grpIt.Get();
	       if (grpLabel >= 0) {
		    currentLabel = subIt.GetCenterPixel();
		    cand = roll_die();

		    // alpha term.
		    denergy = alpha * (int(cand != grpLabel) -  int(currentLabel != grpLabel) );
		    denergy 
			 = denergy
			 + ( int(cand != subIt.GetPixel(xminus)) 
			     - int(currentLabel != subIt.GetPixel(xminus))

			     + int(cand != subIt.GetPixel(xplus))
			     - int(currentLabel != subIt.GetPixel(xplus))

			     + int(cand != subIt.GetPixel(yminus))
			     - int(currentLabel != subIt.GetPixel(yminus))

			     + int(cand != subIt.GetPixel(yplus))
			     - int(currentLabel != subIt.GetPixel(yplus))

			     + int(cand != subIt.GetPixel(zminus))
			     - int(currentLabel != subIt.GetPixel(zminus))

			     + int(cand != subIt.GetPixel(zplus))
			     - int(currentLabel != subIt.GetPixel(zplus)) ) * beta;

		    if (denergy <= 0) {
			 subIt.SetCenterPixel(cand);
		    }
		    else {
			 p_acpt = exp(-denergy);
			 if (uni() < p_acpt) {
			      subIt.SetCenterPixel(cand);
			 }
		    }

	       } // grpLabel >= 0
	  } // for grpIt

	  if ( scanIdx%20 == 0) {
	       printf("scan %i done.\n", scanIdx);
	  }

	  if (verbose >= 1 && scanIdx%20 == 0) {
	       save3dchar(subPtr, submap);
	  }

     } // for scanIdx



     save3dchar(subPtr, submap);
}

