#include "commonalt.h"
using std::cout;
using std::string;

// Init random generator as a global object.
twister_base_gen_type mygenerator(42u);
int main (int argc, char* argv[])
{
     unsigned verbose = 0;
     float beta = 0.7;
     int scanIdx = 0;
     int numScans;
     int numClusters = 1;
     unsigned seed = 0;
     std::string observedimage, truegrpmap, maskimage;

     mygenerator.seed(static_cast<unsigned int>(seed));     

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Given mask image and beta, create synthetic group level probability map.")
	  ("beta,b", po::value<float>(&beta)->default_value(1), 
	   "Set MRF parameter beta.")
	  ("numclusters,k", po::value<int>(&numClusters)->default_value(4),
	       "Number of clusters.")
	  ("scan,s", po::value<int>(&numScans)->default_value(500),
	   "Number of scan on MRF.")
	  ("out,o", po::value<std::string>(&truegrpmap)->default_value("grpprobmap.nii"), "output group probability map, a 4D nii file, with dimension [X, Y, Z, K].")
	  ("mask,m", po::value<std::string>(&maskimage)->default_value("maskimage.nii"), "input maskimage file")
	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("verbose,v", po::value<unsigned int>(&verbose)->default_value(0), 
	   "Verbose level.");

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | argc == 1) {
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

     // read in mask image
     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(maskimage);
     maskReader->Update();
     ImageType3DChar::Pointer maskPtr = maskReader->GetOutput();
     ImageType3DChar::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType maskIdx;     

     // Allocate memory for group level prob map.
     itk::VariableLengthVector<float> grpVector(numClusters);
     itk::VariableLengthVector<float> alpha(numClusters);
     ImageType3DVecF::Pointer grpPtr = ImageType3DVecF::New();
     ImageType3DVecF::IndexType start;
     start.Fill(0);
     ImageType3DVecF::SizeType grpSize = maskSize;
     ImageType3DVecF::RegionType grpRegion;
     grpRegion.SetSize(grpSize);
     grpRegion.SetIndex(start);
     grpPtr->SetRegions(grpRegion);
     grpPtr->SetVectorLength(numClusters);
     grpPtr->SetNumberOfComponentsPerPixel(numClusters);
     grpPtr->Allocate();
     grpVector.Fill(0);
     grpPtr->FillBuffer(grpVector);

     IteratorType3DVecF grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );

     // Define neighborhood iterator

     typedef itk::NeighborhoodIterator< ImageType3DVecF> NeighborhoodIteratorType;
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType neiGrpIt(radius, grpPtr, grpPtr->GetLargestPossibleRegion() );

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1}};

     
     // initialize group prob map  with Dirichlet alpha set to 1.
     alpha.Fill(1); // set alpha to all 1 for uniform Dirichlet.
     float grpVecSum = 0;
     

     boost::gamma_distribution<float> gamma_dist(1); 
     boost::variate_generator<twister_base_gen_type&, boost::gamma_distribution<float>  > gammagnr(mygenerator, gamma_dist);  

     for (grpIt.GoToBegin(), maskIt.GoToBegin(); !grpIt.IsAtEnd(); ++ grpIt, ++ maskIt) {
     	  if (maskIt.Get() > 0) {
     	       grpVecSum = 0;
     	       for (unsigned short clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
     		    grpVector[clsIdx] = gammagnr();
     		    grpVecSum += grpVector[clsIdx];
     	       } // clsIdx
     	       grpVector = grpVector / grpVecSum;
     	       grpIt.Set( grpVector );
     	  } // maskIt > 0
     }

     // grpVector = grpPtr->GetPixel(grpIdx);
     // printf("value is [%f, %f]\n", grpVector[0], grpVector[1]);
     
     if (verbose >= 1) {
	  SaveGrpProb(grpPtr, maskPtr, truegrpmap);
     }

     // real sampling work begin here.
     for (scanIdx = 0; scanIdx < numScans; scanIdx ++) {


	  for (neiGrpIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd();++ neiGrpIt, ++ maskIt) {
	       if (maskIt.Get() > 0) {

		    // compute alpha parameter of Dirichlet distribution.
		    // blabla.
		    alpha = ( neiGrpIt.GetPixel(xminus) + neiGrpIt.GetPixel(xplus)
			      + neiGrpIt.GetPixel(yminus) + neiGrpIt.GetPixel(yplus)
			      + neiGrpIt.GetPixel(zminus) + neiGrpIt.GetPixel(zplus) ) / 6;

		    // get unnormalized alpha.
		    alpha = alpha * beta;

		    // add small epsilon to each component of alpha
		    // for numerical stabiility.
		    alpha = alpha + 1;

		    float grpVecSum = 0;
		    for (unsigned short clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
			 // Uniform real random generator.
			 boost::gamma_distribution<float> gamma_dist(alpha[clsIdx]); 
			 boost::variate_generator<twister_base_gen_type&, boost::gamma_distribution<float>  > gammagnr(mygenerator, gamma_dist);  
			 grpVector[clsIdx] = gammagnr() + 1e-5;
			 grpVecSum += grpVector[clsIdx];
		    } // clsIdx
		    grpVector = grpVector / grpVecSum;
		    neiGrpIt.SetCenterPixel( grpVector );
	       } // neiGrpIt
	  } // maskIt > 0
	  if ( scanIdx%20 == 0) {
	       printf("scan %i done.\n", scanIdx);
	       SaveGrpProb(grpPtr, maskPtr, truegrpmap);
	  }
     }

     SaveGrpProb(grpPtr, maskPtr, truegrpmap);
}
