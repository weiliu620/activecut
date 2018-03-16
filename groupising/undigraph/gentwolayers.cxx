#include "common_mcem.h"

// Init random generator as a global object.
twister_base_gen_type mygenerator(42u);
typedef itk::ConstantBoundaryCondition< ImageType3DChar >  MyBoundCondType;
typedef itk::NeighborhoodIterator< ImageType3DChar, MyBoundCondType > MyNeiItType;
// model parameters.
struct ParType 
{
     unsigned numClusters;
     unsigned numSubs;
     float alpha;
     float betag;
     float betaz;
     unsigned short verbose;
};
	  
int InitRand(ImageType3DChar::Pointer imPtr, 
	     ImageType3DChar::Pointer maskPtr,
	     ParType & par);

int SamplingSub(ImageType3DChar::Pointer grpPtr,
		ImageType3DChar::Pointer subPtr,
		ImageType3DChar::Pointer maskPtr,
		ParType & par);

int SamplingGrp(ImageType3DChar::Pointer grpPtr,
		ImageType3DChar::Pointer maskPtr,
		std::vector<ImageType3DChar::Pointer>  & subVec,
		ParType & par);


int main (int argc, char* argv[])
{
     ParType par;
     int num_scan;

     unsigned seed = 0;

     std::string trueimage, maskimage;
     std::string outGrpLabel, outSubPrefix;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "create group label map.")
	  ("betag", po::value<float>(&par.betag)->default_value(0.5),
	   "group level pairwise interation term")

	  ("betaz", po::value<float>(&par.betaz)->default_value(0.5),
	   "subject level pairwise interation term")
	  ("alpha", po::value<float>(&par.alpha)->default_value(0.5),
	   "connection between group lebel and individual subject level.")

	  ("numClusters,k", po::value<unsigned>(&par.numClusters)->default_value(6),
	       "Number of clusters. Default is 6.")

	  ("numsubs,J", po::value<unsigned>(&par.numSubs)->default_value(10),
	       "Number of subjects.")

	  ("scan,s", po::value<int>(&num_scan)->default_value(500),
	   "Number of scan on MRF.")

	  ("mask,m", po::value<std::string>(&maskimage)->default_value("maskimage.nii"), "input maskimage file")

	  ("grouplabel,g", po::value<std::string>(&outGrpLabel)->default_value("outgrplabel.nii"), "output group label file. nii format.")

	  ("outsubprefix", po::value<std::string>(&outSubPrefix)->default_value("outSub"), "output subject label map file prefix.")

	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")


	  ("verbose,v", po::value<unsigned short>(&par.verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");



     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: generateimage [options]\n";
	       std::cout << mydesc << "\n";
	       return (0);
	  }
     }
     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // init random generator.
     mygenerator.seed(static_cast<unsigned int>(seed));

     // read in mask image
     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(maskimage);
     maskReader->Update();
     ImageType3DChar::Pointer maskPtr = maskReader->GetOutput();
     ImageType3DChar::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());


     // create image buffer for group
     ImageType3DChar::IndexType start;
     start.Fill(0);

     ImageType3DChar::RegionType grpRegion;
     grpRegion.SetSize(maskSize);
     grpRegion.SetIndex(start);
     ImageType3DChar::Pointer grpPtr = ImageType3DChar::New();
     grpPtr->SetRegions(grpRegion);
     grpPtr->Allocate();
     grpPtr->FillBuffer(-1);

     // create buffers for subjects.
     unsigned subIdx = 0;
     std::vector<ImageType3DChar::Pointer> subVec(par.numSubs);
     ImageType3DChar::RegionType subRegion;
     subRegion.SetSize(maskSize);
     subRegion.SetIndex(start);
     for (subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  subVec[subIdx] = ImageType3DChar::New();
	  subVec[subIdx]->SetRegions(subRegion);
	  subVec[subIdx]->Allocate();
	  subVec[subIdx]->FillBuffer(-1);
     }


     // init subject label map.
     for (subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  InitRand(subVec[subIdx], maskPtr, par);
     }
     
     // init group label map.
     InitRand(grpPtr, maskPtr, par);

     for (unsigned scanIdx = 0; scanIdx < num_scan; scanIdx ++) {
	  // sampling sub.
// #pragma omp parallel for
	  for (subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	       SamplingSub(grpPtr, subVec[subIdx], maskPtr, par);
	  }
	  // sampling group.
	  SamplingGrp(grpPtr, maskPtr, subVec, par);

	  if (par.verbose >= 1 && scanIdx%10 == 0) {
	       printf("scan %i done.\n", scanIdx);
	  } // verbose
     }

     // Save final resutls for both group label map and subjects label
     // map. 
     save3dchar(grpPtr, outGrpLabel);
     for (subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  std::string subFileNum = boost::lexical_cast<std::string> (subIdx+1);
	  if (subIdx < 9) {
	       subFileNum.insert(0, "0");
	  }
	  std::string thisSubFilename(outSubPrefix);
	  // add sub number.
	  thisSubFilename.append(subFileNum);
	  thisSubFilename.append(".nii.gz");
	  save3dchar(subVec[subIdx], thisSubFilename);
     }

}

int InitRand(ImageType3DChar::Pointer imPtr, 
	     ImageType3DChar::Pointer maskPtr,
	     ParType & par)
{
     // Uniform integer generator.
     boost::uniform_int<> uni_int(0, par.numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     IteratorType3DChar imIt(imPtr, imPtr->GetRequestedRegion());
     //mask
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );


     for (maskIt.GoToBegin(), imIt.GoToBegin(); !imIt.IsAtEnd(); ++imIt, ++maskIt) {
	  if (maskIt.Get() > 0) {
	       imIt.Set(roll_die());
	  }
     }
}

int SamplingSub(ImageType3DChar::Pointer grpPtr,
		ImageType3DChar::Pointer subPtr,
		ImageType3DChar::Pointer maskPtr,
		ParType & par)

{

     double denergy= 0, denergyPrior = 0, denergyLL=0;
     int cand;
     double p_acpt = 0; // probability to accept the candidate sample.
     signed int currentLabel = 0;
     unsigned grpLabel = 0;

     // Uniform integer generator.
     boost::uniform_int<> uni_int(0, par.numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);


     //mask
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     // group.
     IteratorType3DChar grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );

     // Define neighborhood iterator
     MyBoundCondType constCondition;
     constCondition.SetConstant(-1);
     MyNeiItType::RadiusType radius;
     radius.Fill(1);

     MyNeiItType neiSampleIt( radius, subPtr, subPtr->GetRequestedRegion() );
     neiSampleIt.OverrideBoundaryCondition(&constCondition);

     MyNeiItType::OffsetType xplus = {{1,0, 0}};
     MyNeiItType::OffsetType xminus = {{-1, 0, 0}};
     MyNeiItType::OffsetType yplus = {{0, 1, 0}};
     MyNeiItType::OffsetType yminus = {{0, -1, 0}};
     MyNeiItType::OffsetType zplus = {{0, 0, 1}};
     MyNeiItType::OffsetType zminus = {{0, 0, -1}};



     // sampling Markov Random Field. 
     for (neiSampleIt.GoToBegin(), grpIt.GoToBegin(), maskIt.GoToBegin(); 
	  !neiSampleIt.IsAtEnd(); 
	  ++ neiSampleIt, ++ grpIt, ++ maskIt) {

	  if (maskIt.Get() > 0) {
	       currentLabel = neiSampleIt.GetCenterPixel();
	       cand = roll_die();
	       // alpha (group) term. Negative log of probabililty.
	       grpLabel = grpIt.Get();
	       denergyPrior = par.alpha * ( int(cand != grpLabel) - int(currentLabel != grpLabel) );
	       denergyLL 
		    = int(cand != neiSampleIt.GetPixel(xminus))
		    - int(currentLabel != neiSampleIt.GetPixel(xminus))

		    + int(cand != neiSampleIt.GetPixel(xplus)) 
		    - int(currentLabel != neiSampleIt.GetPixel(xplus))

		    + int(cand != neiSampleIt.GetPixel(yminus)) 
		    - int(currentLabel != neiSampleIt.GetPixel(yminus))

		    + int(cand != neiSampleIt.GetPixel(yplus)) 
		    - int(currentLabel != neiSampleIt.GetPixel(yplus))

		    + int(cand != neiSampleIt.GetPixel(zminus)) 
		    - int(currentLabel != neiSampleIt.GetPixel(zminus))

		    + int(cand != neiSampleIt.GetPixel(zplus)) 
		    - int(currentLabel != neiSampleIt.GetPixel(zplus));

	       // beta (sub level pairwise interation) term.
	       denergyLL = par.betaz * denergyLL;

	       denergy = denergyPrior + denergyLL;

	       // if energy change less than zero, just accept
	       // candidate. otherwise accept with exp(- energy
	       // change).
			 
	       if (denergy <= 0) {
		    neiSampleIt.SetCenterPixel(cand);
	       }
	       else {
		    p_acpt = exp(-denergy);
		    if (uni() < p_acpt) {
			 neiSampleIt.SetCenterPixel(cand);
		    }
	       }
	  } // in mask
     } // iterators.

     
}



int SamplingGrp(ImageType3DChar::Pointer grpPtr,
			 ImageType3DChar::Pointer maskPtr,
			 std::vector<ImageType3DChar::Pointer>  & subVec,
			 ParType & par)
{
     double p_acpt = 0;
     double denergy= 0, denergyPrior = 0, denergyBeta=0;
     int cand = 1;
     int currentLabel = 0;

     //mask
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     // Define neighborhood iterator
     MyBoundCondType constCondition;
     constCondition.SetConstant(-1);
     MyNeiItType::RadiusType radius;
     radius.Fill(1);

     MyNeiItType neiGrpIt(radius, grpPtr, grpPtr->GetLargestPossibleRegion() );
     neiGrpIt.OverrideBoundaryCondition(&constCondition);

     MyNeiItType::OffsetType xplus = {{1,0, 0}};
     MyNeiItType::OffsetType xminus = {{-1, 0, 0}};
     MyNeiItType::OffsetType yplus = {{0, 1, 0}};
     MyNeiItType::OffsetType yminus = {{0, -1, 0}};
     MyNeiItType::OffsetType zplus = {{0, 0, 1}};
     MyNeiItType::OffsetType zminus = {{0, 0, -1}};

     ImageType3DChar::SizeType grpSize = grpPtr->GetLargestPossibleRegion().GetSize();

     // create histogram data block.
     itk::VariableLengthVector<unsigned short> histVector(par.numClusters);

     ImageType3DVecUS::Pointer histPtr = ImageType3DVecUS::New();
     ImageType3DVecUS::RegionType histRegion = grpPtr->GetLargestPossibleRegion();
     histPtr->SetRegions(histRegion);
     histPtr->SetNumberOfComponentsPerPixel(par.numClusters);
     histPtr->SetVectorLength(par.numClusters);
     histPtr->Allocate();

     histVector.Fill(0);
     histPtr->FillBuffer(histVector);
     IteratorType3DVecUS histIt(histPtr, histRegion);


     

     // compute histogram
     for (unsigned short subIdx = 0; subIdx < par.numSubs; subIdx ++) {
	  ConstIteratorType3DChar sampleIt(subVec[subIdx], subVec[subIdx]->GetLargestPossibleRegion() );

	  for (sampleIt.GoToBegin(), histIt.GoToBegin(); 
	       !histIt.IsAtEnd(); ++ histIt, ++sampleIt) {
	       if (sampleIt.Get() >= 0) {
		    histVector = histIt.Get();
		    histVector[sampleIt.Get()] ++;
		    histIt.Set(histVector);
	       }
	  } // for sampelIt.
     } // for subIdx

     // Uniform integer generator.
     boost::uniform_int<> uni_int(0, par.numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     
     for (neiGrpIt.GoToBegin(), histIt.GoToBegin(), maskIt.GoToBegin(); 
	  !histIt.IsAtEnd(); ++ histIt, ++ neiGrpIt, ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       currentLabel = neiGrpIt.GetCenterPixel();
	       cand = roll_die();
	       denergyPrior = par.alpha * (histIt.Get()[currentLabel] - histIt.Get()[cand]);
	       denergyBeta = 
	       
		    int(cand != neiGrpIt.GetPixel(xminus))
		    - int(currentLabel != neiGrpIt.GetPixel(xminus))

		    + int(cand != neiGrpIt.GetPixel(xplus)) 
		    - int(currentLabel != neiGrpIt.GetPixel(xplus))

		    + int(cand != neiGrpIt.GetPixel(yminus)) 
		    - int(currentLabel != neiGrpIt.GetPixel(yminus))

		    + int(cand != neiGrpIt.GetPixel(yplus)) 
		    - int(currentLabel != neiGrpIt.GetPixel(yplus))

		    + int(cand != neiGrpIt.GetPixel(zminus)) 
		    - int(currentLabel != neiGrpIt.GetPixel(zminus))

		    + int(cand != neiGrpIt.GetPixel(zplus)) 
		    - int(currentLabel != neiGrpIt.GetPixel(zplus));	       
	       
	       denergyBeta = denergyBeta * par.betag;

	       denergy = denergyPrior + denergyBeta;


	       if (denergy <= 0) {
		    neiGrpIt.SetCenterPixel( cand );
	       }
	       else {
		    p_acpt = exp(-denergy);
		    if (uni() < p_acpt) {
			 neiGrpIt.SetCenterPixel( cand );
		    }
	       }
	  } // maskIt > 0
     } // neiGrpIt
     return 0;     
}


