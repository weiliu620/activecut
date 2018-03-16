#include "commonalt.h"
#include "MCModel_init.h"


using std::cout;
using std::endl;

twister_base_gen_type mygenerator(42u);

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     unsigned seed = 0;

     unsigned int numSamples = 10;
     unsigned short kmeansIter = 1;
     int numClusters = 4;

     float alpha = 0.5;
     float beta_g = 0.5;
     float beta_z = 0.5;


     std::string fmriPath, grouplabel, maskimage;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "Initialize MRF by Kmeans.")
	  ("kmeansiter", po::value<unsigned short>(&kmeansIter)->default_value(20),
	   "iteration # for kmeans. ")
	  ("seed", po::value<unsigned int>(&seed)->default_value(0),
	   "Random number generator seed.")
	  ("numClusters,k", po::value<int>(&numClusters)->default_value(4),
	       "Number of labels. Default is 5.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.")

	  ("fmripath", po::value<std::string>(&fmriPath)->default_value("."),
	   "noised image file")
	  ("outgrouplabel,o", po::value<std::string>(&grouplabel)->default_value("grouplabel.nii"),
	   "output labeled group lebel image.")
	  ("mask,m", po::value<std::string>(&maskimage)->default_value("maskimage.nii"), "input maskimage file")
	  ("betag", po::value<float>(&beta_g)->default_value(0.5),
	   "group level pairwise interation term")

	  ("betaz", po::value<float>(&beta_z)->default_value(0.5),
	   "subject level pairwise interation term")
	  ("alpha", po::value<float>(&alpha)->default_value(0.5),
	   "connection between group lebel and individual subject level.");


     po::options_description cmdline_options;
     cmdline_options.add(generic);

     po::positional_options_description p;
     p.add("in", -1);

     po::variables_map vm;        
     po::store(po::command_line_parser(argc, argv).
	       options(cmdline_options).positional(p).run(), vm);
     po::notify(vm);    

     try {

	  if (vm.count("help") | argc == 1) {
	       std::cout << "Usage: segimagealt [options]\n";
	       std::cout << cmdline_options << "\n";
	       return 0;
	  }
     }

     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // init random generator.
     mygenerator.seed(static_cast<unsigned int>(seed));

     boost::filesystem::path fmriPathVar(fmriPath);
     boost::filesystem::directory_iterator fmriPathEnd;

     SeriesReaderType5DFloat::Pointer fmriReader = SeriesReaderType5DFloat::New();

     for (boost::filesystem::directory_iterator fmriPathIt(fmriPathVar); fmriPathIt != fmriPathEnd; ++fmriPathIt) {

     	  fmriReader->AddFileName( (*fmriPathIt).path().string());	  
     	  cout <<  "add " << (*fmriPathIt).path().string() << "\n";
     }

     fmriReader->Update();
     ImageType5DFloat::Pointer imagePtr = fmriReader->GetOutput();
     ImageType5DFloat::SizeType fmriSize = imagePtr->GetLargestPossibleRegion().GetSize();

     // read in mask image
     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(maskimage);
     maskReader->Update();
     ImageType3DChar::Pointer maskPtr = maskReader->GetOutput();

     ImageType3DChar::RegionType maskRegion;
     maskRegion = maskReader->GetOutput()->GetLargestPossibleRegion();
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();

     ImageType3DChar::PointType maskOrigin = maskPtr->GetOrigin();
     ImageType3DChar::SpacingType maskSpacing = maskPtr->GetSpacing();
     ImageType3DChar::DirectionType maskDirection = maskPtr->GetDirection();

     if (maskSize[0] != fmriSize[0] 
     	 || maskSize[1] != fmriSize[1] 
     	 || maskSize[2] != fmriSize[2]) {
     	  cout << "mask image and true label image have different size. Need double check before masking. Exit. " << endl;
	  exit(1);
     }


     // Allocate memory for group level label map.
     ImageType3DChar::Pointer grpPtr = ImageType3DChar::New();
     ImageType3DChar::IndexType grpIdx;
     grpIdx.Fill(0);

     ImageType3DChar::SizeType grpSize;
     grpSize[0] = fmriSize[0];
     grpSize[1] = fmriSize[1];
     grpSize[2] = fmriSize[2];

     ImageType3DChar::RegionType grpRegion;
     grpRegion.SetSize(grpSize);
     grpRegion.SetIndex(grpIdx);
     grpPtr->SetRegions(grpRegion);
     grpPtr->Allocate();
     grpPtr->FillBuffer(-1);

     grpPtr->SetOrigin( grpPtr->GetOrigin() );
     grpPtr->SetSpacing( grpPtr->GetSpacing() ); 
     grpPtr->SetDirection( grpPtr->GetDirection () );

     MCModel mcmodel(imagePtr,
		     grpPtr,
		     numClusters,
     		     verbose);
     mcmodel.printSelf("normal");

     /*********************************************/
     unsigned numSubs = fmriSize[4];
     unsigned repeatIdx = 0; // repeat index of K-means.
     double minMeanSSR = 1e10, meanSSR = 0;
     unsigned subIdx = 0, clsIdx = 0;

     // Define a bestcc to save best cc.
     std::vector < vnl_vector<float> > bestcc(numClusters);
     for (subIdx = 0; subIdx < numSubs; subIdx++) {
	  bestvmm[subIdx].comp.resize(numClusters);
	  for (clsIdx = 0; clsIdx < numClusters; clsIdx ++) {
	       bestvmm[subIdx].comp[clsIdx].mu.set_size(fmriSize[3]);
	  }
     }
     

     for (repeatIdx = 0; repeatIdx < kmeansIter; repeatIdx ++) {
     	  printf("kmeans run  %i begin:\n", repeatIdx+1);

     	  meanSSR = mcmodel.kmeans(imagePtr, grpPtr, maskPtr, 1e-6);

	  if (verbose >= 1) {
	       printf("minMeanSSR = %f, current MeanSSR = %f.\n", minMeanSSR, meanSSR);
	       mcmodel.printSelf("normal");
	  }

     	  if (meanSSR < minMeanSSR) {
     	       minMeanSSR = meanSSR;
     	       // Save best mu in to bestvmm
	       mcmodel.GetVMM(bestvmm);
	       
	       save3dchar(grpPtr, grouplabel);
     	  }
     }

     // Found the best vmm. Restore the best vmm into m_vmm.
     mcmodel.SetVMM(bestvmm);

     // Given best vmm (actually, only mu's are used), estimate labels. 
     mcmodel.estLabelsByMean(imagePtr, maskPtr, grpPtr);     
     /**********************************************************/

     printf("main(): Final minMeanSSR = %f.\n", minMeanSSR);
     save3dchar(grpPtr, grouplabel);
}

