#include "commonalt.h"
#include <utilalt.h>

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     unsigned tsIdx = 0;

     std::string maskfile, infmrifile, outfmrifile, regmaskfile, betafile;
     std::string masking, demean;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "Regress out white matter or CSF mean time course from gray matter.")
	  ("mask,m", po::value<std::string>(&maskfile)->default_value("maskimage.nii"), "input maskimage file")
	  ("regressormask,r", po::value<std::string>(&regmaskfile)->default_value("maskimage.nii"), "mask file for regressor (CSF or white matter). Float number in (0,1)")
	  ("infmri,i", po::value<std::string>(&infmrifile)->default_value("."),
	   "input fmri image file name.")

	  ("outfmri,o", po::value<std::string>(&outfmrifile)->default_value("outfmri.nii.gz"),
	   "output fmri image file name.")
	  ("outbeta,b", po::value<std::string>(&betafile)->default_value("beta.nii.gz"),
	   "beta data volume file name.")

	  ("demean,d", po::value<std::string>(&demean)->default_value("yes"), 
	   "remove the intercept (mean) from the residual on each voxel.")
	  ("verbose,v", po::value<unsigned short>(&verbose)->default_value(0), 
	   "verbose level. 0 for minimal output. 3 for most output.");


     po::options_description cmdline_options;
     cmdline_options.add(generic);

     po::positional_options_description p;
     p.add("in", -1);

     po::variables_map vm;        
     po::store(po::command_line_parser(argc, argv).
	       options(cmdline_options).positional(p).run(), vm);
     po::notify(vm);    

     try {
	  if (vm.count("help") | (argc == 1) ) {
	       std::cout << "Usage: cleanmask [options]\n";
	       std::cout << cmdline_options << "\n";
	       return 0;
	  }
     }

     catch(std::exception& e)
     {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // read input fmri file.
     ReaderType4DFloat::Pointer fmriReader = ReaderType4DFloat::New();
     fmriReader->SetFileName(infmrifile);
     fmriReader->Update();
     ImageType4DFloat::Pointer fmriPtr = fmriReader->GetOutput();
     ImageType4DFloat::IndexType fmriIdx;
     ImageType4DFloat::SizeType fmriSize = fmriPtr->GetLargestPossibleRegion().GetSize();

     // read input mask file.
     ReaderType3DChar::Pointer maskReader = ReaderType3DChar::New();
     maskReader->SetFileName(maskfile);
     maskReader->Update();
     ImageType3DChar::Pointer maskPtr = maskReader->GetOutput();
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType maskIdx;
     

     // read regressor's mask file.
     ReaderType3DChar::Pointer rmaskReader = ReaderType3DChar::New();
     rmaskReader->SetFileName(regmaskfile);
     rmaskReader->Update();
     ImageType3DChar::Pointer rmaskPtr = rmaskReader->GetOutput();
     IteratorType3DChar rmaskIt(rmaskPtr, rmaskPtr->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType rmaskIdx;

     // Create beta image.
     ImageType3DFloat::Pointer betaPtr = ImageType3DFloat::New();
     ImageType3DFloat::IndexType start;
     start[0] = 0;
     start[1] = 0;
     start[2] = 0;

     ImageType3DFloat::SizeType betaSize = maskPtr->GetLargestPossibleRegion().GetSize();

     ImageType3DFloat::RegionType betaRegion;
     betaRegion.SetSize(betaSize);
     betaRegion.SetIndex(start);
     betaPtr->SetRegions(betaRegion);
     betaPtr->Allocate();
     betaPtr->FillBuffer( 0 );
     IteratorType3DFloat betaIt(betaPtr, betaPtr->GetLargestPossibleRegion() );

     // Compute mean time course of the regressor.
     unsigned totalPts = 0;
     vnl_vector<float> regressor(fmriSize[3], 0);

     for (rmaskIt.GoToBegin(); !rmaskIt.IsAtEnd(); ++ rmaskIt) {
	  if (rmaskIt.Get() > 0) {
	       totalPts ++;
	       rmaskIdx = rmaskIt.GetIndex();
	       fmriIdx[0] = rmaskIdx[0];
	       fmriIdx[1] = rmaskIdx[1];
	       fmriIdx[2] = rmaskIdx[2];
	       for (tsIdx = 0; tsIdx < fmriSize[3]; tsIdx ++) {
		    fmriIdx[3] = tsIdx;
		    regressor[tsIdx] = regressor[tsIdx] + fmriPtr->GetPixel(fmriIdx);
	       }
	  }
     }

     regressor = regressor / totalPts;
     regressor = regressor - regressor.mean();
     if (verbose >= 1) {
	  printVnlVector(regressor, fmriSize[3]);
	  printf("total number of voxels in regressor mask: %i\n", totalPts);
     }

     // do the regression work.
     double beta = 0;
     double xy = 0, xx = 0, meany = 0;
     for (maskIt.GoToBegin(), betaIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++ betaIt) {

	  maskIdx = maskIt.GetIndex();
	  fmriIdx[0] = maskIdx[0];
	  fmriIdx[1] = maskIdx[1];
	  fmriIdx[2] = maskIdx[2];	       
	  // compute the mean of y_i over all time points i
	  meany = 0;
	  for (tsIdx = 0; tsIdx < fmriSize[3]; tsIdx ++) {
	       fmriIdx[3] = tsIdx;
	       meany = meany + fmriPtr->GetPixel(fmriIdx) / fmriSize[3];
	  }

	  if (maskIt.Get() > 0) {
	       // compute beta.
	       xy = 0, xx = 0;
	       for (tsIdx = 0; tsIdx < fmriSize[3]; tsIdx ++) {
		    fmriIdx[3] = tsIdx;
		    xy = xy + regressor[tsIdx] * (fmriPtr->GetPixel(fmriIdx) - meany);
		    xx = xx + regressor[tsIdx] * regressor[tsIdx];
	       }
	       beta = xy / xx;
	       betaIt.Set( beta );

	       // compute residual.
	       for (tsIdx = 0; tsIdx < fmriSize[3]; tsIdx ++) {
		    fmriIdx[3] = tsIdx;
		    fmriPtr->SetPixel(fmriIdx, fmriPtr->GetPixel(fmriIdx) - regressor[tsIdx] * beta);
	       }
	  } // maskIt

	  if (demean.compare("yes") == 0) {
	       for (tsIdx = 0; tsIdx < fmriSize[3]; tsIdx ++) {
		    fmriIdx[3] = tsIdx;
		    fmriPtr->SetPixel(fmriIdx, fmriPtr->GetPixel(fmriIdx) - meany);
	       }
	  }
	       
     }

     // write the residual to file.
     WriterType4DFloat::Pointer writer = WriterType4DFloat::New();
     writer->SetInput(fmriPtr);
     writer->SetFileName(outfmrifile);
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


     // write the beta to file.
     WriterType3DFloat::Pointer betaWriter = WriterType3DFloat::New();
     betaWriter->SetInput( betaPtr );
     betaWriter->SetFileName(betafile);
     try 
     { 
     	  betaWriter->Update(); 
     } 
     catch( itk::ExceptionObject & err ) 
     { 
     	  std::cerr << "ExceptionObject caught !" << std::endl; 
     	  std::cerr << err << std::endl; 
     	  return EXIT_FAILURE;
     } 

     return 0;
}
