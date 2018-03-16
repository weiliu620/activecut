#include <commonalt.h>
#include <utilalt.h>

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
     unsigned short verbose = 0;
     unsigned tsIdx = 0;
     unsigned short nreg = 0;

     std::string maskfile, infmrifile, outfmrifile, motionparfile, betafile;
     std::string masking, demean;

     // Declare a group of options that will be allowed only on
     // command line
     po::options_description generic("Options can only used at commandline");
     generic.add_options()
	  ("help,h", "Regress out white matter or CSF mean time course from gray matter.")
	  ("mask,m", po::value<std::string>(&maskfile)->default_value("maskimage.nii"), "input maskimage file")
	  ("infmri,i", po::value<std::string>(&infmrifile)->default_value("."),
	   "input fmri image file name.")
	  ("motionpar,p", po::value<std::string>(&motionparfile)->default_value("motionpar.par"),
	   "motion parameter file from motion correction program. Tx6 ascii matrix.")

	  ("nreg,r", po::value<unsigned short>(&nreg)->default_value(6),
	   "number of regressors in motion correction. Should keep this as default value 6.")

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


     // regressor.
     vnl_matrix<float> regressor(fmriSize[3], nreg, 0);
     // read motion parameter file.
     std::string tsValue;
     std::ifstream regStream(motionparfile.c_str() );
     if (regStream.is_open() ) {

	  for (unsigned timeIdx = 0; timeIdx < fmriSize[3]; timeIdx ++) {
	       for (unsigned regIdx = 0; regIdx < nreg; regIdx ++) {
		    regStream >> tsValue;
		    regressor[timeIdx][regIdx] = boost::lexical_cast<float> (tsValue);
	       }

	  }
     }
     else {
	  std::cout << "can not open file "  << motionparfile << std::endl;
     }

     if (verbose >= 1) {
	  printVnlMatrix(regressor, 20);
     }

     // normalize regressors so they have zero mean and unit variance.
     float std = 0;
     for (unsigned regIdx = 0; regIdx < nreg; regIdx ++) {
	  std = regressor.get_column(regIdx).two_norm() / sqrt(fmriSize[3]);
	  regressor.set_column(regIdx, (regressor.get_column(regIdx) - regressor.get_column(regIdx).mean() ) / std);
     }
     
     // Create beta image.
     ImageType4DFloat::Pointer betaPtr = ImageType4DFloat::New();
     ImageType4DFloat::IndexType start;
     start.Fill( 0 );
     ImageType4DFloat::SizeType betaSize = fmriSize;
     betaSize[3] = nreg;
     ImageType4DFloat::RegionType betaRegion;
     betaRegion.SetSize(betaSize);
     betaRegion.SetIndex(start);
     betaPtr->SetRegions(betaRegion);
     betaPtr->Allocate();
     betaPtr->FillBuffer( 0 );
     ImageType4DFloat::IndexType betaIdx;

     // do the regression work.
     double xy = 0, xx = 0, meany = 0;
     vnl_vector<double> beta(nreg, 0);
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {

	  maskIdx = maskIt.GetIndex();
	  fmriIdx[0] = maskIdx[0];
	  fmriIdx[1] = maskIdx[1];
	  fmriIdx[2] = maskIdx[2];
	  betaIdx = fmriIdx;
	  // compute the mean of y_i over all time points i
	  meany = 0;
	  for (tsIdx = 0; tsIdx < fmriSize[3]; tsIdx ++) {
	       fmriIdx[3] = tsIdx;
	       meany = meany + fmriPtr->GetPixel(fmriIdx) / fmriSize[3];
	  }

	  if (maskIt.Get() > 0) {
	       // compute beta.
	       for (unsigned regIdx = 0; regIdx < nreg; regIdx ++) {
		    xx = 0, xy = 0;
		    for (tsIdx = 0; tsIdx < fmriSize[3]; tsIdx ++) {
			 fmriIdx[3] = tsIdx;
			 xy = xy + regressor[tsIdx][regIdx] * (fmriPtr->GetPixel(fmriIdx) - meany);
			 xx = xx + regressor[tsIdx][regIdx] * regressor[tsIdx][regIdx];
		    }
		    beta[regIdx] = xy / xx;

		    // Save beta value.
		    betaIdx[3] = regIdx;
		    betaPtr->SetPixel(betaIdx, beta[regIdx]);
	       } // regIdx
	       
	       // compute residual by removing explained variance from y.
	       for (tsIdx = 0; tsIdx < fmriSize[3]; tsIdx ++) {
		    fmriIdx[3] = tsIdx;
		    float residual = fmriPtr->GetPixel(fmriIdx);
		    for (unsigned regIdx = 0; regIdx < nreg; regIdx ++) {
			 residual -= regressor[tsIdx][regIdx] * beta[regIdx];
		    }
		    fmriPtr->SetPixel(fmriIdx, residual);
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
     WriterType4DFloat::Pointer betaWriter = WriterType4DFloat::New();
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
