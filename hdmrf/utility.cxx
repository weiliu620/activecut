#include <common.h>

int save_corr(const vnl_vector<double> & corr,
	      std::string out_file,
	      ParType par)
{

     // allocate memory for outPtr
     ImageType2DF::Pointer outPtr = ImageType2DF::New();
     ImageType2DF::IndexType start;
     start.Fill(0);
     ImageType2DF::SizeType outSize;
     outSize[0] = par.n_pts;
     outSize[1] = par.n_pts;
     ImageType2DF::RegionType outRegion;
     outRegion.SetSize(outSize);
     outRegion.SetIndex(start);

     outPtr->SetRegions(outRegion);
     outPtr->Allocate();
     outPtr->FillBuffer(0);
     ImageType2DF::IndexType outIdx;

     for (unsigned r = 0; r < par.n_pts; r ++) {
	  outIdx[0] = r;
	  for (unsigned c = 0; c <= r; c ++) {
	       outIdx[1] = c;
	       // linear idx of the conn r.v.
	       unsigned n = r * par.n_pts + c;
	       outPtr->SetPixel(outIdx, corr[n]);
	  } // c
     } // r
     
     WriterType2DF::Pointer writer = WriterType2DF::New();
	  
     writer->SetInput(outPtr);
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

     std::cout << "save_corr(): File " << out_file << " saved.\n";

     return 0;
}


int save_labelmap(const vnl_matrix<double> & Z,
		  std::string out_file,
		  ParType par)
{
     // allocate memory for outPtr
     ImageType3DF::Pointer outPtr = ImageType3DF::New();
     ImageType3DF::IndexType start;
     start.Fill(0);
     ImageType3DF::SizeType outSize;
     outSize[0] = par.n_pts;
     outSize[1] = par.n_pts;
     outSize[2] = Z.cols();
     ImageType3DF::RegionType outRegion;
     outRegion.SetSize(outSize);
     outRegion.SetIndex(start);

     outPtr->SetRegions(outRegion);
     outPtr->Allocate();
     outPtr->FillBuffer(0);
     ImageType3DF::IndexType outIdx;

     for (unsigned r = 0; r < par.n_pts; r ++) {
	  outIdx[0] = r;
	  for (unsigned c = 0; c <= r; c ++) {
	       outIdx[1] = c;
	       // linear idx of the conn r.v.
	       unsigned n = r * par.n_pts + c;
	       for (unsigned k = 0; k < Z.cols(); k ++) {
		    outIdx[2] = k;
		    outPtr->SetPixel(outIdx, Z(n, k));
	       } // k
	  } // c
     } // r
     
     WriterType3DF::Pointer writer = WriterType3DF::New();
	  
     writer->SetInput(outPtr);
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

     std::cout << "save_labelmap(): File " << out_file << " saved.\n";

     return 0;
}

int print_sparse_matrix(vnl_sparse_matrix<float> mat)
{
     mat.reset();
     unsigned r = 0, c = 0;
     while(mat.next()) {
	  r = mat.getrow();
	  c = mat.getcolumn();
	  printf("[%i %i]: %f\n", r, c, mat.value() );
     }
}

int print_par(ParType par)
{
     for (unsigned k = 0; k < NCOMP; k ++) {
	  printf("comp %i: mu = %.3f, var = %.2f, numPts = %.0f, pi = %f\n", k, par.gmm.comp[k].mu, par.gmm.comp[k].var, par.gmm.comp[k].numPts, par.gmm.comp[k].pi);
     }
}

