#include <common.h>
#include <utility.h>
namespace po = boost::program_options;

int compute_corr(vnl_vector<double> & C,
		 ImageType4DF::Pointer fmriPtr, 
		 ImageType3DU::Pointer maskPtr,
		 std::vector<ImageType3DU::IndexType> & to_xyz,
		 ImageType3DU::Pointer to_linPtr,
		 ParType par);

int build_graph(ImageType3DU::Pointer maskPtr,
		vnl_sparse_matrix<float> & G,
		ImageType3DU::Pointer to_linPtr,
		ParType par);

int add_edge(unsigned n1, // voxel linear index.  
	     unsigned n2,
	     unsigned n_pts,
	     vnl_sparse_matrix<float> & G);

int E_step(const vnl_vector<double> & C,
	   vnl_matrix<double> & Z,
	   vnl_sparse_matrix<float> & G,
	   ParType par,
     	   double beta);

int M_step(const vnl_vector<double> & C,
	   const vnl_matrix<double> & Z,
	   ParType & par);

int main(int argc, char* argv[])
{
     std::string fmri_file, mask_file, label_file;
     unsigned maxem = 0;
     
     // program options.
     ParType par;
     
     // init parameters.
     par.gmm.comp.resize(NCOMP);
     par.gmm.n_comp = NCOMP;

     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "connectivity estimation by variation Bayesian and 6D MRF.")
	  ("maxem,e", po::value<unsigned>(&maxem)->default_value(10),
	   "max numbef of EM iteration.")
	  ("beta,b", po::value<double>(&par.beta)->default_value(0.5),
	   "smooth parameter.")

	  ("mu0", po::value<double>(&par.gmm.comp[0].mu)->default_value(0),
	   "mu of comp 0.")
	  ("mu1", po::value<double>(&par.gmm.comp[1].mu)->default_value(0.5),
	   "mu of comp 1.")
	  ("mu2", po::value<double>(&par.gmm.comp[2].mu)->default_value(1.0),
	   "mu of comp 1.")

	  ("fmri,f", po::value<std::string>(&fmri_file)->default_value("fmri.nii.gz"),
	   "fmri file")
	  ("mask,m", po::value<std::string>(&mask_file)->default_value("mask.nii.gz"),
	   "mask file.")

	  ("label,l", po::value<std::string>(&label_file)->default_value("Z.nii.gz"),
	   "hidden labels.")

	  ("verbose,v", po::value<unsigned>(&par.verbose)->default_value(0),
	   "verbose level in [0, 3]. ")
	  ;

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: grabcut [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // read fmri images.
     ReaderType4DF::Pointer fmriReader = ReaderType4DF::New();
     fmriReader->SetFileName(fmri_file);
     fmriReader->Update();
     ImageType4DF::Pointer fmriPtr = fmriReader->GetOutput();

     // read mask.
     ReaderType3DU::Pointer maskReader = ReaderType3DU::New();
     maskReader->SetFileName(mask_file);
     maskReader->Update();
     ImageType3DU::Pointer maskPtr = maskReader->GetOutput();

     // compute num of points in gray matter.
     typedef itk::StatisticsImageFilter<ImageType3DU> StatisticsImageFilterType;
     StatisticsImageFilterType::Pointer statisticsImageFilter
          = StatisticsImageFilterType::New ();
     statisticsImageFilter->SetInput(maskPtr);
     statisticsImageFilter->Update();
     par.n_pts = statisticsImageFilter->GetSum();
     printf("Total number of gray matter voxels: %i\n", par.n_pts);

     // extract BOLD 
     std::vector<ImageType3DU::IndexType> to_xyz(par.n_pts); // linear index to xyz voxel coordinates.
     vnl_vector<double> corr(par.n_pts * par.n_pts, 0);
     ImageType3DU::Pointer to_linPtr = ImageType3DU::New();
     compute_corr(corr, fmriPtr, maskPtr, to_xyz, to_linPtr, par);
     save_corr(corr, "corr.nii.gz", par);

     for (unsigned k = 0; k < par.gmm.n_comp; k ++) {
	  par.gmm.comp[k].var = 0.01;
	  par.gmm.comp[k].pi = 1/double(par.gmm.n_comp);
     }

     vnl_sparse_matrix<float> G;
     build_graph(maskPtr, G, to_linPtr, par);
     // print_sparse_matrix(G);

     // define a hidden connectivity variables.
     vnl_matrix<double> Z(par.n_pts * par.n_pts, NCOMP, 0); 

     E_step(corr, Z, G, par, 0);
     M_step(corr, Z, par);
     print_par(par);

     unsigned em_iter = 0;
     while (em_iter < maxem) {
	  printf("EM iter %i begin.\n", em_iter);
     	  em_iter ++;
     	  E_step(corr, Z, G, par, par.beta);
     	  M_step(corr, Z, par);
	  print_par(par);
     }
     save_labelmap(Z, label_file, par);
     return 0;
}

int compute_corr(vnl_vector<double> & C,
		 ImageType4DF::Pointer fmriPtr, 
		 ImageType3DU::Pointer maskPtr,
		 std::vector<ImageType3DU::IndexType> & to_xyz,
		 ImageType3DU::Pointer to_linPtr,
		 ParType par)
{
     IteratorType4DF fmriIt(fmriPtr, fmriPtr->GetLargestPossibleRegion());
     ImageType4DF::SizeType fmriSize = fmriPtr->GetLargestPossibleRegion().GetSize();
     ImageType4DF::IndexType fmriIdx;

     IteratorType3DU maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());
     ImageType3DU::IndexType maskIdx;

     // allocate memory for to_linPtr
     ImageType3DU::RegionType maskRegion = maskPtr->GetLargestPossibleRegion();
     to_linPtr->SetRegions(maskRegion);
     to_linPtr->Allocate();
     to_linPtr->FillBuffer(0);

     vnl_matrix<double> bold(par.n_pts, fmriSize[3]);
     to_xyz.resize(par.n_pts);
     C.set_size(par.n_pts * par.n_pts);

     unsigned n = 0;
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       maskIdx = maskIt.GetIndex();
	       fmriIdx[0] = maskIdx[0];
	       fmriIdx[1] = maskIdx[1];
	       fmriIdx[2] = maskIdx[2];
	       
	       to_xyz[n] = maskIdx;
	       to_linPtr->SetPixel(maskIdx, n);

	       for (fmriIdx[3] = 0; fmriIdx[3] < fmriSize[3]; fmriIdx[3] ++) {
		    bold(n, fmriIdx[3]) = fmriPtr->GetPixel(fmriIdx);
	       }

	       // substract the mean.
	       bold.set_row(n, bold.get_row(n) - bold.get_row(n).mean () );

	       // divide by the variance.
	       double var = bold.get_row(n).squared_magnitude() / bold.cols();
	       if (var > 0 ) {
		    bold.scale_row(n, sqrt(1/var));
	       }
	       else {
		    // all elements must be all zero. Do nothing.
	       }

	       n++;
	  } // > 0
     } // maskIt
     
     vnl_matrix<double> corr_mat = bold * bold.transpose() / bold.cols();

     // Fisher transform. 
     for (unsigned r = 0; r < corr_mat.rows(); r++) {
	  for (unsigned c = 0; c < corr_mat.cols(); c++) {
	       double corr = corr_mat(r, c);
	       if (corr > 0.999) {
		    corr_mat(r, c) = 3.8;
	       }
	       else if (corr < -0.999) {
		    corr_mat(r, c) = -3.8;
	       }
	       else {
		    corr_mat(r,c) = 0.5 * log( (1+corr) / (1 - corr));
	       }
	  } // c
     }// r

     // covert to vector.
     for (unsigned r = 0; r < corr_mat.rows(); r++) {
	  for (unsigned c = 0; c < corr_mat.cols(); c++) {
	       C[r * par.n_pts + c] = corr_mat(r, c);
	  }
     }
     return 0;
}

int build_graph(ImageType3DU::Pointer maskPtr,
		vnl_sparse_matrix<float> & G,
		ImageType3DU::Pointer to_linPtr,
		ParType par)
{
     // define neighborhood iterator
     typedef itk::ConstantBoundaryCondition<ImageType3DU>  BoundaryConditionType;
     typedef itk::NeighborhoodIterator< ImageType3DU, BoundaryConditionType > NeighborhoodIteratorType;
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);

     NeighborhoodIteratorType maskIt(radius, maskPtr, maskPtr->GetLargestPossibleRegion());
     ImageType3DU::IndexType maskIdx, nbrIdx;     
     unsigned int nei_set_array[] = {4, 10, 12, 14, 16, 22};

     // we need to set the sparse matrix size here. Otherwise the assigment will
     // fail. Not sure why this happens. (ideally it should just add an entry).
     G.set_size(par.n_pts * par.n_pts, par.n_pts * par.n_pts);

     unsigned v1 = 0, v2 = 0; // voxel linear index.
     
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.GetCenterPixel() > 0) {
	       maskIdx = maskIt.GetIndex();
	       v1 = to_linPtr->GetPixel(maskIdx);
	       // std::cout << "build_graph(): maskIdx " << maskIdx << "--> " << v1 << std::endl;
	       // first we need explicitly add edge with itself. Otherwise the
	       // sparse matrix does not have entry at diagonal.
	       // add_edge(v1, v1, par.n_pts, G);
	       for (unsigned nbr_id = 0; nbr_id < 6; nbr_id ++) {
		    unsigned offset = nei_set_array[nbr_id];
		    nbrIdx =  maskIt.GetIndex(offset);
		    v2 = to_linPtr->GetPixel(nbrIdx);
		    // make sure (v1, v2) is unordered pair and only count once.
		    if ( (maskIt.GetPixel(offset) > 0) && (v2 > v1) ) {
			 // add edge to the graph for all connectivity variables.
			 // std::cout << "    build_graph(): nbrIdx " << nbrIdx << "--> " << v2 << std::endl;
			 add_edge(v1, v2, par.n_pts, G);

			 } // neighbor is a good one.
	       } // nbr_id
	  } // > 0
     }

     return 0;
}

int add_edge(unsigned v1, // voxel linear index.  
	     unsigned v2,
	     unsigned n_pts,
	     vnl_sparse_matrix<float> & G)
{
     unsigned n1 = 0, n2 = 0, k = 0;
     unsigned max_v = v1>v2?v1:v2;
     unsigned min_v = v1<v2?v1:v2;


     // when n1, n2 are the first voxel, the second voxel k should be smaller
     // than min_v
     for (k = 0; k < min_v; k ++) {
	  // compute the linear index of the connectivity variables.
	  n1 = v1 * n_pts + k;
	  n2 = v2 * n_pts + k;
	  G(n1, n2) = 1;
	  G(n2, n1) = 1;
	  // printf("    edge_edge(): n1 = %i, n2 = %i\n", n1, n2);
     }

     // when v1, v2 are the second voxel, the first voxel k should be larger
     // than max_v (since we're using the lower triangular matrix)
     for (k = max_v; k < n_pts; k ++) {
	  n1 = k * n_pts + v1;
	  n2 = k * n_pts + v2;
	  G(n1, n2) = 1;
	  G(n2, n1) = 1;
	  // printf("    edge_edge(): n1 = %i, n2 = %i\n", n1, n2);
     }

     return 0;
}

int E_step(const vnl_vector<double> & C,
	   vnl_matrix<double> & Z,
	   vnl_sparse_matrix<float> & G,
	   ParType par,
	   double beta)
{
     vnl_vector<double> exp_term (NCOMP, 0);     
     vnl_sparse_matrix<float>::row G_row;
     unsigned n1 = 0, k = 0, n2 = 0;
     unsigned v1 = 0, v2 = 0; // voxel linear index.
     for (v1 = 0; v1 < par.n_pts; v1 ++) {
	  // only update lower triangular, which means the vox2 should always be
	  // smaller than voxel 1.
	  for (v2 = 0; v2 < v1; v2 ++) {
	       // convert the linear voxel id to linear connnectivty variable
	       // id.
	       n1 = v1 * par.n_pts + v2;
	       G_row = G.get_row(n1);
	       for (k = 0; k < NCOMP; k ++) {
	       	    exp_term[k] = (-0.5) * log(2*PI)
	       		 - 0.5 * log(par.gmm.comp[k].var)
	       		 - pow(C[n1] - par.gmm.comp[k].mu, 2) / (2 * par.gmm.comp[k].var);
	       	    exp_term[k] += log(par.gmm.comp[k].pi);
	       }


	       for (k = 0; k < NCOMP; k++) {
		    vnl_vector<double> z(Z.cols(), 0);
		    z[k] = 1;
		    for (vnl_sparse_matrix<float>::row::const_iterator col_iter = G_row.begin(); col_iter != G_row.end();  ++ col_iter) {
			 n2 = (*col_iter).first;
			 // if z = 1
			 exp_term[1] += beta * dot_product(z, Z.get_row(n2));
		    } 
	       } // k

	       double m = exp_term.max_value();
	       exp_term = exp_term - m;

	       for (k = 0; k < NCOMP; k++) {
		    exp_term[k] = exp(exp_term[k]);
	       }
	       exp_term /= exp_term.sum();
	       Z.set_row(n1, exp_term);
	  } // v2

	  // deal with diagonal 
	  n1 = v1 * par.n_pts + v1;
	  Z.set_row(n1, vnl_vector<double>(Z.cols(), 0));
	  Z(n1, 2) = 1; // comp 2 is for perfect correlation. 
     } // v1

     return 0;
}

int M_step(const vnl_vector<double> & C,
	   const vnl_matrix<double> & Z,
	   ParType & par)
{
     // clear to zero before computation.
     for (unsigned k = 0; k < NCOMP; k ++) {
	  par.gmm.comp[k].mu = 0;
	  par.gmm.comp[k].numPts = 0;
	  par.gmm.comp[k].var = 0;
     }

     // estimate mu, first compute sum.
     unsigned v1 = 0, v2 = 0, n = 0, k = 0;
     for (v1 = 0; v1 < par.n_pts; v1 ++) {
	  // only update lower triangular.
	  for (v2 = 0; v2 < v1; v2 ++) {
	       n = v1 * par.n_pts + v2;	       
	       for ( k = 0; k < NCOMP; k ++) {
		    par.gmm.comp[k].numPts += Z(n, k);
		    par.gmm.comp[k].mu += Z(n, k) * C[n];
	       }
	  }// v2
     } // v1

     // compute mean from the sum
     for ( k = 0; k < NCOMP; k ++) {
	  par.gmm.comp[k].mu /= par.gmm.comp[k].numPts;
     }

     // compute variance.
     for (v1 = 0; v1 < par.n_pts; v1 ++) {
	  for (v2 = 0; v2 < v1; v2 ++) {
	       n = v1 * par.n_pts + v2;	       
	       for ( k = 0; k < NCOMP; k ++) {
		    par.gmm.comp[k].var += Z(n, k) * pow(C[n] - par.gmm.comp[k].mu, 2);
	       }
	  } // v2
     } // v1
     
     for ( k = 0; k < NCOMP; k ++) {
	  par.gmm.comp[k].var /= par.gmm.comp[k].numPts;
     }

     // estimate N
     double N = 0;
     for ( k = 0; k < NCOMP; k ++) {
	  N += par.gmm.comp[k].numPts;
     }

     // estimate pi
     for ( k = 0; k < NCOMP; k ++) {
	  par.gmm.comp[k].pi = par.gmm.comp[k].numPts / N;
     }


}
