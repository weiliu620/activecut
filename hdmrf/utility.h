int save_corr(const vnl_vector<double> & corr,
	      std::string out_file,
	      ParType par);

// save hidden label map.
int save_labelmap(const vnl_matrix<double> & Z,
		  std::string out_file,
		  ParType par);
// print each element of the sparse matrix.
int print_sparse_matrix(vnl_sparse_matrix<float> mat);

int print_par(ParType Par);


