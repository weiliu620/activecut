int logistic(vnl_vector<double> & score_map,
	     const vnl_matrix<double> & data,
	     vnl_sparse_matrix <double> & con_map, // can not set const due to a vnl bug.
	     const vnl_vector<unsigned> & alpha,
	     const ParType & par);
