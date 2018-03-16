int saveimage2d(ImageType2D::Pointer ip, std::string fname);

int saveimage3d(ImageType3D::Pointer ip, std::string fname);

int saveimage4dchar(ImageType4DChar::Pointer ip, std::string fname);

int save3dcharInc(ImageType3DChar::Pointer ip, std::string fname);

int saveExtractSlice(ImageType3D::Pointer inPtr,
		     unsigned sliceIdx, std::string fname);

int saveExtractVolume(ImageType4DChar::Pointer inPtr,
		      unsigned volumeIdx, std::string fname);

int printVnlVector(vnl_vector<float> vec, unsigned numElements);

int printVnlMatrix(vnl_matrix<float> mat, unsigned numElements);

double logBesselI(float nu, double x);

int printvmm(std::vector <VMM> & vmm);

int SaveSample(ImageType3DChar::Pointer samplePtr, 		 
	       std::string filename);

// used by segimagealt. for allocating memory, initialize, and also used as mask.
int InitSamples(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec, 
		ImageType3DChar::Pointer grpPtr);

int SaveSamples(std::vector< std::vector<ImageType3DChar::Pointer> > &  sampleVec,
		ImageType3DChar::Pointer maskPtr,
		std::string basename);

int SaveGrpSamples( std::vector<ImageType3DChar::Pointer> &  sampleVec,
		    ImageType3DChar::Pointer maskPtr,
		    std::string filename);

int saveimage(ImageType3DS::Pointer ip, std::string fname);
int just_test(int a);
