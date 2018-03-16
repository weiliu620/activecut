double logBesselI(float nu, double x);

int PrintPar(unsigned prlevel, ParStruct & par);

int save3dcharInc(ImageType3DChar::Pointer ip, std::string fname);
int save3dchar(ImageType3DChar::Pointer ip, std::string fname);
int printVnlVector(vnl_vector<float> vec, unsigned numElements);
int printVnlVector(vnl_vector<double> vec, unsigned numElements);

