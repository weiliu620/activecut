#include "commonalt.h"
typedef itk::Image< vnl_vector< float >, 3 > VnlVectorImageType;

#ifndef __COMMONALT_H__
#define __COMMONALT_H__

class VMModel
{
private:
     unsigned int m_numClusters;
     unsigned int m_numSubs;

     VMM m_vmm;
     double m_beta_g;
     float m_gamma;
     unsigned short m_verbose;
     unsigned m_numSamples;
     unsigned m_tsLength;
     float m_temperature;

public:
     VMModel(std::vector<VnlVectorImageType::Pointer> & fmriVec,
	     std::vector<ImageType3DChar::Pointer> & grpVec,
	     ImageType3DChar::Pointer maskPtr,
	     unsigned numClusters,
	     float beta_g,
	     float gamma,
	     unsigned short verbose );


     ~VMModel(); 

     int printSelf(std::string format);

     int printMeanKappa();


     int estimateKappa(float timeSeriesLength);


     int mcsampling( std::vector<ImageType3DChar::Pointer> & grpVec,
			 std::vector<VnlVectorImageType::Pointer> & fmriVec,
			 ImageType3DChar::Pointer maskPtr,
			 unsigned burnin);

     int GrpSampling(ImageType3DChar::Pointer grpPtr,
			 std::vector<VnlVectorImageType::Pointer> & fmriVec,
			 ImageType3DChar::Pointer maskPtr);

     int EstimateMu( std::vector<ImageType3DChar::Pointer> & grpVec,
			 std::vector<VnlVectorImageType::Pointer> & fmriVec,
			 ImageType3DChar::Pointer maskPtr);

     int estimateMuSub(std::vector<ImageType3DChar::Pointer> & grpVec,
		       VnlVectorImageType::Pointer fmriPtr,
		       ImageType3DChar::Pointer maskPtr,
		       unsigned subIdx);

     int SetTemperature(float newTemp);

     int ICM( std::vector<ImageType3DChar::Pointer> & grpVec,
	      std::vector<VnlVectorImageType::Pointer> & fmriVec,
	      ImageType3DChar::Pointer maskPtr,
	      unsigned burnin);

     int ICMInternal(ImageType3DChar::Pointer grpPtr,
		     std::vector<VnlVectorImageType::Pointer> & fmriVec,
		     ImageType3DChar::Pointer maskPtr);
};



#endif
