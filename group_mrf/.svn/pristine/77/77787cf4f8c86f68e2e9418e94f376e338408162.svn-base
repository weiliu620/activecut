#include "commonalt.h"

#ifndef __COMMONALT_H__
#define __COMMONALT_H__

class MCModel
{
private:
     unsigned int m_numClusters;
     unsigned int m_numSubs;
     /* std::vector < CompType >  m_cls,  m_oldCls; */
     std::vector <VMM> m_vmm;
     double m_alpha;
     double m_beta_g;
     double m_beta_z;
     float m_temperature;
     unsigned short m_verbose;

     // neighborSum pointer, for log likelihood evaluation and also
     // for estimating beta.
     ImageType6DChar::Pointer m_neighborSumPtr;

public:
     MCModel(ImageType5DFloat::Pointer imagePtr,
	     ImageType3DChar::Pointer groupPtr,
	     unsigned numClusters,
	     unsigned short verbose);


     ~MCModel(); 

     int printSelf(std::string format);

     int testClusterCenterMove(float eps);

     // Init cluster center by K-means.
     double kmeans(ImageType5DFloat::Pointer imagePtr, 
		   ImageType3DChar::Pointer groupPtr,
		   ImageType3DChar::Pointer maskPtr,
		   float converge_eps);

     // init kmeans by kmeans++ method.
     int kmeanspp(ImageType5DFloat::Pointer imagePtr,
			   ImageType3DChar::Pointer maskPtr);

     int printMeanKappa();

     // compute kmeanspp's pdf distance for single subject. (can be
     // used in openmp for parallelization).
     int minDistPerSub(ImageType3DFloat::Pointer pdfPtr, 
		       ImageType5DFloat::Pointer imagePtr,
		       ImageType3DChar::Pointer maskPtr,
		       unsigned int clsIdx,
		       unsigned int subIdx);
     
     double estLabelsByMean(ImageType5DFloat::Pointer imagePtr,
			    ImageType3DChar::Pointer maskPtr,
			    ImageType3DChar::Pointer groupPtr);

     double estLabelsByMeanSub(ImageType5DFloat::Pointer imagePtr,
			       ImageType3DChar::Pointer maskPtr,
			       ImageType4DFloat::Pointer allSubDistPtr,
			       unsigned int subIdx);

     int estimateMu(ImageType3DChar::Pointer groupPtr,
		    ImageType5DFloat::Pointer imagePtr,
		    ImageType3DChar::Pointer maskPtr);


     int estimateMuSub(ImageType3DChar::Pointer groupPtr,
		       ImageType5DFloat::Pointer imagePtr,
		       ImageType3DChar::Pointer maskPtr,
		       unsigned subIdx);

     int KMeansConverge(std::vector<VMM> & curvmm, std::vector<VMM> & prevmm);
     int SetVMM(std::vector <VMM> & vmm);
     int GetVMM(std::vector <VMM> & vmm);
};





#endif
