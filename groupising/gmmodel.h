#include "common_gmmodel.h"

class GMModel
{
private:
     unsigned int m_numClusters;
     unsigned int m_numSubs;
     GMM m_gmm;
     double m_alpha;
     unsigned short m_verbose;


public:
     GMModel(VnlVectorImageType::Pointer obsPtr,
		 ImageType3DVecF::Pointer grpPtr,
		 ImageType3DChar::Pointer maskPtr,
		 unsigned numClusters,
		 unsigned short verbose );

     ~GMModel(); 

     int printSelf(std::string format);

     int Estep(ImageType3DVecF::Pointer grpPtr,
		   VnlVectorImageType::Pointer obsPtr,
		   ImageType3DChar::Pointer maskPtr);

     int EstimateMu(ImageType3DVecF::Pointer grpPtr,
		    VnlVectorImageType::Pointer obsPtr,
		    ImageType3DChar::Pointer maskPtr);

     int estimateSigma(ImageType3DVecF::Pointer grpPtr,
			   VnlVectorImageType::Pointer obsPtr,
			   ImageType3DChar::Pointer maskPtr);

     double llhood(ImageType3DVecF::Pointer grpPtr,
		   VnlVectorImageType::Pointer obsPtr,
		   ImageType3DChar::Pointer maskPtr);

};


