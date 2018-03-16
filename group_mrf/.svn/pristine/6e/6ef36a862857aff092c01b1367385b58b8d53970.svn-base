#include "commonalt.h"



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
     double m_gamma;
     float m_temperature;
     unsigned int m_numSamples;
     unsigned int m_tsLength;
     unsigned short m_verbose;

public:
     MCModel(std::vector<VnlVectorImageType::Pointer> & fmriVec,
	     std::vector< std::vector<ImageType3DChar::Pointer> > &  sampleVec,
	     std::vector<ImageType3DChar::Pointer> &  grpVec,
	     ImageType3DChar::Pointer maskPtr,
	     ParStruct par);

     ~MCModel(); 

     int printSelf(std::string format);

     int printMeanKappa();

     int estimateMu(ImageType3DChar::Pointer groupPtr,
		    std::vector<VnlVectorImageType::Pointer> & fmriVec,
		    ImageType3DChar::Pointer maskPtr);

     int estimateMuSub(ImageType3DChar::Pointer groupPtr,
		       VnlVectorImageType::Pointer fmriPtr,
		       ImageType3DChar::Pointer maskPtr,
		       unsigned subIdx);

     int estimateKappa();

     int SamplingSub(ImageType3DChar::Pointer groupPtr,
		     VnlVectorImageType::Pointer fmriPtr,
		     ImageType3DChar::Pointer  samplePtr,
		     unsigned subIdx);


     int mcsampling(std::vector<ImageType3DChar::Pointer> & grpVec,
		    std::vector<VnlVectorImageType::Pointer> & fmriVec,
		    ImageType3DChar::Pointer maskPtr,
		    std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		    unsigned burnin);

     int GrpSampling(ImageType3DChar::Pointer groupPtr,
		     std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		     ImageType3DChar::Pointer maskPtr);

     int estimateMuSubSample(std::vector<ImageType3DChar::Pointer> & sampleSubVec,
			     VnlVectorImageType::Pointer imagePtr,
			     ImageType3DChar::Pointer maskPtr,
			     unsigned subIdx);

     int EstimateMuFromSample(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
			      std::vector<VnlVectorImageType::Pointer> & fmriVec,
			      ImageType3DChar::Pointer maskPtr);


     float GetTemperature();
     int SetTemperature(float newTemp);

     double ComputeQ3(VnlVectorImageType::Pointer fmriPtr,
		      ImageType3DChar::Pointer samplePtr,
		      ImageType3DChar::Pointer maskPtr,
		      unsigned subIdx);

     double ComputeQ1(ImageType3DChar::Pointer &  grpPtr,
		      ImageType3DVecUS::Pointer & histPtr,
		      ImageType3DChar::Pointer maskPtr);

     double ComputeQ2(ImageType3DChar::Pointer &  grpPtr,
		      ImageType3DChar::Pointer & samplePtr,
		      ImageType3DChar::Pointer maskPtr,
		      unsigned subIdx);

     double ComputeEnergy(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
			  std::vector<ImageType3DChar::Pointer> &  grpVec,
			  ImageType3DChar::Pointer maskPtr,
			  std::vector<VnlVectorImageType::Pointer> & fmriVec);

     // Estimate alpha and beta, the parameters in prior P(G,Z).
     double EstimatePriorPar(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
			     std::vector<ImageType3DChar::Pointer> &  grpVec,
			     ImageType3DChar::Pointer maskPtr,
			     float alpha0,
			     float sigma2);

     // Compute 1st and 2nd derivative of E(G) on beta_g.
     double GrpBetaDrv(std::vector<ImageType3DChar::Pointer> &  grpVec,
		       std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		       std::vector <ImageType3DVecUS::Pointer> & histVec,
		       ImageType3DChar::Pointer maskPtr,
		       double & drv1, double & drv2);
     // Compute 1st and 2nd derivatives of E(Z|G) on beta_z.
     double SubBetaDrv(std::vector<ImageType3DChar::Pointer> &  grpVec,
		       std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		       ImageType3DChar::Pointer maskPtr,
		       double & drv1, double & drv2);

    // Compute 1st and 2nd derivatives of Q1 (-log P(G|Z) ) with beta, for each
    // subject.
     double SubBetaDrvPerSub(std::vector<ImageType3DChar::Pointer> &  grpVec,
				      std::vector<ImageType3DChar::Pointer> & sampleSubVec,
				      ImageType3DChar::Pointer maskPtr,
				      double & drv1, double & drv2);

     // Compute 1st and 2nd derivatives of E(G) on alpha.
     double GrpAlphaDrv(std::vector<ImageType3DChar::Pointer> &  grpVec,
			std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
			std::vector <ImageType3DVecUS::Pointer> & histVec,
			ImageType3DChar::Pointer maskPtr,
			double & drv1, double & drv2);
     // Compute 1st and 2nd derivatives of E(Z|G) on alpha.
     double SubAlphaDrv(std::vector<ImageType3DChar::Pointer> &  grpVec,
			std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
			ImageType3DChar::Pointer maskPtr,
			double & drv1, double & drv2);

    // Compute 1st and 2nd derivatives of Q1 (-log P(G|Z) ) with alpha, for each
    // subject.
     double SubAlphaDrvPerSub(std::vector<ImageType3DChar::Pointer> &  grpVec,
			      std::vector<ImageType3DChar::Pointer>  & sampleSubVec,
			      ImageType3DChar::Pointer maskPtr,
			      double & drv1, double & drv2);

     double PlotQvsBeta(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
			std::vector<ImageType3DChar::Pointer> &  grpVec,
			ImageType3DChar::Pointer maskPtr,
			std::vector<VnlVectorImageType::Pointer> & fmriVec);

};


