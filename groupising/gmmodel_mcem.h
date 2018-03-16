#include <common_mcem.h>
#include<gmmodel_graph.h>

class GMModel
{
private:
     unsigned int m_numClusters;
     unsigned int m_numSubs;
     GMM m_gmm;
     double m_alpha;
     unsigned short m_verbose;
     double m_beta_g;
     double m_beta_z;
     double m_gamma;
     double m_temperature;
     unsigned m_numSamples;
     GCoptimizationGeneralGraph * m_gc;
     BiMapType m_biMap;

public:
     GMModel(std::vector<ImageType3DFloat::Pointer> & obsVec,
	     std::vector<ImageType3DChar::Pointer> &  grpVec,
	     std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
	     ImageType3DChar::Pointer maskPtr,
	     ParStruct par);

     ~GMModel(); 

     // print model parameter, beta, alpha, gamma, mu, sigma, etc. Print human
     // readable format if "format" set to "normal", and print table format for
     // ploting if "format" is set to "table".
     int printSelf(std::string format);

     // Sampling for both group and subject label map. The sweep repeat burnin +
     // numSamples times.
     int Sampling(std::vector<ImageType3DChar::Pointer> & grpPtr,
		  std::vector<ImageType3DFloat::Pointer> & obsVec,
		  ImageType3DChar::Pointer maskPtr,
		  std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		  unsigned burnin);

     // Sampling for a single subject. 
     int SamplingSub(ImageType3DChar::Pointer  grpPtr,
		     ImageType3DFloat::Pointer obsPtr,
		     ImageType3DChar::Pointer maskPtr,
		     ImageType3DChar::Pointer sampleVec,
		     unsigned subIdx);


     // Sampling of group.
     int SamplingGrp(ImageType3DChar::Pointer grpPtr,
		     ImageType3DChar::Pointer maskPtr,
		     std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec);

     // Estimation of likelihood parameters mu and sigma (Gaussian
     // distribution).
     int estimatePar(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
			 std::vector<ImageType3DFloat::Pointer> & obsVec,
			 ImageType3DChar::Pointer maskPtr);

     // Estimate mu for a single subject.
     int estimateMuSub(std::vector<ImageType3DChar::Pointer> & sampleSubVec,
		       ImageType3DFloat::Pointer  obsVec,
		       ImageType3DChar::Pointer maskPtr,
		       unsigned subIdx);
     // Estimate sigma for a single subject.
     int estimateSigmaSub(std::vector<ImageType3DChar::Pointer> & sampleSubVec,
			  ImageType3DFloat::Pointer  obsVec,
			  ImageType3DChar::Pointer maskPtr,
			  unsigned subIdx);

     // obsolete.
     int estimateParSingleSample(
	  std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
	  std::vector<ImageType3DFloat::Pointer>  & obsVec,
	  ImageType3DChar::Pointer maskPtr);

     // Estimate alpha, beta_g, beta_z in the MRF prior.
     double EstimatePriorPar(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
			  std::vector<ImageType3DChar::Pointer> &  grpVec,
			  ImageType3DChar::Pointer maskPtr,
			  std::string position); // which location to run this fun.);

     // Compute first and second derivative of log P(G). (G is the group label
     // map).
     double GrpBetaDrv(std::vector<ImageType3DChar::Pointer> &  grpVec,
			     std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
			     std::vector <ImageType3DVecUS::Pointer> & histVec,
			     ImageType3DChar::Pointer maskPtr,
			     double & drv1, double & drv2);
 
     // Compute 1st and 2nd derivatives for subject label map. All subjects
     // share same beta_z.
     double SubBetaDrv(std::vector<ImageType3DChar::Pointer> &  grpVec,
		    std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		    ImageType3DChar::Pointer maskPtr,
		    double & drv1, double & drv2);

     // Compute 1st and 2nd deravatives for alpha (the connection parameter
     // between group and subject labels.
     double GrpAlphaDrv(std::vector<ImageType3DChar::Pointer> &  grpVec,
		     std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		     std::vector <ImageType3DVecUS::Pointer> & histVec,
		     ImageType3DChar::Pointer maskPtr,
		     double & drv1, double & drv2);


     double SubAlphaDrv(std::vector<ImageType3DChar::Pointer> &  grpVec,
		     std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		     ImageType3DChar::Pointer maskPtr,
		     double & drv1, double & drv2);


     int InitGraph(std::vector<ImageType3DFloat::Pointer> & obsVec,
			std::vector<ImageType3DChar::Pointer> &  grpVec,
			std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
			ImageType3DChar::Pointer maskPtr);

     int graphcuts(std::vector<ImageType3DFloat::Pointer> & obsVec,
		   std::vector<ImageType3DChar::Pointer> &  grpVec,
		   std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		   ImageType3DChar::Pointer maskPtr);

     // Compute -log P(X|Z) for all subjects.
     double Q3(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
	      std::vector<ImageType3DFloat::Pointer> & obsVec,
	      ImageType3DChar::Pointer maskPtr);

     // Iterated Conditional Mode (ICM) approximation on single subject.
     float ICMSub(ImageType3DChar::Pointer grpPtr,
		  ImageType3DFloat::Pointer obsPtr,
		  ImageType3DChar::Pointer maskPtr,
		  ImageType3DChar::Pointer samplePtr,
		  unsigned subIdx);

     // Iterated Conditional Mode (ICM) approximation on group labels.
     int ICMGrp(ImageType3DChar::Pointer grpPtr,
		ImageType3DChar::Pointer maskPtr,
		std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec);

     // ICM wrapper, to estimation both group and subjects label map.
     int ICM(std::vector<ImageType3DChar::Pointer> &  grpVec,
	     std::vector<ImageType3DFloat::Pointer> & obsVec,
	     ImageType3DChar::Pointer maskPtr,
	     std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
	     float tol); // max proportion of voxels changed for convergence.

     float GetTemperature();

     int SetTemperature(float newTemp);
};


