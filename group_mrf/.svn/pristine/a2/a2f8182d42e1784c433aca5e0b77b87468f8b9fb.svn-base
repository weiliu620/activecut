#include "common_vmmodel.h"
#include "vmmodel.h"
#include <boost/math/distributions/students_t.hpp>

extern twister_base_gen_type mygenerator; 

VMModel::VMModel(std::vector<VnlVectorImageType::Pointer> & fmriVec,
		 std::vector<ImageType3DChar::Pointer> & grpVec,
		 ImageType3DChar::Pointer maskPtr,
		 unsigned numClusters,
		 float beta_g,
		 float gamma,
		 unsigned short verbose )
{
     unsigned subIdx = 0, clsIdx = 0;

     // init members.
     m_verbose = verbose;
     m_beta_g = beta_g;


     VnlVectorImageType::SizeType fmriSize = 
	  fmriVec[0]->GetLargestPossibleRegion().GetSize();

     m_numSubs = fmriVec.size();
     m_gamma = gamma; // internal gamma
     m_numClusters = numClusters;
     m_numSamples = grpVec.size();

     VnlVectorImageType::IndexType fmriIdx;
     fmriIdx.Fill(0);
     m_tsLength = fmriVec[0]->GetPixel(fmriIdx).size();

     // init m_vmm
     m_vmm.comp.resize(m_numClusters);
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_vmm.comp[clsIdx].sub.resize(m_numSubs);
	  for (subIdx = 0; subIdx < m_numSubs; subIdx++) {
	       m_vmm.comp[clsIdx].sub[subIdx].mu.set_size(m_tsLength);
	       m_vmm.comp[clsIdx].sub[subIdx].mu.fill(0);
	       m_vmm.comp[clsIdx].sub[subIdx].meanNorm = 0;
	       m_vmm.comp[clsIdx].sub[subIdx].kappa = 0;
	       m_vmm.comp[clsIdx].sub[subIdx].numPoints = 0;
	  }
     }

     // initial total number of data points (voxels) in the mask.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());
     ImageType3DChar::IndexType maskIdx;     
     m_vmm.totalPts = 0;
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       m_vmm.totalPts ++;
	  }  // maskIt > 0
     } // maskIt
     if(m_verbose >= 1) {
	  printf("total number of voxels within mask: %i.\n", m_vmm.totalPts);
     }
}

VMModel::~VMModel()
{

}

int VMModel::printSelf(std::string format)
{
     unsigned int clsIdx = 0, subIdx = 0;

     if (format.compare("normal") == 0) {

	  for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
	       printf("comp[%i]: prop = %f\n", clsIdx + 1, m_vmm.comp[clsIdx]. prop);
	       for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
		    printf("     sub[%i]: # Pts: %4.2f, meanNorm = %f, kappa = %4.2f,  ",
				subIdx + 1,
				m_vmm.comp[clsIdx].sub[subIdx].numPoints,
				m_vmm.comp[clsIdx].sub[subIdx].meanNorm,
				m_vmm.comp[clsIdx].sub[subIdx].kappa);
		    printf("mu = ");
		    printVnlVector(m_vmm.comp[clsIdx].sub[subIdx].mu, 3);
	       } // subIdx
	  } // clsIdx.
     }
     else if (format.compare("table") == 0) {
     }
     else {
	  printf("VMModel::print(): print format not recognized.\n");
     }
     fflush(stdout);
     
     return 0;
}



int VMModel::mcsampling( std::vector<ImageType3DChar::Pointer> & grpVec,
			 std::vector<VnlVectorImageType::Pointer> & fmriVec,
			 ImageType3DChar::Pointer maskPtr,
			 unsigned burnin)
{
     unsigned subIdx = 0;
     unsigned scanIdx = 0;
     for (scanIdx = 0; scanIdx < burnin + m_numSamples; scanIdx ++) {
	  GrpSampling(grpVec[m_numSamples -1], fmriVec, maskPtr);

	  if (m_verbose >= 1) {
	       printf("      mcsampling(): scan %i done.\n", scanIdx);
	  }

	  // After burnin period, it's time to save it to correct place.
	  if (scanIdx >= burnin) {
	       IteratorType3DChar srcIt(grpVec[m_numSamples - 1], 
					grpVec[m_numSamples - 1]->GetLargestPossibleRegion() );

	       IteratorType3DChar destIt(grpVec[scanIdx - burnin], 
					 grpVec[scanIdx - burnin]->GetLargestPossibleRegion() );

	       for (srcIt.GoToBegin(), destIt.GoToBegin(); !srcIt.IsAtEnd(); ++ srcIt, ++ destIt) {
		    destIt.Set(srcIt.Get());
	       }
	  } // scanIdx >= burnin
     } // scanIdx
}


int VMModel::GrpSampling(ImageType3DChar::Pointer grpPtr,
			 std::vector<VnlVectorImageType::Pointer> & fmriVec,
			 ImageType3DChar::Pointer maskPtr)
{
     double p_acpt = 0;
     double denergy= 0, dllenergy =0;
     int cand = 1;
     int currentLabel = 0;
     
     unsigned clsIdx = 0, subIdx = 0;
     // Uniform integer generator.
     boost::uniform_int<> uni_int(0, m_numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     //mask
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType maskIdx;     


     // images
     VnlVectorImageType::SizeType fmriSize = 
	  fmriVec[0]->GetLargestPossibleRegion().GetSize();
     VnlVectorImageType::IndexType fmriIdx;    
     fmriIdx.Fill(0);


     // Define neighborhood iterator
     typedef itk::NeighborhoodIterator< ImageType3DChar> NeighborhoodIteratorType;
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType neiGrpIt(radius, grpPtr, grpPtr->GetLargestPossibleRegion() );

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1}};

     ImageType3DChar::SizeType grpSize = grpPtr->GetLargestPossibleRegion().GetSize();
     vnl_vector <float> timeSeries(m_tsLength, 0);

     // compute normalization constant of von Mises
     // Fisher. vmfLogConst[i] is the log (c_d (kappa) ) for the i'th
     // clusters. See page 1350 of "Clustering on the Unit-Sphere
     // using von Mises Fisher..."
     vnl_matrix<double> vmfLogConst(m_numClusters, m_numSubs);
     float myD = m_tsLength;
     double const Pi = 4 * atan(1);
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       vmfLogConst[clsIdx][subIdx] = (myD/2 - 1) * log (m_vmm.comp[clsIdx].sub[subIdx].kappa) - myD/2 * log(2*Pi) -  logBesselI(myD/2 - 1, m_vmm.comp[clsIdx].sub[subIdx].kappa);
	  }
     }

     for (neiGrpIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ neiGrpIt, ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       currentLabel = neiGrpIt.GetCenterPixel();
	       cand = roll_die();
	       denergy = 

		     int(cand != neiGrpIt.GetPixel(xminus))
		    - int(currentLabel != neiGrpIt.GetPixel(xminus))

		    + int(cand != neiGrpIt.GetPixel(xplus)) 
		    - int(currentLabel != neiGrpIt.GetPixel(xplus))

		    + int(cand != neiGrpIt.GetPixel(yminus)) 
		    - int(currentLabel != neiGrpIt.GetPixel(yminus))

		    + int(cand != neiGrpIt.GetPixel(yplus)) 
		    - int(currentLabel != neiGrpIt.GetPixel(yplus))

		    + int(cand != neiGrpIt.GetPixel(zminus)) 
		    - int(currentLabel != neiGrpIt.GetPixel(zminus))

		    + int(cand != neiGrpIt.GetPixel(zplus)) 
		    - int(currentLabel != neiGrpIt.GetPixel(zplus));	       

	       denergy = denergy * m_beta_g;


	       // data energy term.

	       // change of likelihood energy.

	       maskIdx = maskIt.GetIndex();
	       dllenergy = 0;
// #pragma omp parallel for private(subIdx) private(timeSeries)  reduction(+:dllenergy)
	       for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {

		    timeSeries = fmriVec[subIdx]->GetPixel(maskIdx);
		    dllenergy += ( - vmfLogConst[cand][subIdx] - m_vmm.comp[cand].sub[subIdx].kappa * inner_product(timeSeries, m_vmm.comp[cand].sub[subIdx].mu) +  vmfLogConst[currentLabel][subIdx] + m_vmm.comp[currentLabel].sub[subIdx].kappa * inner_product(timeSeries, m_vmm.comp[currentLabel].sub[subIdx].mu) );
	       }

	       // normalize likelihood so it does not swamp prior because of
	       // high dimension and multi subjects.
	       dllenergy =  m_gamma * dllenergy ;


	       denergy = denergy + dllenergy;

	       // temperature
	       denergy = denergy / m_temperature;
	       
	       if (denergy <= 0) {
		    neiGrpIt.SetCenterPixel( cand );
	       }
	       else {
		    p_acpt = exp(-denergy);
		    if (uni() < p_acpt) {
			 neiGrpIt.SetCenterPixel( cand );
		    }
	       }
	  } // maskIt > 0
     } // neiGrpIt
     return 0;
}


// Given group's labels, estimate mu and numPoints.
int VMModel::EstimateMu( std::vector<ImageType3DChar::Pointer> & grpVec,
			 std::vector<VnlVectorImageType::Pointer> & fmriVec,
			 ImageType3DChar::Pointer maskPtr)
{
     unsigned int clsIdx = 0, subIdx = 0;

     // reset all mu to zero.
     for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	  for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	       m_vmm.comp[clsIdx].sub[subIdx].numPoints = 0;
	       m_vmm.comp[clsIdx].sub[subIdx].mu = 0;
	  }
     }

     // compute mu and numPoints for all sub, all clusters.
#pragma omp parallel for
     for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	  estimateMuSub(grpVec, fmriVec[subIdx], maskPtr, subIdx);
	  if (m_verbose >= 2) {
	       printf("estimateMu(): estimabeMuSub() done for subject %i.\n", subIdx);
	  }
     } // for subIdx

     // compute proportion for clusters from sub's numPoints. Since all subject
     // share one label map z, they have same numPoints for one clusters. So
     // here I just sum over all subjects' numPoints and divide by num of subs,
     // which is same to single sub's numPoints. Then I divide it by total
     // number of points.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_vmm.comp[clsIdx].prop = 0;
	  for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       m_vmm.comp[clsIdx].prop += float (m_vmm.comp[clsIdx].sub[subIdx].numPoints) / float(m_vmm.totalPts * m_numSubs);
	  }
     }

     return 0;
}


int VMModel::estimateMuSub(std::vector<ImageType3DChar::Pointer> & grpVec,
			   VnlVectorImageType::Pointer fmriPtr,
			   ImageType3DChar::Pointer maskPtr,
			   unsigned subIdx)
{
     unsigned int clsIdx = 0;

     // images
     VnlVectorImageType::IndexType imageIdx;    
     VnlVectorImageType::SizeType imageSize = fmriPtr->GetLargestPossibleRegion().GetSize();
     IteratorTypeVnlVector fmriIt(fmriPtr, 
				   fmriPtr->GetLargestPossibleRegion().GetSize() );

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     ImageType3DChar::RegionType maskRegion = maskPtr->GetLargestPossibleRegion();
     IteratorType3DChar maskIt(maskPtr, maskRegion);
     ImageType3DChar::SizeType maskSize = maskRegion.GetSize();


     vnl_vector <float> timeSeries(m_tsLength, 0);

     // reset all mu to zero.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_vmm.comp[clsIdx].sub[subIdx].numPoints = 0;
	  m_vmm.comp[clsIdx].sub[subIdx].mu = 0;
     }


     imageIdx[4] = subIdx;
     for (unsigned short sampleIdx = 0; sampleIdx < m_numSamples; sampleIdx ++) {
	  ConstIteratorType3DChar sampleIt(grpVec[sampleIdx], 
					   grpVec[sampleIdx]->GetLargestPossibleRegion() );

	  for (fmriIt.GoToBegin(), maskIt.GoToBegin(), sampleIt.GoToBegin(); 
	       !maskIt.IsAtEnd(); 
	       ++ fmriIt, ++ maskIt, ++ sampleIt) {

	       if (maskIt.Get() > 0) {
		    clsIdx = sampleIt.Get();
		    m_vmm.comp[clsIdx].sub[subIdx].mu += fmriIt.Get();
		    m_vmm.comp[clsIdx].sub[subIdx].numPoints ++;
	       }  // maskIt > 0
	  } // for maskIt

     } // sampleIdx

     // compute mean time series and meanNorm.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  // meanNorm see page 1350 of "Clustering on the Unit
	  // Hypersphere using von Mises-Fisher distributions" by
	  // Banerjee.
	  m_vmm.comp[clsIdx].sub[subIdx].meanNorm = m_vmm.comp[clsIdx].sub[subIdx].mu.two_norm() / m_vmm.comp[clsIdx].sub[subIdx].numPoints;
	  m_vmm.comp[clsIdx].sub[subIdx].mu.normalize();
     }
     return 0;
}


int VMModel::estimateKappa(float timeSeriesLength)
{
     unsigned clsIdx = 0, subIdx = 0;

     float RBar = 0, kappa_mle = 0;
     float Dim = timeSeriesLength;
     unsigned N = 0; // num of points for this cluster.

     for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	  if (m_verbose >= 1) {
	       printf("estimateKappaWithPrior(), subIdx[%u]. \n", subIdx);
	  }
	  for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	       N = m_vmm.comp[clsIdx].sub[subIdx].numPoints;
	       RBar = m_vmm.comp[clsIdx].sub[subIdx].meanNorm;
	       kappa_mle = RBar * (Dim - RBar * RBar) / (1 - RBar * RBar);
	       m_vmm.comp[clsIdx].sub[subIdx].kappa = kappa_mle;
	  } // clsIdx
     } // subIdx
     return 0;
}


int VMModel::SetTemperature(float newTemp)
{
     m_temperature = newTemp;
     return 0;
}


int VMModel::ICM( std::vector<ImageType3DChar::Pointer> & grpVec,
		  std::vector<VnlVectorImageType::Pointer> & fmriVec,
		  ImageType3DChar::Pointer maskPtr,
		  unsigned burnin)
{
     unsigned subIdx = 0;
     unsigned scanIdx = 0;
     for (scanIdx = 0; scanIdx < burnin; scanIdx ++) {
	  ICMInternal(grpVec[m_numSamples -1], fmriVec, maskPtr);
	  if (m_verbose >= 1) {
	       printf("      ICM(): scan %i done.\n", scanIdx);
	  }

     } // scanIdx
}

int VMModel::ICMInternal(ImageType3DChar::Pointer grpPtr,
		 std::vector<VnlVectorImageType::Pointer> & fmriVec,
		 ImageType3DChar::Pointer maskPtr)
{
     vnl_vector<double> denergy(m_numClusters, 0);
     double dllenergy = 0;
     unsigned clsIdx = 0, subIdx = 0;
     int currentLabel = 0, cand = 0;

     //mask
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType maskIdx;     


     // images
     VnlVectorImageType::SizeType fmriSize = 
	  fmriVec[0]->GetLargestPossibleRegion().GetSize();
     VnlVectorImageType::IndexType fmriIdx;    
     fmriIdx.Fill(0);


     // Define neighborhood iterator
     typedef itk::NeighborhoodIterator< ImageType3DChar> NeighborhoodIteratorType;
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType neiGrpIt(radius, grpPtr, grpPtr->GetLargestPossibleRegion() );

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1}};

     ImageType3DChar::SizeType grpSize = grpPtr->GetLargestPossibleRegion().GetSize();
     vnl_vector <float> timeSeries(m_tsLength, 0);

     // compute normalization constant of von Mises
     // Fisher. vmfLogConst[i] is the log (c_d (kappa) ) for the i'th
     // clusters. See page 1350 of "Clustering on the Unit-Sphere
     // using von Mises Fisher..."
     vnl_matrix<double> vmfLogConst(m_numClusters, m_numSubs);
     float myD = m_tsLength;
     double const Pi = 4 * atan(1);
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       vmfLogConst[clsIdx][subIdx] = (myD/2 - 1) * log (m_vmm.comp[clsIdx].sub[subIdx].kappa) - myD/2 * log(2*Pi) -  logBesselI(myD/2 - 1, m_vmm.comp[clsIdx].sub[subIdx].kappa);
	  }
     }

     
     for (neiGrpIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ neiGrpIt, ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       currentLabel = neiGrpIt.GetCenterPixel();

	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
		    denergy[clsIdx] = 
			 int(clsIdx != neiGrpIt.GetPixel(xminus))
			 - int(currentLabel != neiGrpIt.GetPixel(xminus))

			 + int(clsIdx != neiGrpIt.GetPixel(xplus)) 
			 - int(currentLabel != neiGrpIt.GetPixel(xplus))

			 + int(clsIdx != neiGrpIt.GetPixel(yminus)) 
			 - int(currentLabel != neiGrpIt.GetPixel(yminus))

			 + int(clsIdx != neiGrpIt.GetPixel(yplus)) 
			 - int(currentLabel != neiGrpIt.GetPixel(yplus))

			 + int(clsIdx != neiGrpIt.GetPixel(zminus)) 
			 - int(currentLabel != neiGrpIt.GetPixel(zminus))

			 + int(clsIdx != neiGrpIt.GetPixel(zplus)) 
			 - int(currentLabel != neiGrpIt.GetPixel(zplus));	       

		    denergy[clsIdx] = denergy[clsIdx] * m_beta_g;


		    // data energy term.

		    // change of likelihood energy.

		    maskIdx = maskIt.GetIndex();
		    dllenergy = 0;
// #pragma omp parallel for private(subIdx) private(timeSeries)  reduction(+:dllenergy)
		    for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {

			 timeSeries = fmriVec[subIdx]->GetPixel(maskIdx);
			 dllenergy += ( - vmfLogConst[clsIdx][subIdx] - m_vmm.comp[clsIdx].sub[subIdx].kappa * inner_product(timeSeries, m_vmm.comp[clsIdx].sub[subIdx].mu) +  vmfLogConst[currentLabel][subIdx] + m_vmm.comp[currentLabel].sub[subIdx].kappa * inner_product(timeSeries, m_vmm.comp[currentLabel].sub[subIdx].mu) );
		    }

		    // normalize likelihood so it does not swamp prior because of
		    // high dimension and multi subjects.
		    dllenergy = dllenergy / m_gamma;


		    denergy[clsIdx] = denergy[clsIdx] + dllenergy;
	       
	       } // clsIdx
	       neiGrpIt.SetCenterPixel( denergy.arg_min() );	       
	  } // maskIt > 0
     } // neiGrpIt
     return 0;
}
