#include "commonalt.h"
#include "MCModel.h"
#include <utilalt.h>
extern twister_base_gen_type mygenerator;

int MCModel::SamplingSub(ImageType3DChar::Pointer grpPtr,
			 VnlVectorImageType::Pointer fmriPtr,
			 ImageType3DChar::Pointer  samplePtr,
			 unsigned subIdx)
{
     using boost::math::cyl_bessel_i;

     double p_acpt = 0;
     double denergy= 0, denergyPrior = 0, denergyLL=0, denergyVMF = 0;
     int cand = 1;
     int currentLabel = 0;
     unsigned clsIdx = 0;
     unsigned grpLabel = 0;

     // images
     VnlVectorImageType::SizeType fmriSize = 
	  fmriPtr->GetLargestPossibleRegion().GetSize();
     VnlVectorImageType::IndexType fmriIdx;    
     fmriIdx.Fill(0);
     IteratorTypeVnlVector fmriIt(fmriPtr, 
				  fmriPtr->GetLargestPossibleRegion().GetSize() );

     // group label
     IteratorType3DChar grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );

     // samples
     ImageType3DChar::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3DChar::IndexType sampleIdx;    
     sampleIdx.Fill(0);
  
     // Define neighborhood iterator
     typedef itk::NeighborhoodIterator< ImageType3DChar> NeighborhoodIteratorType;
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);

     NeighborhoodIteratorType neiSampleIt( radius, samplePtr, 
					   samplePtr->GetLargestPossibleRegion() );

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1}};

     vnl_vector <float> timeSeries(m_tsLength, 0);

     // Uniform integer generator.
     boost::uniform_int<> uni_int(0, m_numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     // debug. Print generator state.
// #pragma omp critical(printstate)
//      {
// 	  std::cout << "mygenerator state: " << mygenerator << "for sub " << subIdx <<  std::endl << std::endl;
//      }


     // compute normalization constant of von Mises
     // Fisher. vmfLogConst[i] is the log (c_d (kappa) ) for the i'th
     // clusters. See page 1350 of "Clustering on the Unit-Sphere
     // using von Mises Fisher..."
     vnl_vector<double> vmfLogConst(m_numClusters);
     float myD = m_tsLength;

     double const Pi = 4 * atan(1);

     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
     	  vmfLogConst[clsIdx] = (myD/2 - 1) * log (m_vmm[subIdx].comp[clsIdx].kappa) - myD/2 * log(2*Pi) -  logBesselI(myD/2 - 1, m_vmm[subIdx].comp[clsIdx].kappa);
     }

     // sampling Markov Random Field. 
     for (neiSampleIt.GoToBegin(), grpIt.GoToBegin(), fmriIt.GoToBegin(); 
	  !neiSampleIt.IsAtEnd(); 
	  ++ neiSampleIt, ++ grpIt, ++fmriIt) {
	  currentLabel = neiSampleIt.GetCenterPixel();
	  if (currentLabel >= 0) {

	       cand = roll_die();

	       // alpha (group) term. Negative log of probabililty.
	       grpLabel = grpIt.Get();

	       denergyPrior = m_alpha * ( int(cand != grpLabel) - int(currentLabel != grpLabel) );
	       denergyLL 
		    = int(cand != neiSampleIt.GetPixel(xminus))
		    - int(currentLabel != neiSampleIt.GetPixel(xminus))

		    + int(cand != neiSampleIt.GetPixel(xplus)) 
		    - int(currentLabel != neiSampleIt.GetPixel(xplus))

		    + int(cand != neiSampleIt.GetPixel(yminus)) 
		    - int(currentLabel != neiSampleIt.GetPixel(yminus))

		    + int(cand != neiSampleIt.GetPixel(yplus)) 
		    - int(currentLabel != neiSampleIt.GetPixel(yplus))

		    + int(cand != neiSampleIt.GetPixel(zminus)) 
		    - int(currentLabel != neiSampleIt.GetPixel(zminus))

		    + int(cand != neiSampleIt.GetPixel(zplus)) 
		    - int(currentLabel != neiSampleIt.GetPixel(zplus));

	       // beta (sub level pairwise interation) term.
	       denergyLL = m_beta_z * denergyLL;


	       // The energy = - Log(likelihood).
	       timeSeries = fmriIt.Get();
	       denergyVMF = (- vmfLogConst[cand] - m_vmm[subIdx].comp[cand].kappa * inner_product(timeSeries, m_vmm[subIdx].comp[cand].mu)) - (- vmfLogConst[currentLabel] - m_vmm[subIdx].comp[currentLabel].kappa * inner_product(fmriIt.Get(), m_vmm[subIdx].comp[currentLabel].mu));

	       // all change of energy. VMF energy weighted by gamma.
	       denergy = denergyPrior + denergyLL + m_gamma * denergyVMF;

	       // system temperature.
	       denergy = denergy / m_temperature;

	       // if energy change less than zero, just accept
	       // candidate. otherwise accept with exp(- energy
	       // change).
			 
	       if (denergy <= 0) {
		    neiSampleIt.SetCenterPixel(cand);
	       }
	       else {
		    p_acpt = exp(-denergy);
		    if (uni() < p_acpt) {
			 neiSampleIt.SetCenterPixel(cand);
		    }
	       }
	  } // in mask
     } // iterators.

     if (m_verbose >= 2) {
	  printf("        mcsampling(): mcSamplingSub() done for subject %i.\n", subIdx);
     }

     return 0;
}


int MCModel::estimateMuSubSample(std::vector<ImageType3DChar::Pointer> & sampleSubVec,
				 VnlVectorImageType::Pointer fmriPtr,
				 ImageType3DChar::Pointer maskPtr,
				 unsigned subIdx)
{
     unsigned int clsIdx = 0;
     
     // images
     VnlVectorImageType::SizeType fmriSize = fmriPtr->GetLargestPossibleRegion().GetSize();
     IteratorTypeVnlVector fmriIt(fmriPtr, 
				   fmriPtr->GetLargestPossibleRegion().GetSize() );

     // mask image.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DChar::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     ImageType3DChar::IndexType maskIdx;

     vnl_vector <float> timeSeries(m_tsLength, 0);

     // reset all mu to zero.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_vmm[subIdx].comp[clsIdx].numPoints = 0;
	  m_vmm[subIdx].comp[clsIdx].mu = 0;
     }

     for (unsigned short sampleIdx = 0; sampleIdx < m_numSamples; sampleIdx ++) {
	  ConstIteratorType3DChar sampleIt(sampleSubVec[sampleIdx], 
					   sampleSubVec[sampleIdx]->GetLargestPossibleRegion() );
	  for (fmriIt.GoToBegin(), maskIt.GoToBegin(), sampleIt.GoToBegin(); 
	       !maskIt.IsAtEnd(); 
	       ++ fmriIt, ++ maskIt, ++ sampleIt) {
	       if (maskIt.Get() > 0) {
		    maskIdx = maskIt.GetIndex();
		    clsIdx = sampleIt.Get();
		    m_vmm[subIdx].comp[clsIdx].mu += fmriIt.Get();
		    m_vmm[subIdx].comp[clsIdx].numPoints ++;
	       }  // maskIt > 0
	  } // for maskIt

     } // sampleIdx

     // compute mean time series and meanNorm.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  // meanNorm see page 1350 of "Clustering on the Unit
	  // Hypersphere using von Mises-Fisher distributions" by
	  // Banerjee.
	  if (m_vmm[subIdx].comp[clsIdx].numPoints > 0) {
	       m_vmm[subIdx].comp[clsIdx].meanNorm = m_vmm[subIdx].comp[clsIdx].mu.two_norm() / m_vmm[subIdx].comp[clsIdx].numPoints;
	       m_vmm[subIdx].comp[clsIdx].mu.normalize();
	  }
	  else {
	       m_vmm[subIdx].comp[clsIdx].meanNorm = 0;
	  }
     }
     return 0;
}

// Compute 1st and 2nd derivatives of Q1 (-log P(G|Z) ) with beta, for each subject.
double MCModel::SubBetaDrvPerSub(std::vector<ImageType3DChar::Pointer> &  grpVec,
				 std::vector<ImageType3DChar::Pointer> & sampleSubVec,
				 ImageType3DChar::Pointer maskPtr,
				 double & drv1, double & drv2)
{
     // mask image.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     MyBoundCondType constCondition;
     constCondition.SetConstant(-1); 
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1}};

     double a = 0, b = 0, M0 = 0, M1 = 0, M2 = 0, acur = 0, bcur = 0;
     unsigned curSubLabel = 0, grpLabel = 0;
     unsigned clsIdx = 0;
     drv1 = 0, drv2 = 0;
     double Q2 = 0;
     for (unsigned short sampleIdx = 0; sampleIdx < m_numSamples; sampleIdx ++) {
	  IteratorType3DChar grpIt(grpVec[sampleIdx], grpVec[sampleIdx]->GetLargestPossibleRegion() );

	  NeighborhoodIteratorType neiSampleIt( radius, sampleSubVec[sampleIdx], 
						sampleSubVec[sampleIdx]->GetLargestPossibleRegion() );
	  neiSampleIt.OverrideBoundaryCondition(&constCondition);

	  for (neiSampleIt.GoToBegin(), grpIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ neiSampleIt, ++ grpIt, ++ maskIt) {
	       if (maskIt.Get() > 0) {
		    curSubLabel = neiSampleIt.GetCenterPixel();

		    M0 = 0, M1 = 0, M2 = 0;
		    for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
			 a = - int (grpIt.Get() != clsIdx);
			 b = - (int(neiSampleIt.GetPixel(xminus) >= 0 && clsIdx != neiSampleIt.GetPixel(xminus) )
				+ int(neiSampleIt.GetPixel(xplus) >= 0 && clsIdx != neiSampleIt.GetPixel(xplus))
				+ int(neiSampleIt.GetPixel(yminus) >= 0 && clsIdx != neiSampleIt.GetPixel(yminus))
				+ int(neiSampleIt.GetPixel(yplus) >= 0 && clsIdx != neiSampleIt.GetPixel(yplus))
				+ int(neiSampleIt.GetPixel(zminus) >= 0 && clsIdx != neiSampleIt.GetPixel(zminus))
				+ int(neiSampleIt.GetPixel(zplus) >= 0 && clsIdx != neiSampleIt.GetPixel(zplus)) );

			 if (clsIdx == curSubLabel) {
			      acur = a;
			      bcur = b;
			 }
			 M0 += exp (a * m_alpha + b * m_beta_z);
			 M1 += b * exp(a *m_alpha + b * m_beta_z);
			 M2 += b * b * exp(a * m_alpha + b * m_beta_z);
		    } // for clsIdx
		    drv1 += (- bcur + M1 / M0);
		    drv2 += (M2 * M0 - M1 * M1)/ (M0 * M0);
		    Q2 = Q2 + (- acur * m_alpha - bcur * m_beta_z + log (M0) )/ m_numSamples;
	       }  // maskIt > 0
	  } // for maskIt
     } // sampleIdx

     return Q2;
}



double MCModel::SubAlphaDrvPerSub(std::vector<ImageType3DChar::Pointer> &  grpVec,
				  std::vector<ImageType3DChar::Pointer>  & sampleSubVec,
				  ImageType3DChar::Pointer maskPtr,
				  double & drv1, double & drv2)
{
     // mask image.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     MyBoundCondType constCondition;
     constCondition.SetConstant(-1);
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1}};

     double a = 0, b = 0, M0 = 0, M1 = 0, M2 = 0, acur = 0, bcur = 0;
     unsigned curSubLabel = 0;
     unsigned clsIdx = 0;
     drv1 = 0, drv2 = 0;
     double Q2 = 0;
     for (unsigned short mcIdx = 0; mcIdx < m_numSamples; mcIdx ++) {
	  IteratorType3DChar grpIt(grpVec[mcIdx], grpVec[mcIdx]->GetLargestPossibleRegion() );
	  NeighborhoodIteratorType neiSampleIt( radius, sampleSubVec[mcIdx], 
						sampleSubVec[mcIdx]->GetLargestPossibleRegion() );
	  neiSampleIt.OverrideBoundaryCondition(&constCondition);

	  for (neiSampleIt.GoToBegin(), grpIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ neiSampleIt, ++ grpIt, ++ maskIt) {
	       if (maskIt.Get() > 0) {
		    curSubLabel = neiSampleIt.GetCenterPixel();

		    M0 = 0, M1 = 0, M2 = 0;
		    for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
			 a = - int (grpIt.Get() != clsIdx);
			 b = - (int(clsIdx != neiSampleIt.GetPixel(xminus))
				+ int(clsIdx != neiSampleIt.GetPixel(xplus))
				+ int(clsIdx != neiSampleIt.GetPixel(yminus))
				+ int(clsIdx != neiSampleIt.GetPixel(yplus))
				+ int(clsIdx != neiSampleIt.GetPixel(zminus))
				+ int(clsIdx != neiSampleIt.GetPixel(zplus)) );
			 if (clsIdx == curSubLabel) {
			      acur = a;
			      bcur = b;
			 }
			 M0 += exp (a * m_alpha + b * m_beta_z);
			 M1 += a * exp(a * m_alpha + b * m_beta_z);
			 M2 += a * a * exp(a * m_alpha + b * m_beta_z);
		    } // for clsIdx
		    drv1 += (- acur + M1 / M0);
		    drv2 += (M2 * M0 - M1 * M1)/ (M0 * M0);
		    Q2 = Q2 + (-acur * m_alpha - bcur * m_beta_z + log (M0) ) / m_numSamples;
	       }  // maskIt > 0
	  } // for maskIt
     } // mcIdx
     return Q2;
}
