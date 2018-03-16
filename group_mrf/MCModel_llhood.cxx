#include "commonalt.h"
#include "MCModel.h"
#include <boost/math/special_functions/gamma.hpp>
#include <utilalt.h>
extern twister_base_gen_type mygenerator;


double MCModel::ComputeQ1(ImageType3DChar::Pointer &  grpPtr,
		   ImageType3DVecUS::Pointer & histPtr,
		   ImageType3DChar::Pointer maskPtr)

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

     double a = 0, b = 0, M0 = 0, acur = 0, bcur = 0;
     unsigned curLabel = 0;
     unsigned clsIdx = 0;
     double Q1 = 0;

     NeighborhoodIteratorType neiGrpIt(radius, grpPtr, grpPtr->GetLargestPossibleRegion() );
     neiGrpIt.OverrideBoundaryCondition(&constCondition);
     IteratorType3DVecUS histIt(histPtr, histPtr->GetLargestPossibleRegion() );

     for (neiGrpIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ neiGrpIt, ++maskIt) {
	  if (maskIt.Get() > 0) {
	       curLabel = neiGrpIt.GetCenterPixel();
	       M0 = 0;
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
		    a = (double)(histIt.Get()[clsIdx]) - m_numSubs;
		    b = - (int(neiGrpIt.GetPixel(xminus) >= 0 && clsIdx != neiGrpIt.GetPixel(xminus) )
			   + int(neiGrpIt.GetPixel(xplus) >= 0 && clsIdx != neiGrpIt.GetPixel(xplus))
			   + int(neiGrpIt.GetPixel(yminus) >= 0 && clsIdx != neiGrpIt.GetPixel(yminus))
			   + int(neiGrpIt.GetPixel(yplus) >= 0 && clsIdx != neiGrpIt.GetPixel(yplus))
			   + int(neiGrpIt.GetPixel(zminus) >= 0 && clsIdx != neiGrpIt.GetPixel(zminus))
			   + int(neiGrpIt.GetPixel(zplus) >= 0 && clsIdx != neiGrpIt.GetPixel(zplus)) );

		    if (clsIdx == curLabel) {
			 acur = a;
			 bcur = b;
		    }
		    M0 += exp (a * m_alpha + b * m_beta_g);
	       } // for clsIdx
	       Q1 = Q1 + (- acur * m_alpha - bcur * m_beta_g + log(M0));
	  }  // maskIt > 0
     } // for maskIt
     
     return Q1;
}

double MCModel::ComputeQ2(ImageType3DChar::Pointer &  grpPtr,
		   ImageType3DChar::Pointer & samplePtr,
		   ImageType3DChar::Pointer maskPtr,
		   unsigned subIdx)

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

     double a = 0, b = 0, M0 = 0, acur = 0, bcur = 0;
     unsigned curSubLabel = 0;
     unsigned clsIdx = 0;
     double Q2 = 0;
     IteratorType3DChar grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );
     NeighborhoodIteratorType neiSampleIt( radius, samplePtr, samplePtr->GetLargestPossibleRegion() );
     neiSampleIt.OverrideBoundaryCondition(&constCondition);

     for (neiSampleIt.GoToBegin(), grpIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ neiSampleIt, ++ grpIt, ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       curSubLabel = neiSampleIt.GetCenterPixel();

	       M0 = 0;
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
	       } // for clsIdx
	       Q2 = Q2 + (- acur * m_alpha - bcur * m_beta_z + log (M0) );
	  }  // maskIt > 0
     } // for maskIt
     return Q2;
}

// - log P(X|Z) for single sample, single subject.
double MCModel::ComputeQ3(VnlVectorImageType::Pointer fmriPtr,
		   ImageType3DChar::Pointer samplePtr,
		   ImageType3DChar::Pointer maskPtr,
		   unsigned subIdx)
{
     using boost::math::cyl_bessel_i;
     unsigned clsIdx = 0;
     double logPXZ = 0;

     // images
     IteratorTypeVnlVector fmriIt(fmriPtr, 
				  fmriPtr->GetLargestPossibleRegion().GetSize() );
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     // samples
     ImageType3DChar::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3DChar::IndexType sampleIdx;    
     sampleIdx.Fill(0);
     ConstIteratorType3DChar sampleIt(samplePtr, 
				      samplePtr->GetLargestPossibleRegion() );

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

     vnl_vector <float> timeSeries(m_tsLength, 0);
     for (maskIt.GoToBegin(), sampleIt.GoToBegin(), fmriIt.GoToBegin(); !sampleIt.IsAtEnd(); ++ maskIt, ++ sampleIt, ++ fmriIt) {
	  if (maskIt.Get() > 0 ) {
	       timeSeries = fmriIt.Get();
	       clsIdx = sampleIt.Get();
	       logPXZ += ( vmfLogConst[clsIdx] + m_vmm[subIdx].comp[clsIdx].kappa * inner_product(timeSeries, m_vmm[subIdx].comp[clsIdx].mu) );
	  }
     }

     return (- m_gamma * logPXZ);
}

double MCModel::ComputeEnergy(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
			      std::vector<ImageType3DChar::Pointer> &  grpVec,
			      ImageType3DChar::Pointer maskPtr,
			      std::vector<VnlVectorImageType::Pointer> & fmriVec)
{
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());

     // create histogram data block.
     itk::VariableLengthVector<unsigned short> histVector(m_numClusters);
     histVector.Fill(0);
     std::vector <ImageType3DVecUS::Pointer> histVec(m_numSamples);
     ImageType3DVecUS::RegionType histRegion = maskPtr->GetLargestPossibleRegion();
     for (unsigned mcIdx = 0; mcIdx < m_numSamples; mcIdx ++) {
	  histVec[mcIdx] = ImageType3DVecUS::New();
	  histVec[mcIdx]->SetRegions(histRegion);
	  histVec[mcIdx]->SetNumberOfComponentsPerPixel(m_numClusters);
	  histVec[mcIdx]->SetVectorLength(m_numClusters);
	  histVec[mcIdx]->Allocate();
	  histVec[mcIdx]->FillBuffer(histVector);
     }

     // compute histogram
     for (unsigned short mcIdx = 0; mcIdx < m_numSamples; mcIdx ++) {
	  IteratorType3DVecUS histIt(histVec[mcIdx], histRegion);
	  for (unsigned subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       ConstIteratorType3DChar sampleIt(sampleVec[subIdx][mcIdx], sampleVec[subIdx][mcIdx]->GetLargestPossibleRegion() );

	       for (sampleIt.GoToBegin(), histIt.GoToBegin(), maskIt.GoToBegin(); 
		    !maskIt.IsAtEnd(); ++ histIt, ++sampleIt, ++ maskIt) {
		    if (maskIt.Get() > 0) {
			 histVector = histIt.Get();
			 histVector[sampleIt.Get()] ++;
			 histIt.Set(histVector);
		    }
	       } // for sampelIt.
	  } // for subIdx
     } // mcIdx

     vnl_vector<float> Q1(m_numSamples, 0);
     vnl_vector<float> Q2(m_numSamples, 0);
     vnl_vector<float> Q3(m_numSamples, 0);

#pragma omp parallel for
     for (unsigned mcIdx = 0; mcIdx < m_numSamples; mcIdx ++) {
	  Q1[mcIdx] = ComputeQ1(grpVec[mcIdx], histVec[mcIdx], maskPtr);
     }

     double subQ2 = 0;
#pragma omp parallel for
     for (unsigned mcIdx = 0; mcIdx < m_numSamples; mcIdx ++) {
	  for (unsigned subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       subQ2 = ComputeQ2(grpVec[mcIdx], sampleVec[subIdx][mcIdx], maskPtr, subIdx); 
	       Q2[mcIdx] += subQ2;
	  }
     }

     // Q3. -log P(X|Z)
     double subQ3 = 0;
#pragma omp parallel for
     for (unsigned mcIdx = 0; mcIdx < m_numSamples; mcIdx ++) {
	  for (unsigned subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       subQ3 = ComputeQ3(fmriVec[subIdx], sampleVec[subIdx][mcIdx], maskPtr, subIdx);
	       Q3[mcIdx] += subQ3;
	  }
     }

     printf("Q1[mcIdx]: ");
     printVnlVector(Q1, 100);
     printf("Q2[mcIdx]: ");
     printVnlVector(Q2, 100);
     printf("Q3[mcIdx]: ");
     printVnlVector(Q3, 100);

     printf("Q1/E(G)  = %E, Q2/E(Z|G) = %E, Q3/E(X|Z) = %E, Q2+Q3 = %E, Q/E(X,G,Z) = %E\n", Q1.mean(), Q2.mean(), Q3.mean(), Q2.mean() + Q3.mean(), Q1.mean() + Q2.mean() + Q3.mean() );
     return (Q1.mean() + Q2.mean() + Q3.mean());
}
