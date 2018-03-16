#include "common_mcem.h"
#include "gmmodel_mcem.h"
#include <boost/math/distributions/students_t.hpp>
extern twister_base_gen_type mygenerator; 
typedef itk::ConstantBoundaryCondition< ImageType3DChar >  MyBoundCondType;

GMModel::GMModel(std::vector<ImageType3DFloat::Pointer> & obsVec,
		 std::vector<ImageType3DChar::Pointer> &  grpVec,
		 std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		 ImageType3DChar::Pointer maskPtr,
		 ParStruct par)
{
     unsigned clsIdx = 0;

     // init members.
     m_verbose = par.verbose;
     m_beta_g = par.betag;
     m_beta_z = par.betaz;
     m_alpha = par.alpha;
     m_gamma = par.gamma;
     m_temperature = 1; // init to 1.
     m_numSamples = par.numSamples;
     m_numSubs = sampleVec.size();
     m_numClusters = par.numClusters;
     m_verbose = par.verbose;

     // initial total number of data points (voxels) in the mask.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());
     ImageType3DChar::IndexType maskIdx;     
     m_gmm.totalPts = 0;
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       m_gmm.totalPts ++;
	  }  // maskIt > 0
     } // maskIt

     if(m_verbose >= 1) {
	  printf("total number of voxels within mask: %i.\n", m_gmm.totalPts);
     }

     // init m_gmm
     m_gmm.comp.resize(m_numClusters);
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_gmm.comp[clsIdx].prop.set_size(m_numSubs);
	  m_gmm.comp[clsIdx].prop = 1 / float (m_numClusters);

	  m_gmm.comp[clsIdx].numPoints.set_size(m_numSubs);
	  m_gmm.comp[clsIdx].numPoints = m_gmm.totalPts / float (m_numClusters);

	  m_gmm.comp[clsIdx].mu.set_size(m_numSubs);
	  m_gmm.comp[clsIdx].sigma.set_size(m_numSubs);
     }
     
}

GMModel::~GMModel()
{
     // delete m_gc;
}

int GMModel::printSelf(std::string format)
{     unsigned int clsIdx = 0, subIdx = 0;

     if (format.compare("normal") == 0) {

	  printf("beta_g = %f, beta_z = %f, alpha = %f, temperature = %f\n", m_beta_g, m_beta_z, m_alpha, m_temperature);

	  for (clsIdx = 0; clsIdx < m_numClusters; clsIdx++) {
	       printf("comp[%i]: \n", clsIdx+1);
	       printf("  numPoints = ");
	       for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
		    printf("%.1f ", m_gmm.comp[clsIdx].numPoints[subIdx]);
	       }
	       printf("\n");
	       
	       printf("  prop = ");	       
	       for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
		    printf("%2.2f ", m_gmm.comp[clsIdx].prop[subIdx]);
	       }
	       printf("\n");


	       printf("  mu = ");
	       for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
		    printf("%2.2f ", m_gmm.comp[clsIdx].mu[subIdx]);
	       }
	       printf("\n");
	       printf("  sigma = ");
	       for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
		    printf("%2.2f ", m_gmm.comp[clsIdx].sigma[subIdx]);
	       }
	       printf("\n");
	  } // clsIdx.
     }
     else if (format.compare("colname") == 0) {
	  printf("%10s %10s %10s %10s %10s", "TB:",   "betag",   "betaz", "alpha", "temp");
	  for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
		    printf("  mu_c%i_s%i", clsIdx+1, subIdx+1);
		    printf("  si_c%i_s%i", clsIdx+1, subIdx+1);
	       }
	  }
	  printf("\n");
     }
	  
     else if (format.compare("table") == 0) {
	  printf("%10s %10.3f %10.3f %10.3f %10.3f", "TB:", m_beta_g, m_beta_z, m_alpha, m_temperature);
	  for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
		    printf("%10.3f", m_gmm.comp[clsIdx].mu[subIdx]);
		    printf("%10.3f", m_gmm.comp[clsIdx].sigma[subIdx]);
	       }
	  }
	  printf("\n");
     }
     else {
	  printf("GMModel::print(): print format not recognized.\n");
     }
     fflush(stdout);
     
     return 0;
}


int GMModel::Sampling(std::vector<ImageType3DChar::Pointer> & grpVec,
		      std::vector<ImageType3DFloat::Pointer> & obsVec,
		      ImageType3DChar::Pointer maskPtr,
		      std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		      unsigned burnin)
{
     unsigned subIdx = 0, scanIdx = 0;;
     for (scanIdx = 0; scanIdx < burnin + m_numSamples; scanIdx ++) {
	  // flip between sampling of group and individual subs.
	  // sampling Z (subject's label)
#pragma omp parallel for
	  for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       SamplingSub(grpVec[m_numSamples - 1], obsVec[subIdx], maskPtr, sampleVec[subIdx][m_numSamples - 1], subIdx);
	  }
	  if (m_verbose >= 3) {
	       printf("      Sampling(): scan %i subject level done.\n", scanIdx + 1);
	  }

	  // sampling G (group label)
	  SamplingGrp(grpVec[m_numSamples - 1], maskPtr, sampleVec);

	  if (m_verbose >= 3) {
	       printf("      Sampling(): scan %i group level done.\n", scanIdx + 1);
	  }

	  // After burnin period, save it to correct place.
	  if (scanIdx >= burnin) {
	       for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
		    // deinfe src iterator for copy mc samples.
		    IteratorType3DChar sampleSrcIt(sampleVec[subIdx][m_numSamples - 1], 
						   sampleVec[subIdx][m_numSamples - 1]->GetLargestPossibleRegion() );

		    IteratorType3DChar sampleDestIt(sampleVec[subIdx][scanIdx - burnin], 
						    sampleVec[subIdx][scanIdx - burnin]->GetLargestPossibleRegion() );

		    for (sampleSrcIt.GoToBegin(), sampleDestIt.GoToBegin(); !sampleSrcIt.IsAtEnd(); ++ sampleSrcIt, ++ sampleDestIt) {
			 sampleDestIt.Set(sampleSrcIt.Get());
		    }
	       } // end for.

	       // also save grp label map.
	       IteratorType3DChar srcGrpIt(grpVec[m_numSamples-1], grpVec[m_numSamples-1]->GetLargestPossibleRegion() );
	       IteratorType3DChar destGrpIt(grpVec[scanIdx-burnin], grpVec[scanIdx-burnin]->GetLargestPossibleRegion() );	       
	       for (srcGrpIt.GoToBegin(), destGrpIt.GoToBegin(); !srcGrpIt.IsAtEnd(); ++ srcGrpIt, ++ destGrpIt) {
			 destGrpIt.Set(srcGrpIt.Get());
	       }
	       
	  } // scanIdx >= burnin
	  if (m_verbose >= 1 && (scanIdx%10 == 0) ) {
	       printf("mcsampling(): scan %i done.\n", scanIdx);
	  }
     } // for scanIdx.

}

int GMModel::SamplingSub(ImageType3DChar::Pointer grpPtr,
			 ImageType3DFloat::Pointer obsPtr,
			 ImageType3DChar::Pointer maskPtr,
			 ImageType3DChar::Pointer samplePtr,
			 unsigned subIdx)
{

     double p_acpt = 0;
     double denergy= 0, denergyPrior = 0, denergyLL=0, denergyVMF = 0;
     int cand = 1;
     int currentLabel = 0;
     unsigned clsIdx = 0;
     double dlogll = 0;
     unsigned grpLabel = 0;



     ImageType3DFloat::IndexType obsIdx;    
     ImageType3DFloat::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DFloat obsIt(obsPtr, obsPtr->GetLargestPossibleRegion() );

     //mask
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType maskIdx;     

     // group map.
     IteratorType3DChar grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );
     itk::VariableLengthVector<float> grpVector(m_numClusters);

     // samples
     ImageType3DChar::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3DChar::IndexType sampleIdx;    
     sampleIdx.Fill(0);

     // Define neighborhood iterator
     MyBoundCondType constCondition;
     constCondition.SetConstant(-1);
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType neiSampleIt( radius, samplePtr, 
					   samplePtr->GetLargestPossibleRegion() );
     neiSampleIt.OverrideBoundaryCondition(&constCondition);

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1}};

     // Uniform integer generator.
     boost::uniform_int<> uni_int(0, m_numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);

     // sampling Markov Random Field. 
     for (neiSampleIt.GoToBegin(), grpIt.GoToBegin(), obsIt.GoToBegin(); 
	  !neiSampleIt.IsAtEnd(); 
	  ++ neiSampleIt, ++ grpIt, ++obsIt) {
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
	       denergyVMF = log (m_gmm.comp[cand].sigma[subIdx]) + pow( (obsIt.Get() - m_gmm.comp[cand].mu[subIdx]), 2) / ( 2 * pow(m_gmm.comp[cand].sigma[subIdx], 2))
		    - log (m_gmm.comp[currentLabel].sigma[subIdx]) - pow( (obsIt.Get() - m_gmm.comp[currentLabel].mu[subIdx]), 2) / ( 2 * pow(m_gmm.comp[currentLabel].sigma[subIdx], 2));

	       denergy = denergyPrior + denergyLL + m_gamma * denergyVMF;


	       // temperature take effect.
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
}


int GMModel::SamplingGrp(ImageType3DChar::Pointer grpPtr,
			 ImageType3DChar::Pointer maskPtr,
			 std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec)
{
     double p_acpt = 0;
     double denergy= 0;
     int cand = 1;
     int currentLabel = 0;

     //mask
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     // Define neighborhood iterator
     MyBoundCondType constCondition;
     constCondition.SetConstant(-1);     
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType neiGrpIt(radius, grpPtr, grpPtr->GetLargestPossibleRegion() );
     neiGrpIt.OverrideBoundaryCondition(&constCondition);

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1}};

     // create histogram data block.
     itk::VariableLengthVector<unsigned short> histVector(m_numClusters);

     ImageType3DVecUS::Pointer histPtr = ImageType3DVecUS::New();
     ImageType3DVecUS::RegionType histRegion = grpPtr->GetLargestPossibleRegion();
     histPtr->SetRegions(histRegion);
     histPtr->SetNumberOfComponentsPerPixel(m_numClusters);
     histPtr->SetVectorLength(m_numClusters);
     histPtr->Allocate();

     histVector.Fill(0);
     histPtr->FillBuffer(histVector);
     IteratorType3DVecUS histIt(histPtr, histRegion);

     // compute histogram
     for (unsigned short subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	  ConstIteratorType3DChar sampleIt(sampleVec[subIdx][m_numSamples - 1], 
					   sampleVec[subIdx][m_numSamples - 1]->GetLargestPossibleRegion() );

	  for (sampleIt.GoToBegin(), histIt.GoToBegin(); 
	       !histIt.IsAtEnd(); ++ histIt, ++sampleIt) {
	       if (sampleIt.Get() >= 0) {
		    histVector = histIt.Get();
		    histVector[sampleIt.Get()] ++;
		    histIt.Set(histVector);
	       }
	  } // for sampelIt.
     } // for subIdx

     // Uniform integer generator.
     boost::uniform_int<> uni_int(0, m_numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);
     
     for (neiGrpIt.GoToBegin(), histIt.GoToBegin(), maskIt.GoToBegin(); 
	  !histIt.IsAtEnd(); ++ histIt, ++ neiGrpIt, ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       currentLabel = neiGrpIt.GetCenterPixel();
	       cand = roll_die();
	       denergy = m_alpha * (histIt.Get()[currentLabel] - histIt.Get()[cand])

		    + (int(cand != neiGrpIt.GetPixel(xminus))
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
		       - int(currentLabel != neiGrpIt.GetPixel(zplus)) ) * m_beta_g;	       


	       // temperature take effect.
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



int GMModel::estimatePar(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
			 std::vector<ImageType3DFloat::Pointer> & obsVec,
			 ImageType3DChar::Pointer maskPtr)
{
     unsigned int clsIdx = 0, subIdx = 0;
     // compute mu and numPoints for all sub, all clusters.
     for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	  estimateMuSub(sampleVec[subIdx], obsVec[subIdx], maskPtr, subIdx);	  
     } // for subIdx


     // compute sigma for all sub, all clusters.
     for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	  estimateSigmaSub(sampleVec[subIdx], obsVec[subIdx], maskPtr, subIdx);	  
     } // for subIdx
     return 0;
}


int GMModel::estimateMuSub(std::vector<ImageType3DChar::Pointer> & sampleSubVec,
			   ImageType3DFloat::Pointer  obsPtr,
			   ImageType3DChar::Pointer maskPtr,
			   unsigned subIdx)
{
     unsigned int clsIdx = 0;
     unsigned grpLabel = 0;

     // images
     VnlVectorImageType::IndexType obsIdx;    
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DFloat obsIt(obsPtr, 
				   obsPtr->GetLargestPossibleRegion().GetSize() );

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );


     // reset all mu to zero.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_gmm.comp[clsIdx].numPoints[subIdx] = 0;
	  m_gmm.comp[clsIdx].mu[subIdx] = 0;
     }

     for (unsigned short sampleIdx = 0; sampleIdx < m_numSamples; sampleIdx ++) {
	  ConstIteratorType3DChar sampleIt(sampleSubVec[sampleIdx], 
					   sampleSubVec[sampleIdx]->GetLargestPossibleRegion() );
	  for (obsIt.GoToBegin(), maskIt.GoToBegin(), sampleIt.GoToBegin(); 
	       !maskIt.IsAtEnd(); 
	       ++ obsIt, ++ maskIt, ++ sampleIt) {
	       if (maskIt.Get() > 0) {
		    clsIdx = sampleIt.Get();
		    m_gmm.comp[clsIdx].mu[subIdx] += obsIt.Get();
		    m_gmm.comp[clsIdx].numPoints[subIdx] ++;
	       }  // maskIt > 0
	  } // for maskIt

     } // sampleIdx


     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_gmm.comp[clsIdx].mu[subIdx] = m_gmm.comp[clsIdx].mu[subIdx] / m_gmm.comp[clsIdx].numPoints[subIdx];
	  m_gmm.comp[clsIdx].prop[subIdx] = m_gmm.comp[clsIdx].numPoints[subIdx] / (m_gmm.totalPts * m_numSamples);
     }
     return 0;
}


int GMModel::estimateSigmaSub(std::vector<ImageType3DChar::Pointer> & sampleSubVec,
			      ImageType3DFloat::Pointer  obsPtr,
			      ImageType3DChar::Pointer maskPtr,
			      unsigned subIdx)
{
     unsigned int clsIdx = 0;
     unsigned grpLabel = 0;

     // images
     VnlVectorImageType::IndexType obsIdx;    
     VnlVectorImageType::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DFloat obsIt(obsPtr, 
				   obsPtr->GetLargestPossibleRegion().GetSize() );

     // mask image.
     ImageType3DChar::IndexType maskIdx;     
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );


     // reset all sigma to zero.
     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_gmm.comp[clsIdx].sigma[subIdx] = 0;
     }

     for (unsigned short sampleIdx = 0; sampleIdx < m_numSamples; sampleIdx ++) {
	  ConstIteratorType3DChar sampleIt(sampleSubVec[sampleIdx], 
					   sampleSubVec[sampleIdx]->GetLargestPossibleRegion() );
	  for (obsIt.GoToBegin(), maskIt.GoToBegin(), sampleIt.GoToBegin(); 
	       !maskIt.IsAtEnd(); 
	       ++ obsIt, ++ maskIt, ++ sampleIt) {
	       if (maskIt.Get() > 0) {
		    clsIdx = sampleIt.Get();
		    m_gmm.comp[clsIdx].sigma[subIdx] += pow( (obsIt.Get() - m_gmm.comp[clsIdx].mu[subIdx]), 2);
	       }  // maskIt > 0
	  } // for maskIt

     } // sampleIdx


     for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_gmm.comp[clsIdx].sigma[subIdx] = sqrt ( m_gmm.comp[clsIdx].sigma[subIdx] / m_gmm.comp[clsIdx].numPoints[subIdx]);
     }
     return 0;
}


double GMModel::EstimatePriorPar(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
				 std::vector<ImageType3DChar::Pointer> &  grpVec,
				 ImageType3DChar::Pointer maskPtr,
				 std::string position) // which location to run this fun.
{

     double beta_g_old = m_beta_g, beta_z_old = m_beta_z, alpha_old = m_alpha;
     double drv1 = 0, drv2 = 0;
     double drv1g = 0, drv1z = 0, drv2g = 0, drv2z = 0; // for alpha.
     
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

     // Estimate alpha, beta.
     double Q1 = 0, Q2 = 0, Q1_old = 0, Q2_old = 0;;
     unsigned numIter = 0;
     do {
	  Q1_old = GrpBetaDrv(grpVec, sampleVec, histVec, maskPtr, drv1, drv2);
	  beta_g_old = m_beta_g;
	  m_beta_g = beta_g_old - drv1/drv2;
	  Q1 = GrpBetaDrv(grpVec, sampleVec, histVec, maskPtr, drv1, drv2);

	  if (m_verbose >= 2) {
	       printf("    mbeta_g: %5.4f -> %5.4f. Q1: %.3f -> %.3f. drv1 = %.4f, drv2 = %.4f\n", beta_g_old, m_beta_g, Q1_old, Q1, drv1, drv2);
	  }

     	  Q2_old = SubBetaDrv(grpVec, sampleVec, maskPtr, drv1, drv2);
     	  beta_z_old = m_beta_z;
     	  m_beta_z = beta_z_old - drv1/drv2;
     	  Q2 = SubBetaDrv(grpVec, sampleVec, maskPtr, drv1, drv2);

	  if (m_verbose >= 2) {
	       printf("    beta_z: %5.4f -> %5.4f. Q2: %.3f -> %.3f. drv1 = %.4f, drv2 = %.4f\n", beta_z_old, m_beta_z, Q2_old, Q2, drv1, drv2);
	  }

	  if (position.compare("init") == 0) {
	       // alpha unchanged. Same with initial value.
	  }
	  else if (position.compare("normal") == 0) {
	       // GrpAlphaDrv(grpVec, sampleVec, histVec, maskPtr, drv1g, drv2g);
	       // SubAlphaDrv(grpVec, sampleVec, maskPtr, drv1z, drv2z);
	       // alpha_old = m_alpha;

	       // // the derivative of (q1+q2) is equal to the sum of their
	       // // derivatives.
	       // m_alpha = alpha_old - (drv1g + drv1z)/(drv2g+drv2z);
	       // Q1 = GrpAlphaDrv(grpVec, sampleVec, histVec, maskPtr, drv1g, drv2g);
	       // Q2 = SubAlphaDrv(grpVec, sampleVec, maskPtr, drv1z, drv2z);
	  }
	  else {
	       printf("EstimatePriorPar(): position variable have wrong value.\n");
	       exit(1);
	  }

	  printf("beta_g = %5.4f, beta_z = %5.4f, alpha = %5.4f, Q1 (E(G)) = %.3f, Q2 (E(Z|G) = %.3f, (Q1+Q2) = %.3f\n", m_beta_g, m_beta_z, m_alpha, Q1, Q2,(Q1+Q2));
	  numIter ++;
     } 
     while (fabs(m_beta_z - beta_z_old) + fabs(m_beta_g - beta_g_old) + fabs(m_alpha - alpha_old) > 1e-5 && numIter <= 2); 

     return (Q1+Q2);
}




double GMModel::GrpBetaDrv(std::vector<ImageType3DChar::Pointer> &  grpVec,
			std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
			std::vector <ImageType3DVecUS::Pointer> & histVec,
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
     unsigned curLabel = 0;
     unsigned clsIdx = 0;
     drv1 = 0, drv2 = 0;
     double Q1 = 0;
     for (unsigned short sampleIdx = 0; sampleIdx < m_numSamples; sampleIdx ++) {
	  NeighborhoodIteratorType neiGrpIt(radius, grpVec[sampleIdx], grpVec[sampleIdx]->GetLargestPossibleRegion() );
	  neiGrpIt.OverrideBoundaryCondition(&constCondition);
	  IteratorType3DVecUS histIt(histVec[sampleIdx], histVec[sampleIdx]->GetLargestPossibleRegion() );

	  for (neiGrpIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ neiGrpIt, ++maskIt) {
	       if (maskIt.Get() > 0) {
		    curLabel = neiGrpIt.GetCenterPixel();
		    M0 = 0, M1 = 0, M2 = 0;
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
			 M1 += b * exp(a * m_alpha + b * m_beta_g);
			 M2 += b * b * exp(a * m_alpha + b * m_beta_g);

		    } // for clsIdx
		    drv1 = drv1 - bcur + M1 / M0;
		    drv2 += (M2 * M0 - M1 * M1)/ (M0 * M0);
		    Q1 = Q1 + (- acur * m_alpha - bcur * m_beta_g + log(M0)) / m_numSamples;
	       }  // maskIt > 0
	  } // for maskIt

     } // sampleIdx
     
     return Q1;
}

double GMModel::SubBetaDrv(std::vector<ImageType3DChar::Pointer> &  grpVec,
			std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
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
	  for (unsigned short subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       NeighborhoodIteratorType neiSampleIt( radius, sampleVec[subIdx][sampleIdx], 
						     sampleVec[subIdx][sampleIdx]->GetLargestPossibleRegion() );
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
	  } // subIdx	       
     } // sampleIdx
     
     return Q2;
}





double GMModel::GrpAlphaDrv(std::vector<ImageType3DChar::Pointer> &  grpVec,
			 std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
			 std::vector <ImageType3DVecUS::Pointer> & histVec,
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
     double Q1 = 0;
     for (unsigned short sampleIdx = 0; sampleIdx < m_numSamples; sampleIdx ++) {
	  NeighborhoodIteratorType neiGrpIt(radius, grpVec[sampleIdx], grpVec[sampleIdx]->GetLargestPossibleRegion() );
	  neiGrpIt.OverrideBoundaryCondition(&constCondition);
	  IteratorType3DVecUS histIt(histVec[sampleIdx], histVec[sampleIdx]->GetLargestPossibleRegion() );

	  for (neiGrpIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ neiGrpIt, ++ maskIt) {
	       if (maskIt.Get() > 0) {
		    curSubLabel = neiGrpIt.GetCenterPixel();

		    M0 = 0, M1 = 0, M2 = 0;
		    for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
			 a = (double)(histIt.Get()[clsIdx]) - m_numSubs;
			 b = - (int(neiGrpIt.GetPixel(xminus) >= 0 && clsIdx != neiGrpIt.GetPixel(xminus) )
				+ int(neiGrpIt.GetPixel(xplus) >= 0 && clsIdx != neiGrpIt.GetPixel(xplus))
				+ int(neiGrpIt.GetPixel(yminus) >= 0 && clsIdx != neiGrpIt.GetPixel(yminus))
				+ int(neiGrpIt.GetPixel(yplus) >= 0 && clsIdx != neiGrpIt.GetPixel(yplus))
				+ int(neiGrpIt.GetPixel(zminus) >= 0 && clsIdx != neiGrpIt.GetPixel(zminus))
				+ int(neiGrpIt.GetPixel(zplus) >= 0 && clsIdx != neiGrpIt.GetPixel(zplus)) );

			 if (clsIdx == curSubLabel) {
			      acur = a;
			      bcur = b;
			 }
			 M0 += exp (a * m_alpha + b * m_beta_g);
			 M1 += a * exp(a * m_alpha + b * m_beta_g);
			 M2 += a * a * exp(a * m_alpha + b * m_beta_g);
		    } // for clsIdx
		    drv1 += (- acur + M1 / M0);
		    drv2 += (M2 * M0 - M1 * M1)/ (M0 * M0);
		    Q1 = Q1 + (-acur * m_alpha - bcur * m_beta_g + log(M0)) / m_numSamples;
	       }  // maskIt > 0
	  } // for maskIt
     } // sampleIdx

     return Q1;
}


double GMModel::SubAlphaDrv(std::vector<ImageType3DChar::Pointer> &  grpVec,
			 std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
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
     for (unsigned short mcIdx = 0; mcIdx < m_numSamples; mcIdx ++) {
	  IteratorType3DChar grpIt(grpVec[mcIdx], grpVec[mcIdx]->GetLargestPossibleRegion() );
	  for (unsigned short subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       NeighborhoodIteratorType neiSampleIt( radius, sampleVec[subIdx][mcIdx], 
						     sampleVec[subIdx][mcIdx]->GetLargestPossibleRegion() );
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
	  } // subIdx	       
     } // mcIdx

     return Q2;
}



// ICM on subject j.
float GMModel::ICMSub(ImageType3DChar::Pointer grpPtr,
		      ImageType3DFloat::Pointer obsPtr,
		      ImageType3DChar::Pointer maskPtr,
		      ImageType3DChar::Pointer samplePtr,
		      unsigned subIdx)
{
     vnl_vector<double> denergy(m_numClusters, 0);
     double denergyPrior = 0, denergyLL=0, denergyVMF = 0;

     int currentLabel = 0;
     unsigned clsIdx = 0;
     double dlogll = 0;
     unsigned grpLabel = 0;

     ImageType3DFloat::IndexType obsIdx;    
     ImageType3DFloat::SizeType obsSize = obsPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DFloat obsIt(obsPtr, obsPtr->GetLargestPossibleRegion() );

     //mask
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );
     ImageType3DChar::IndexType maskIdx;     

     // group map.
     IteratorType3DChar grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );
     itk::VariableLengthVector<float> grpVector(m_numClusters);

     // samples
     ImageType3DChar::SizeType sampleSize = 
	  samplePtr->GetLargestPossibleRegion().GetSize();
     ImageType3DChar::IndexType sampleIdx;    
     sampleIdx.Fill(0);

     // Define neighborhood iterator
     MyBoundCondType constCondition;
     constCondition.SetConstant(-1);
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType neiSampleIt( radius, samplePtr, 
					   samplePtr->GetLargestPossibleRegion() );
     neiSampleIt.OverrideBoundaryCondition(&constCondition);

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1}};


     unsigned numFlips = 0; // number of voxed changed in this scan.
     unsigned newLabel = 0;
     for (neiSampleIt.GoToBegin(), grpIt.GoToBegin(), obsIt.GoToBegin(); 
	  !neiSampleIt.IsAtEnd(); 
	  ++ neiSampleIt, ++ grpIt, ++obsIt) {
	  currentLabel = neiSampleIt.GetCenterPixel();
	  if (currentLabel >= 0) {
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {

		    // alpha (group) term. Negative log of probabililty.
		    grpLabel = grpIt.Get();
		    denergyPrior = m_alpha * ( int(clsIdx != grpLabel) - int(currentLabel != grpLabel) );
		    denergyLL 
			 = int(clsIdx != neiSampleIt.GetPixel(xminus))
			 - int(currentLabel != neiSampleIt.GetPixel(xminus))

			 + int(clsIdx != neiSampleIt.GetPixel(xplus)) 
			 - int(currentLabel != neiSampleIt.GetPixel(xplus))

			 + int(clsIdx != neiSampleIt.GetPixel(yminus)) 
			 - int(currentLabel != neiSampleIt.GetPixel(yminus))

			 + int(clsIdx != neiSampleIt.GetPixel(yplus)) 
			 - int(currentLabel != neiSampleIt.GetPixel(yplus))

			 + int(clsIdx != neiSampleIt.GetPixel(zminus)) 
			 - int(currentLabel != neiSampleIt.GetPixel(zminus))

			 + int(clsIdx != neiSampleIt.GetPixel(zplus)) 
			 - int(currentLabel != neiSampleIt.GetPixel(zplus));

		    // beta (sub level pairwise interation) term.
		    denergyLL = m_beta_z * denergyLL;


		    // The energy = - Log(likelihood).
		    denergyVMF = log (m_gmm.comp[clsIdx].sigma[subIdx]) + pow( (obsIt.Get() - m_gmm.comp[clsIdx].mu[subIdx]), 2) / ( 2 * pow(m_gmm.comp[clsIdx].sigma[subIdx], 2))
			 - log (m_gmm.comp[currentLabel].sigma[subIdx]) - pow( (obsIt.Get() - m_gmm.comp[currentLabel].mu[subIdx]), 2) / ( 2 * pow(m_gmm.comp[currentLabel].sigma[subIdx], 2));

		    denergy[clsIdx] = denergyPrior + denergyLL + denergyVMF;
	       }

	       newLabel = denergy.arg_min();
	       if (newLabel != currentLabel) {
		    neiSampleIt.SetCenterPixel( newLabel );
		    numFlips ++;
	       }

	  } // in mask
     } // iterators.
     return ( (float)(numFlips)/(float)(m_gmm.totalPts));
}


int GMModel::ICMGrp(ImageType3DChar::Pointer grpPtr,
		    ImageType3DChar::Pointer maskPtr,
		    std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec)
{
     vnl_vector<double> denergy(m_numClusters, 0);
     int currentLabel = 0;
     unsigned clsIdx = 0;

     //mask
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     // Define neighborhood iterator
     MyBoundCondType constCondition;
     constCondition.SetConstant(-1);
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType neiGrpIt(radius, grpPtr, grpPtr->GetLargestPossibleRegion() );
     neiGrpIt.OverrideBoundaryCondition(&constCondition);

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1}};


     // create histogram data block.
     itk::VariableLengthVector<unsigned short> histVector(m_numClusters);

     ImageType3DVecUS::Pointer histPtr = ImageType3DVecUS::New();
     ImageType3DVecUS::RegionType histRegion = grpPtr->GetLargestPossibleRegion();
     histPtr->SetRegions(histRegion);
     histPtr->SetNumberOfComponentsPerPixel(m_numClusters);
     histPtr->SetVectorLength(m_numClusters);
     histPtr->Allocate();

     histVector.Fill(0);
     histPtr->FillBuffer(histVector);
     IteratorType3DVecUS histIt(histPtr, histRegion);

     // compute histogram
     for (unsigned short subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	  ConstIteratorType3DChar sampleIt(sampleVec[subIdx][m_numSamples - 1], 
					   sampleVec[subIdx][m_numSamples - 1]->GetLargestPossibleRegion() );

	  for (sampleIt.GoToBegin(), histIt.GoToBegin(); 
	       !histIt.IsAtEnd(); ++ histIt, ++sampleIt) {
	       if (sampleIt.Get() >= 0) {
		    histVector = histIt.Get();
		    histVector[sampleIt.Get()] ++;
		    histIt.Set(histVector);
	       }
	  } // for sampelIt.
     } // for subIdx

     unsigned numFlips = 0; // number of voxed changed in this scan.
     unsigned newLabel = 0;     
     for (neiGrpIt.GoToBegin(), histIt.GoToBegin(), maskIt.GoToBegin(); 
	  !histIt.IsAtEnd(); ++ histIt, ++ neiGrpIt, ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       currentLabel = neiGrpIt.GetCenterPixel();
	       for (clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
		    denergy[clsIdx] = m_alpha * (histIt.Get()[currentLabel] - histIt.Get()[clsIdx])

			 + (int(clsIdx != neiGrpIt.GetPixel(xminus))
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
			    - int(currentLabel != neiGrpIt.GetPixel(zplus)) ) * m_beta_g;
	       }
	       
	       newLabel = denergy.arg_min();
	       if (newLabel != currentLabel) {
		    neiGrpIt.SetCenterPixel( newLabel );
		    numFlips ++;
	       }
	       
	  } // maskIt > 0
     } // neiGrpIt

     // printf("    group: numFlips is %i.\n", numFlips);
     return ((float)(numFlips) / (float)(m_gmm.totalPts) );
}

int GMModel::ICM(std::vector<ImageType3DChar::Pointer> &  grpVec,
		 std::vector<ImageType3DFloat::Pointer> & obsVec,
		 ImageType3DChar::Pointer maskPtr,
		 std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		 float tol) // max proportion of voxels changed for convergence.
{
     unsigned subIdx = 0;
     unsigned scanIdx = 0;
     
     // proportion of voxels changed for each subject and group. Last element
     // saves the value for group.
     vnl_vector<float> propFlip(m_numSubs + 1, 0); 

     do {

#pragma omp parallel for
	  for (subIdx = 0; subIdx < m_numSubs; subIdx ++) {
	       propFlip[subIdx] = ICMSub(grpVec[m_numSamples -1], obsVec[subIdx], maskPtr, sampleVec[subIdx][m_numSamples - 1], subIdx);
	  }
	  
	  propFlip[m_numSubs] = ICMGrp(grpVec[m_numSamples - 1], maskPtr, sampleVec);

	  if (m_verbose >= 1) {
	       printf("ICM(): scan %i done. proportion of voxels changed: %f.\n", scanIdx ++, propFlip.max_value());
	  }
	  scanIdx ++;

     } while (propFlip.max_value() > tol && scanIdx < 10);

     return 0;
}
 
float GMModel::GetTemperature()
{
     return m_temperature;
}

int GMModel::SetTemperature(float newTemp)
{
     m_temperature = newTemp;
     return 0;
}


// - sum_j log (P(X|Z)). energy function.
double GMModel::Q3(std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		      std::vector<ImageType3DFloat::Pointer> & obsVec,
		  ImageType3DChar::Pointer maskPtr)
{
     unsigned clsIdx = 0;

     //mask
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     double Q3 = 0;
     for (unsigned subIdx = 0; subIdx < m_numSubs; subIdx ++) {

	  std::vector<ImageType3DChar::Pointer> & sampleSubVec = sampleVec[subIdx];
	  IteratorType3DFloat obsIt(obsVec[subIdx], obsVec[subIdx]->GetLargestPossibleRegion().GetSize() );

	  for (unsigned short sampleIdx = 0; sampleIdx < m_numSamples; sampleIdx ++) {
	       ConstIteratorType3DChar sampleIt(sampleSubVec[sampleIdx], 
						sampleSubVec[sampleIdx]->GetLargestPossibleRegion() );
	       for (obsIt.GoToBegin(), maskIt.GoToBegin(), sampleIt.GoToBegin(); 
		    !maskIt.IsAtEnd(); 
		    ++ obsIt, ++ maskIt, ++ sampleIt) {
		    if (maskIt.Get() > 0) {
			 clsIdx = sampleIt.Get();
			 Q3 += ( log(m_gmm.comp[clsIdx].sigma[subIdx]) + pow( (obsIt.Get() - m_gmm.comp[clsIdx].mu[subIdx]), 2) / (2 * pow( (m_gmm.comp[clsIdx].sigma[subIdx]), 2) ) ) / m_numSamples;
		    }  // maskIt > 0
	       } // for maskIt

	  } // sampleIdx
     } // subIdx

     return (Q3);
}

