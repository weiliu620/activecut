#include <common_mcem.h>
#include <gmmodel_mcem.h>
#include <gmmodel_graph.h>
#include <GCoptimization.h>

SiteID FindSiteID(ImageType3DChar::IndexType maskIdx,
		  unsigned subIdx,
		  const BiMapType & biMap);

template < typename T >
T **Allocate2DArray( int nRows, int nCols)
{
    T **ppi;
    T *pool;
    T *curPtr;
    //(step 1) allocate memory for array of elements of column

    ppi = new T*[nRows];

    //(step 2) allocate memory for array of elements of each row
    pool = new T [nRows * nCols];

    // Now point the pointers in the right place
    curPtr = pool;
    for( int i = 0; i < nRows; i++)
    {
        *(ppi + i) = curPtr;
         curPtr += nCols;
    }
    return ppi;
}

template < typename T >
void Free2DArray(T** Array)
{
    delete [] *Array;
    delete [] Array;
}

int GMModel::InitGraph(std::vector<ImageType3DFloat::Pointer> & obsVec,
		 std::vector<ImageType3DChar::Pointer> &  grpVec,
		 std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		 ImageType3DChar::Pointer maskPtr)
{
     // build a big graph including all subjects and group lebel.
     unsigned long numSites = m_gmm.totalPts * (m_numSubs + 1);
     m_gc = new GCoptimizationGeneralGraph(numSites, m_numClusters);

     // build a numNeighbors.
     SiteID * numNeighbors = new SiteID [numSites];
     

     // Create a neighborsIndexes. max number of neighbors happens on a voxel at
     // group level, where the neighbors include all 6 group voxels and all
     // subjects voxels at this location.
     unsigned maxNumNeighbors = m_numSubs + 6;
     SiteID ** neighborsIndexes = Allocate2DArray<SiteID>(numSites, maxNumNeighbors);

     // Create a neighborsWeights.
     EnergyTermType ** neighborsWeights = Allocate2DArray<EnergyTermType>(numSites, maxNumNeighbors);

     // assign values
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());
     ImageType3DChar::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     ImageType3DChar::IndexType maskIdx;     


     // creat a boost::bimap for the mapping between super index (x,y,z,subIdx)
     // and site ID.  now iterate the image to fill the biMap.
     SiteID siteID = 0;

     // we think of group map as (J+1)th subject, and put it at the end of all
     // other subjects. So if we found (x,y,z,j) has a j that larger than total
     // number of subjects, that must be a voxel in group.
     SuperIdxType superIdx;
     for (int subIdx = 0; subIdx < m_numSubs + 1; subIdx ++) {
     	  for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
     	       if (maskIt.Get() > 0) {
		    superIdx.idx = maskIt.GetIndex();
		    superIdx.subIdx = subIdx;
		    m_biMap.insert(position(superIdx, siteID));

		    // if (siteID < 100 ) {
		    // 	 std::cout << "siteID = " << siteID <<". size of the bimap is: " << m_biMap.size();

		    // 	 printf("  [%i][%i][%i],[%i] <----> site %i\n", superIdx.idx[0], superIdx.idx[1], superIdx.idx[2], superIdx.subIdx, siteID );
		    // }

     		    siteID ++;		    
     	       }  // maskIt > 0
     	  } // maskIt
     } // subIdx 


     printf("last site is %i\n", m_biMap.left.at(superIdx));

//      for (BiMapType::const_iterator i = m_biMap.begin(), iend = m_biMap.end(); i != iend; ++i) {
// 	  printf("[%i][%i][%i],[%i] <----> site %i\n", i->left.idx[0], i->left.idx[1], i->left.idx[2], i->left.subIdx, i->right );
//      }

//      std::cout << "size of the bimap is: " << m_biMap.size() << std::endl;
//      for (SiteID tmpsite = 2000; tmpsite < 2010; tmpsite ++) {
// 	  printf("site %i ---> [%i][%i][%i],[%i]\n", tmpsite, m_biMap.right.at(tmpsite).idx[0], m_biMap.right.at(tmpsite).idx[1], m_biMap.right.at(tmpsite).idx[2], m_biMap.right.at(tmpsite).subIdx);
// }
     


     typedef itk::ConstantBoundaryCondition<ImageType3DChar>  BoundaryConditionType;
     typedef itk::NeighborhoodIterator< ImageType3DChar, BoundaryConditionType> NeighborhoodIteratorType;
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType neiMaskIt( radius, maskPtr, 
					   maskPtr->GetLargestPossibleRegion() );

     NeighborhoodIteratorType::OffsetType xplus = {{1,0, 0}};
     NeighborhoodIteratorType::OffsetType xminus = {{-1, 0, 0}};
     NeighborhoodIteratorType::OffsetType yplus = {{0, 1, 0}};
     NeighborhoodIteratorType::OffsetType yminus = {{0, -1, 0}};
     NeighborhoodIteratorType::OffsetType zplus = {{0, 0, 1}};
     NeighborhoodIteratorType::OffsetType zminus = {{0, 0, -1}};
     
     SiteID curSite = 0, neiSite = 0;

     // This part does three things. 1) add edge within subjects (like z_r and
     // z_s), for all subjects. This is when subIdx = 0 to m_numSubs - 1. 2) add
     // edge within groups (like g_r and g_s), which happens when subIdx =
     // m_numSubs, since group label map, by my definition in this code, is the
     // j+1 th subject. 3) add edge between group and subjects.

     float beta = 0; // save either beta_g or beta_z, depending on value of
		     // subIdx.
     for (int subIdx = 0; subIdx < m_numSubs + 1; subIdx ++) {
	  if (subIdx < m_numSubs) beta = m_beta_z;
	  else beta = m_beta_g;

	  for (neiMaskIt.GoToBegin(); !neiMaskIt.IsAtEnd(); ++ neiMaskIt) {
	       if (neiMaskIt.GetCenterPixel() > 0) {

		    curSite = FindSiteID(neiMaskIt.GetIndex(), subIdx, m_biMap);

		    // if curSite is in subject instead of in group map, add the
		    // edge between group and sub. Do this for both subject and
		    // group labels.
		    if (subIdx < m_numSubs + 1) {
			 // Get group label.
			 neiSite = FindSiteID(neiMaskIt.GetIndex(), m_numSubs, m_biMap);
			 neighborsIndexes[curSite][ numNeighbors[curSite] ] = neiSite;
			 neighborsWeights[curSite][ numNeighbors[curSite] ] = m_alpha;
			 numNeighbors[curSite] ++;

			 // also add neighbors for group label.
			 neighborsIndexes[neiSite][ numNeighbors[neiSite] ] = curSite;
			 neighborsWeights[neiSite][ numNeighbors[neiSite] ] = m_alpha;
			 numNeighbors[neiSite] ++;
		    }
			 

		    if (neiMaskIt.GetPixel(xminus) > 0) {
			 neiSite = FindSiteID(neiMaskIt.GetIndex(xminus), subIdx, m_biMap);
			 neighborsIndexes[curSite][ numNeighbors[curSite] ] = neiSite;
			 neighborsWeights[curSite][ numNeighbors[curSite] ] = beta;
			 numNeighbors[curSite] ++;
		    }

		    if (neiMaskIt.GetPixel(xplus) > 0) {
			 neiSite = FindSiteID(neiMaskIt.GetIndex(xplus), subIdx, m_biMap);
			 neighborsIndexes[curSite][ numNeighbors[curSite] ] = neiSite;
			 neighborsWeights[curSite][ numNeighbors[curSite] ] = beta;
			 numNeighbors[curSite] ++;
		    }
		    
		    if (neiMaskIt.GetPixel(yminus) > 0) {
			 neiSite = FindSiteID(neiMaskIt.GetIndex(yminus), subIdx, m_biMap);
			 neighborsIndexes[curSite][ numNeighbors[curSite] ] = neiSite;
			 neighborsWeights[curSite][ numNeighbors[curSite] ] = beta;
			 numNeighbors[curSite] ++;
		    }

		    if (neiMaskIt.GetPixel(yplus) > 0) {
			 neiSite = FindSiteID(neiMaskIt.GetIndex(yplus), subIdx, m_biMap);
			 neighborsIndexes[curSite][ numNeighbors[curSite] ] = neiSite;
			 neighborsWeights[curSite][ numNeighbors[curSite] ] = beta;
			 numNeighbors[curSite] ++;
		    }


		    if (neiMaskIt.GetPixel(zminus) > 0) {
			 neiSite = FindSiteID(neiMaskIt.GetIndex(zminus), subIdx, m_biMap);
			 neighborsIndexes[curSite][numNeighbors[curSite] ] = neiSite;
			 neighborsWeights[curSite][ numNeighbors[curSite] ] = beta;
			 numNeighbors[curSite] ++;
		    }

		    if (neiMaskIt.GetPixel(zplus) > 0) {
			 neiSite = FindSiteID(neiMaskIt.GetIndex(zplus), subIdx, m_biMap);
			 neighborsIndexes[curSite][ numNeighbors[curSite] ] = neiSite;
			 neighborsWeights[curSite][ numNeighbors[curSite] ] = beta;
			 numNeighbors[curSite] ++;
		    }

	       }  // maskIt > 0
	  } // maskIt

     } // subIdx
     m_gc->setAllNeighbors(numNeighbors, neighborsIndexes, neighborsWeights);
}


// Given volume index (x,y,z) and subject ID (may be real subject or can be
// group map which has ID m_numSubs) and mapping ptr, compute the site ID of
// this point.
SiteID FindSiteID(ImageType3DChar::IndexType maskIdx,
		  unsigned subIdx,
		  const BiMapType & biMap)
{
     SuperIdxType superIdx;
     superIdx.idx = maskIdx;
     superIdx.subIdx = subIdx;
     return biMap.left.at(superIdx);
}


int GMModel::graphcuts(std::vector<ImageType3DFloat::Pointer> & obsVec,
		 std::vector<ImageType3DChar::Pointer> &  grpVec,
		 std::vector< std::vector<ImageType3DChar::Pointer> > & sampleVec,
		 ImageType3DChar::Pointer maskPtr)

{
     unsigned long numSites = m_gmm.totalPts * (m_numSubs + 1);
     EnergyTermType * dataCostArray = new EnergyTermType[numSites * m_numClusters];
     SuperIdxType superIdx;
     SiteID  siteID = 0;

     // assign data cost.
     float obsValue = 0;
     for (siteID = 0; siteID < numSites; siteID ++) {
     	  for (unsigned clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
     	       superIdx = m_biMap.right.at(siteID);
     	       if (superIdx.subIdx < m_numSubs) {
     		    obsValue = obsVec[superIdx.subIdx]->GetPixel(superIdx.idx);
     		    dataCostArray[siteID * m_numClusters + clsIdx] = log (m_gmm.comp[clsIdx].sigma[superIdx.subIdx]) +  pow((obsValue - m_gmm.comp[clsIdx].mu[superIdx.subIdx]), 2) / (2 * pow(m_gmm.comp[clsIdx].sigma[superIdx.subIdx], 2));

     	       }
     	       else {
     		    // this siteID is group label. Since group label is not
     		    // connected to observed data, the data cost is set to a
     		    // small number.
     		    dataCostArray[siteID * m_numClusters + clsIdx] = 0;
     	       }
     	  } // clsIdx
     }// siteID

     m_gc -> setDataCost(dataCostArray);


     // assign smoothness cost.
     EnergyTermType * smoothCostArray = new EnergyTermType [m_numClusters * m_numClusters];
     for (unsigned label1 = 0; label1 < m_numClusters; label1++) {
	  for (unsigned label2 = 0; label2 < m_numClusters; label2++) {
	       smoothCostArray[label1 + m_numClusters * label2] = (float) (label1 != label2);
	  }
     }
     m_gc->setSmoothCost(smoothCostArray);
     
     // init the graph.
     ImageType3D::IndexType idx;
     unsigned clsIdx = 0;
     for (siteID = 0; siteID < numSites; siteID ++) {
	  superIdx = m_biMap.right.at(siteID);
	  if (superIdx.subIdx < m_numSubs) {
	       clsIdx = sampleVec[superIdx.subIdx][m_numSamples - 1]->GetPixel(superIdx.idx);
	  }
	  else {
	       // must be in group.
	       clsIdx = grpVec[m_numSamples - 1]->GetPixel(superIdx.idx);
	  }
	  m_gc->setLabel(siteID, clsIdx);
     }

     // do the optimization work.
     for (unsigned clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_gc->alpha_expansion(clsIdx);
	  
	  if (m_verbose >= 1) {
	       printf("Total energy = %f, data energy = %f, smooth energy = %f\n", m_gc->compute_energy(), m_gc->giveDataEnergy(), m_gc->giveSmoothEnergy() );
	  }
     }

     // do the optimization work.
     for (unsigned clsIdx = 0; clsIdx < m_numClusters; clsIdx ++) {
	  m_gc->alpha_expansion(clsIdx);
	  
	  if (m_verbose >= 1) {
	       printf("Total energy = %f, data energy = %f, smooth energy = %f\n", m_gc->compute_energy(), m_gc->giveDataEnergy(), m_gc->giveSmoothEnergy() );
	  }
     }

     
     // write labels back to images.
     for (siteID = 0; siteID < numSites; siteID ++) {
	  superIdx = m_biMap.right.at(siteID);
	  if (superIdx.subIdx < m_numSubs) {
	       sampleVec[superIdx.subIdx][m_numSamples - 1]->SetPixel(superIdx.idx, m_gc->whatLabel(siteID));
	  }
	  else {
	       // this siteID is group label.
	       grpVec[m_numSamples - 1]->SetPixel(superIdx.idx, m_gc->whatLabel(siteID));
	  }
     }// siteID

     delete [] dataCostArray;
     delete [] smoothCostArray;



}
