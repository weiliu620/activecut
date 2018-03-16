#ifndef __SINGLEHEADER_H__
#define __SINGLEHEADER_H__

#include "GCoptimization.h"
#include<boost/bimap.hpp> 
#include <GCoptimization.h>


#ifdef GCO_ENERGYTYPE32
	typedef double EnergyType;        // 32-bit energy total
#else
	typedef double EnergyType;  // 64-bit energy total
#endif
	typedef float EnergyTermType;    // 32-bit energy terms
	typedef Energy<EnergyTermType,EnergyTermType,EnergyType> EnergyT;
	typedef EnergyT::Var VarID;
	typedef int LabelID;                     // Type for labels
	typedef VarID SiteID;                    // Type for sites
	typedef EnergyTermType (*SmoothCostFn)(SiteID s1, SiteID s2, LabelID l1, LabelID l2);
	typedef EnergyTermType (*DataCostFn)(SiteID s, LabelID l);
	typedef EnergyTermType (*SmoothCostFnExtra)(SiteID s1, SiteID s2, LabelID l1, LabelID l2,void *);
	typedef EnergyTermType (*DataCostFnExtra)(SiteID s, LabelID l,void *);


struct  SuperIdxType
{
     ImageType3DChar::IndexType idx;
     unsigned short subIdx;
};


struct SuperIdxTypeCompare {
     bool operator() (const SuperIdxType& lhs, const SuperIdxType& rhs) const
	  {
	       if (lhs.subIdx < rhs.subIdx) return true;
	       else if (lhs.subIdx > rhs.subIdx) return false;

	       else if (lhs.idx[2] < rhs.idx[2]) return true;
	       else if (lhs.idx[2] > rhs.idx[2]) return false;

	       else if (lhs.idx[1] < rhs.idx[1]) return true;
	       else if (lhs.idx[1] > rhs.idx[1]) return false;

	       else if (lhs.idx[0] < rhs.idx[0]) return true;
	       else if (lhs.idx[0] > rhs.idx[0]) return false;
	       else return false;

	  }
};

typedef boost::bimap<boost::bimaps::set_of<SuperIdxType, SuperIdxTypeCompare>, SiteID> BiMapType;
typedef BiMapType::value_type position;

#endif
