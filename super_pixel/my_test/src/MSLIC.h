// July 31, 2013
// header file for MSLIC.cxx
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;


void NormalizeVolumeData(
     int**&		          	        VolumeData,
     const int&                                 maxvalue,
     const int&                                width,
     const int&                                height,
     const int&                                depth);

void GetKValues_FiveChannelsAndXYZ(
     vector<double>&				kseedsT1,
     vector<double>&				kseedsT2,
     vector<double>&				kseedsGRE,
     vector<double>&				kseedsFLAIR,
     vector<double>&				kseedsSWI,
     vector<double>&				kseedsx,
     vector<double>&				kseedsy,
     vector<double>&				kseedsz,
     int**&		          	        T1vecvec,
     int**&		          	        T2vecvec,
     int**&		          	        GREvecvec,
     int**&		          	        FLAIRvecvec,
     int**&		          	        SWIvecvec,
     const int&				STEP,
     const int&                                width,
     const int&                                height,
     const int&                                depth);



void PerformSupervoxelMSLIC(
	vector<double>&				kseedsT1,
	vector<double>&				kseedsT2,
	vector<double>&				kseedsGRE,
	vector<double>&				kseedsFLAIR,
	vector<double>&				kseedsSWI,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	vector<double>&				kseedsz,
        int**&					klabels,
	int**&		          	        T1vecvec,
	int**&		          	        T2vecvec,
	int**&		          	        GREvecvec,
	int**&		          	        FLAIRvecvec,
	int**&		          	        SWIvecvec,
        const int&				STEP,
	const double&				compactness,
	const int&					width,
	const int&					height,
	const int&					depth);



void EnforceSupervoxelLabelConnectivity(
	int**&						labels,//input - previous labels, output - new labels
	const int&					width,
	const int&					height,
	const int&					depth,
	int&						numlabels,
	const int&					STEP);


void DoSupervoxelSegmentation(
     int**&				ubuff1,
     int**&				ubuff2,
     int**&				ubuff3,
     int**&				ubuff4,
     int**&				ubuff5,
     const int&					width,
     const int&					height,
     const int&					depth,
     int**&						klabels,
     int&						numlabels,
     const int&					supervoxelsize,
     const double&               compactness);
