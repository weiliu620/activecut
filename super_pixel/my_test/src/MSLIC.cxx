// July 31, 2013
// c++ code for multimodal super pixel

#include <cfloat>
#include <math.h>       /* pow */
#include <MSLIC.h>

//===========================================================================
///	GetKValues_FiveChannelsAndXYZ
///
/// The k seed values are taken as uniform spatial pixel samples.
//===========================================================================

void NormalizeVolumeData(
     int**&		          	        VolumeData,
     const int&                                 maxvalue,
     const int&                                width,
     const int&                                height,
     const int&                                depth)
{
     int i,j,k,s;
     int MaxIntensity = 0;
     int MinIntensity = 1000000;

     int SliceSize = width*height;
     //int** NormlizedVolume;
     //NormlizedVolume = new int*[depth];

     // find the maximum & minimum intensity values
     for(k=0; k<depth; k++)
     {
	  for(j=0; j< height; j++)
	  {
	       for(i=0; i< width; i++)
	       {
		    int CurIndex = j*width + i;
		    if(VolumeData[k][CurIndex] > MaxIntensity)
		    {
			 MaxIntensity = VolumeData[k][CurIndex];
		    }
		    if(VolumeData[k][CurIndex] < MinIntensity)
		    {
			 MinIntensity = VolumeData[k][CurIndex];
		    }
	       }
	  }
     }
     std::cout << "Before normlization, Max value: " << MaxIntensity << std::endl;
     std::cout << "Before normlization, Min value: " << MinIntensity << std::endl;

     // normalize the data
     float DiffMaxMin = MaxIntensity - MinIntensity;
     for(k=0; k<depth; k++)
     {
	  for(j=0; j< height; j++)
	  {
	       for(i=0; i< width; i++)
	       {
		    int CurIndex = j*width + i;
		    float CurNewValue = ((VolumeData[k][CurIndex] - MinIntensity)/DiffMaxMin)*maxvalue;
		    VolumeData[k][CurIndex] = (int) CurNewValue;
	       }
	  }
     }

     MaxIntensity = 0;
     MinIntensity = 1000000;
     for(k=0; k<depth; k++)
     {
	  //NormlizedVolume[k] = new int[SliceSize];
	  for(j=0; j< height; j++)
	  {
	       for(i=0; i< width; i++)
	       {
		    int CurIndex = j*width + i;
		    if(VolumeData[k][CurIndex] > MaxIntensity)
		    {
			 MaxIntensity = VolumeData[k][CurIndex];
		    }
		    if(VolumeData[k][CurIndex] < MinIntensity)
		    {
			 MinIntensity = VolumeData[k][CurIndex];
		    }
		    
	       }
	  }
     }
     
     std::cout << "After normlization, Max value: " << MaxIntensity << std::endl;
     std::cout << "After normlization, Min value: " << MinIntensity << std::endl;
}

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
     const int&                                depth)
{

     const bool hexgrid = false;
     int numseeds(0);
     int n(0);

     int xstrips = (0.5+double(width)/double(STEP));
     int ystrips = (0.5+double(height)/double(STEP));
     int zstrips = (0.5+double(depth)/double(STEP));

     int xerr = width  - STEP*xstrips;if(xerr < 0){xstrips--;xerr = width - STEP*xstrips;}
     int yerr = height - STEP*ystrips;if(yerr < 0){ystrips--;yerr = height- STEP*ystrips;}
     int zerr = depth  - STEP*zstrips;if(zerr < 0){zstrips--;zerr = depth - STEP*zstrips;}

     double xerrperstrip = double(xerr)/double(xstrips);
     double yerrperstrip = double(yerr)/double(ystrips);
     double zerrperstrip = double(zerr)/double(zstrips);

     int xoff = STEP/2;
     int yoff = STEP/2;
     int zoff = STEP/2;

     //-------------------------
     numseeds = xstrips*ystrips*zstrips;
     //-------------------------
     kseedsT1.resize(numseeds);
     kseedsT2.resize(numseeds);
     kseedsGRE.resize(numseeds);
     kseedsFLAIR.resize(numseeds);
     kseedsSWI.resize(numseeds);
     kseedsx.resize(numseeds);
     kseedsy.resize(numseeds);
     kseedsz.resize(numseeds);

     for( int z = 0; z < zstrips; z++ )
     {
	  int ze = z*zerrperstrip;
	  int d = (z*STEP+zoff+ze);
	  for( int y = 0; y < ystrips; y++ )
	  {
	       int ye = y*yerrperstrip;
	       for( int x = 0; x < xstrips; x++ )
	       {
		    int xe = x*xerrperstrip;
		    int i = (y*STEP+yoff+ye)*width + (x*STEP+xoff+xe);
				
		    kseedsT1[n] = T1vecvec[d][i];
		    kseedsT2[n] = T2vecvec[d][i];
		    kseedsGRE[n] = GREvecvec[d][i];
		    kseedsFLAIR[n] = FLAIRvecvec[d][i];
		    kseedsSWI[n] = SWIvecvec[d][i];
		    kseedsx[n] = (x*STEP+xoff+xe);
		    kseedsy[n] = (y*STEP+yoff+ye);
		    kseedsz[n] = d;
		    n++;
	       }
	  }
     }

}


//===========================================================================
///	PerformSupervoxelMSLIC
///
///	Performs k mean segmentation. It is fast because it searches locally, not
/// over the entire image.
//===========================================================================
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
	const int&					depth)
{
	int sz = width*height;
	const int numk = kseedsT1.size();
        //int numitr(0);

	//----------------
	int offset = STEP;
        //if(STEP < 8) offset = STEP*1.5;//to prevent a crash due to a very small step size
	//----------------

	vector<double> clustersize(numk, 0);
	vector<double> inv(numk, 0);//to store 1/clustersize[k] values

	vector<double> sigmaT1(numk, 0);
	vector<double> sigmaT2(numk, 0);
	vector<double> sigmaGRE(numk, 0);
	vector<double> sigmaFLAIR(numk, 0);
	vector<double> sigmaSWI(numk, 0);
	vector<double> sigmax(numk, 0);
	vector<double> sigmay(numk, 0);
	vector<double> sigmaz(numk, 0);

	vector< double > initdouble(sz, DBL_MAX);
	vector< vector<double> > distvec(depth, initdouble);
	//vector<double> distvec(sz, DBL_MAX);

	double invwt = 1.0/((STEP/compactness)*(STEP/compactness));//compactness = 20.0 is usually good.

	int x1, y1, x2, y2, z1, z2;
	double valueT1, valueT2, valueGRE, valueFLAIR, valueSWI;
	double dist;
	double distxyz;
	for( int itr = 0; itr < 5; itr++ )
	{
		distvec.assign(depth, initdouble);
		for( int n = 0; n < numk; n++ )
		{
                        y1 = max(0.0,			kseedsy[n]-offset);
                        y2 = min((double)height,	kseedsy[n]+offset);
                        x1 = max(0.0,			kseedsx[n]-offset);
                        x2 = min((double)width,	kseedsx[n]+offset);
                        z1 = max(0.0,			kseedsz[n]-offset);
                        z2 = min((double)depth,	kseedsz[n]+offset);


			for( int z = z1; z < z2; z++ )
			{
				for( int y = y1; y < y2; y++ )
				{
					for( int x = x1; x < x2; x++ )
					{
						int i = y*width + x;

						valueT1 = T1vecvec[z][i];
						valueT2 = T2vecvec[z][i];
						valueGRE = GREvecvec[z][i];
						valueFLAIR = FLAIRvecvec[z][i];
						valueSWI = SWIvecvec[z][i];

						dist =			(valueT1 - kseedsT1[n])*(valueT1 - kseedsT1[n]) +
										(valueT2 - kseedsT2[n])*(valueT2 - kseedsT2[n]) +
										(valueGRE - kseedsGRE[n])*(valueGRE - kseedsGRE[n]) +
						                                (valueFLAIR - kseedsFLAIR[n])*(valueFLAIR - kseedsFLAIR[n]) + 
						                                (valueSWI - kseedsSWI[n])*(valueSWI - kseedsSWI[n]);

						distxyz =		(x - kseedsx[n])*(x - kseedsx[n]) +
										(y - kseedsy[n])*(y - kseedsy[n]) +
										(z - kseedsz[n])*(z - kseedsz[n]);
						//------------------------------------------------------------------------
						dist += distxyz*invwt;
						//------------------------------------------------------------------------
						if( dist < distvec[z][i] )
						{
							distvec[z][i] = dist;
							klabels[z][i]  = n;
						}
					}
				}
			}
		}
		//-----------------------------------------------------------------
		// Recalculate the centroid and store in the seed values
		//-----------------------------------------------------------------
		//instead of reassigning memory on each iteration, just reset.
	
		sigmaT1.assign(numk, 0);
		sigmaT2.assign(numk, 0);
		sigmaGRE.assign(numk, 0);
		sigmaFLAIR.assign(numk, 0);
		sigmaSWI.assign(numk, 0);
		sigmax.assign(numk, 0);
		sigmay.assign(numk, 0);
		sigmaz.assign(numk, 0);
		clustersize.assign(numk, 0);

		for( int d = 0; d < depth; d++  )
		{
			int ind(0);
			for( int r = 0; r < height; r++ )
			{
				for( int c = 0; c < width; c++ )
				{
					sigmaT1[klabels[d][ind]] += T1vecvec[d][ind];
					sigmaT2[klabels[d][ind]] += T2vecvec[d][ind];
					sigmaGRE[klabels[d][ind]] += GREvecvec[d][ind];
					sigmaFLAIR[klabels[d][ind]] += FLAIRvecvec[d][ind];
					sigmaSWI[klabels[d][ind]] += SWIvecvec[d][ind];
					sigmax[klabels[d][ind]] += c;
					sigmay[klabels[d][ind]] += r;
					sigmaz[klabels[d][ind]] += d;

					clustersize[klabels[d][ind]] += 1.0;
					ind++;
				}
			}
		}

		{for( int k = 0; k < numk; k++ )
		{
			if( clustersize[k] <= 0 ) clustersize[k] = 1;
			inv[k] = 1.0/clustersize[k];//computing inverse now to multiply, than divide later
		}}
		
		{for( int k = 0; k < numk; k++ )
		{
			kseedsT1[k] = sigmaT1[k]*inv[k];
			kseedsT2[k] = sigmaT2[k]*inv[k];
			kseedsGRE[k] = sigmaGRE[k]*inv[k];
			kseedsFLAIR[k] = sigmaFLAIR[k]*inv[k];
			kseedsSWI[k] = sigmaSWI[k]*inv[k];
			kseedsx[k] = sigmax[k]*inv[k];
			kseedsy[k] = sigmay[k]*inv[k];
			kseedsz[k] = sigmaz[k]*inv[k];
		}}
	}
}


//===========================================================================
///	RelabelStraySupervoxels
//===========================================================================
void EnforceSupervoxelLabelConnectivity(
	int**&						labels,//input - previous labels, output - new labels
	const int&					width,
	const int&					height,
	const int&					depth,
	int&						numlabels,
	const int&					STEP)
{
	const int dx10[10] = {-1,  0,  1,  0, -1,  1,  1, -1,  0, 0};
	const int dy10[10] = { 0, -1,  0,  1, -1, -1,  1,  1,  0, 0};
	const int dz10[10] = { 0,  0,  0,  0,  0,  0,  0,  0, -1, 1};

	int sz = width*height;
	const int SUPSZ = STEP*STEP*STEP;

	int adjlabel(0);//adjacent label
        int* xvec = new int[SUPSZ*10];//a large enough size
        int* yvec = new int[SUPSZ*10];//a large enough size
        int* zvec = new int[SUPSZ*10];//a large enough size
	//------------------
	// memory allocation
	//------------------
	int** nlabels = new int*[depth];
	{for( int d = 0; d < depth; d++ )
	{
		nlabels[d] = new int[sz];
		for( int i = 0; i < sz; i++ ) nlabels[d][i] = -1;
	}}
	//------------------
	// labeling
	//------------------
	int lab(0);
	{for( int d = 0; d < depth; d++ )
	{
		int i(0);
		for( int h = 0; h < height; h++ )
		{
			for( int w = 0; w < width; w++ )
			{
				if(nlabels[d][i] < 0)
				{
					nlabels[d][i] = lab;
					//-------------------------------------------------------
					// Quickly find an adjacent label for use later if needed
					//-------------------------------------------------------
					{for( int n = 0; n < 10; n++ )
					{
						int x = w + dx10[n];
						int y = h + dy10[n];
						int z = d + dz10[n];
						if( (x >= 0 && x < width) && (y >= 0 && y < height) && (z >= 0 && z < depth) )
						{
							int nindex = y*width + x;
							if(nlabels[z][nindex] >= 0)
							{
								adjlabel = nlabels[z][nindex];
							}
						}
					}}
					
					xvec[0] = w; yvec[0] = h; zvec[0] = d;
					int count(1);
					for( int c = 0; c < count; c++ )
					{
						for( int n = 0; n < 10; n++ )
						{
							int x = xvec[c] + dx10[n];
							int y = yvec[c] + dy10[n];
							int z = zvec[c] + dz10[n];

							if( (x >= 0 && x < width) && (y >= 0 && y < height) && (z >= 0 && z < depth))
							{
								int nindex = y*width + x;

								if( 0 > nlabels[z][nindex] && labels[d][i] == labels[z][nindex] )
								{
									xvec[count] = x;
									yvec[count] = y;
									zvec[count] = z;
									nlabels[z][nindex] = lab;
									count++;
								}
							}

						}
					}
					//-------------------------------------------------------
					// If segment size is less then a limit, assign an
					// adjacent label found before, and decrement label count.
					//-------------------------------------------------------
					if(count <= (SUPSZ >> 2))//this threshold can be changed according to needs
					{
						for( int c = 0; c < count; c++ )
						{
							int ind = yvec[c]*width+xvec[c];
							nlabels[zvec[c]][ind] = adjlabel;
						}
						lab--;
					}
					//--------------------------------------------------------
					lab++;
				}
				i++;
			}
		}
	}}
	//------------------
	// mem de-allocation
	//------------------
	{for( int d = 0; d < depth; d++ )
	{
		for( int i = 0; i < sz; i++ ) labels[d][i] = nlabels[d][i];
	}}
	{for( int d = 0; d < depth; d++ )
	{
		delete [] nlabels[d];
	}}
	delete [] nlabels;
	//------------------
	if(xvec) delete [] xvec;
	if(yvec) delete [] yvec;
	if(zvec) delete [] zvec;
	//------------------
	numlabels = lab;
	//------------------
}

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
     const double&               compactness)
{
     const int STEP = 0.5 + pow(double(supervoxelsize),1.0/3.0);
     std::cout << "STEP is: " << STEP << std::endl;

     vector<double> kseedsT1(0);
     vector<double> kseedsT2(0);
     vector<double> kseedsGRE(0);
     vector<double> kseedsFLAIR(0);
     vector<double> kseedsSWI(0);
     vector<double> kseedsx(0);
     vector<double> kseedsy(0);
     vector<double> kseedsz(0);

     //--------------------------------------------------
     int sz = width*height;

     //--------------------------------------------------
     int** T1vecvec = new int*[depth];
     int** T2vecvec = new int*[depth];
     int** GREvecvec = new int*[depth];
     int** FLAIRvecvec = new int*[depth];
     int** SWIvecvec = new int*[depth];
     for( int d = 0; d < depth; d++ )
     {
	  T1vecvec[d] = new int[sz];
	  T2vecvec[d] = new int[sz];
	  GREvecvec[d] = new int[sz];
	  FLAIRvecvec[d] = new int[sz];
	  SWIvecvec[d] = new int[sz];
	  for( int s = 0; s < sz; s++ )
	  {
	       klabels[d][s] = -1;
	  }
     }

     T1vecvec = ubuff1;
     T2vecvec = ubuff2;
     GREvecvec = ubuff3;
     FLAIRvecvec = ubuff4;
     SWIvecvec = ubuff5;

     std::cout << "T1 normalization " << std::endl;
     NormalizeVolumeData(T1vecvec, 255, width, height, depth);
     std::cout << "T2 normalization " << std::endl;
     NormalizeVolumeData(T2vecvec, 255, width, height, depth);
     std::cout << "GRE normalization " << std::endl;
     NormalizeVolumeData(GREvecvec, 255, width, height, depth);
     std::cout << "FLAIR normalization " << std::endl;
     NormalizeVolumeData(FLAIRvecvec, 255, width, height, depth);
     std::cout << "SWI normalization " << std::endl;
     NormalizeVolumeData(SWIvecvec, 255, width, height, depth);

     // get intial postion of seeds
     std::cout << "Step1: Get intial postions of each patch. " << std::endl;
     GetKValues_FiveChannelsAndXYZ(kseedsT1, kseedsT2, kseedsGRE, kseedsFLAIR, kseedsSWI, kseedsx, kseedsy, kseedsz, T1vecvec, T2vecvec, GREvecvec, FLAIRvecvec, SWIvecvec, STEP, width, height, depth);

     // do supervoxel computation
     std::cout << "Step2: Do super voxel computation. " << std::endl;
     PerformSupervoxelMSLIC(kseedsT1, kseedsT2, kseedsGRE, kseedsFLAIR, kseedsSWI, kseedsx, kseedsy, kseedsz, klabels, T1vecvec, T2vecvec, GREvecvec, FLAIRvecvec, SWIvecvec, STEP, compactness, width, height, depth);

     // enforce supervoxel connectivity
     std::cout << "Step3: enforce super voxel connectivity. " << std::endl;
     EnforceSupervoxelLabelConnectivity(klabels, width, height, depth, numlabels, STEP);

     std::cout << "Super voxel computation is done! " << std::endl;
}
