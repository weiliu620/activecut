//
// C++ Interface: PixelCorr
//
// Description: 
//
//
// Author:  <Wei Liu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PixelCorr_H
#define PixelCorr_H
#endif /* PixelCorr_H */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <itkAnalyzeImageIO.h>
#include <itkNiftiImageIO.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#ifndef SaveCorr_H
#define SaveCorr_H
#endif /* PixelCorr_H */

#ifndef SaveHighlight_H
#define SaveHighlight_H
#endif

#ifndef M_PI	
#define M_PI 3.1415926535897932385	
#endif

#include <boost/math/distributions/normal.hpp> // for normal_distribution.
#include <iterator>
#include <algorithm>
#include <iostream>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <stdlib.h>
#include <time.h>

double GetNormalSample(short * ranArray, int size);
//double GetNormalSample(int maxRand);




	


