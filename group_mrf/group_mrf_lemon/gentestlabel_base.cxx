#include <common.h>
#include <utility.h>

using namespace lemon;
twister_base_gen_type mygenerator(42u);

int SamplingSub(ImageType3DChar::Pointer grpPtr,
		ImageType3DChar::Pointer subPtr,
		ImageType3DChar::Pointer maskPtr,
		ParStruct & par)

{
     typedef itk::ConstantBoundaryCondition< ImageType3DChar >  MyBoundCondType;
     typedef itk::NeighborhoodIterator< ImageType3DChar, MyBoundCondType > MyNeiItType;
     double denergy= 0, denergyPrior = 0, denergyLL=0;
     int cand;
     double p_acpt = 0; // probability to accept the candidate sample.
     signed int currentLabel = 0;
     unsigned grpLabel = 0;

     // Uniform integer generator.
     boost::uniform_int<> uni_int(0, par.numClusters - 1); // Uniform distribution of integers.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_int<> > roll_die(mygenerator, uni_int);

     // Uniform real random generator.
     boost::uniform_real<> uni_dist(0,1); // uniform distribution.
     boost::variate_generator<twister_base_gen_type&, boost::uniform_real<> > uni(mygenerator, uni_dist);


     //mask
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     // group.
     IteratorType3DChar grpIt(grpPtr, grpPtr->GetLargestPossibleRegion() );

     // Define neighborhood iterator
     MyBoundCondType constCondition;
     constCondition.SetConstant(-1);
     MyNeiItType::RadiusType radius;
     radius.Fill(1);



     MyNeiItType neiSampleIt( radius, subPtr, subPtr->GetRequestedRegion() );
     neiSampleIt.OverrideBoundaryCondition(&constCondition);

     MyNeiItType::OffsetType xplus = {{1,0, 0}};
     MyNeiItType::OffsetType xminus = {{-1, 0, 0}};
     MyNeiItType::OffsetType yplus = {{0, 1, 0}};
     MyNeiItType::OffsetType yminus = {{0, -1, 0}};
     MyNeiItType::OffsetType zplus = {{0, 0, 1}};
     MyNeiItType::OffsetType zminus = {{0, 0, -1}};



     // sampling Markov Random Field. 
     for (neiSampleIt.GoToBegin(), grpIt.GoToBegin(), maskIt.GoToBegin(); 
	  !neiSampleIt.IsAtEnd(); 
	  ++ neiSampleIt, ++ grpIt, ++ maskIt) {

	  if (maskIt.Get() > 0) {
	       currentLabel = neiSampleIt.GetCenterPixel();
	       cand = roll_die();
	       // alpha (group) term. Negative log of probabililty.
	       grpLabel = grpIt.Get();
	       denergyPrior = par.alpha * ( int(cand != grpLabel) - int(currentLabel != grpLabel) );
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
	       denergyLL = par.beta * denergyLL;

	       denergy = denergyPrior + denergyLL;

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


