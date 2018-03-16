unsigned ReadObsToVec(std::string fmripath,
		      std::vector<VnlVectorImageType::Pointer> & fmriVec,
		      ImageType3DChar::Pointer maskPtr)
{
     // some of the code below refers to boost filesystem tutorial tut04.cpp.
     IteratorType3DChar maskIt(maskPtr, maskPtr->GetLargestPossibleRegion() );

     // default constructed iterator acts as a end iterator.
     boost::filesystem::path fmripathVar(fmripath);
     boost::filesystem::directory_iterator fmripathEnd;

     // prepare for sorting directory entries.
     typedef std::vector<boost::filesystem::path> PathVec;
     PathVec  sortedEntries;
     copy(boost::filesystem::directory_iterator(fmripath), boost::filesystem::directory_iterator(), std::back_inserter(sortedEntries) );
     sort(sortedEntries.begin(), sortedEntries.end() );

     // Get number of subjects.
     unsigned numSubs  = sortedEntries.size();
     fmriVec.resize(numSubs);

     // get time series length.
     boost::filesystem::directory_iterator obspathIt(obspathVar);
     ReaderType4DFloat srcReader = ReaderType4DFloat::New();
     srcReader->SetFileName( (*obspathIt).path().string() );
     srcReader->Update();
     ImageType4DF::SizeType srcSize = srcReader->GetOutput()->GetLargestPossibleRegion().GetSize();
     ImageType4DF::IndexType srcIdx;

     VnlVectorImageType::SizeType destSize = maskPtr->GetLargestPossibleRegion().GetSize();
     VnlVectorImageType::IndexType destIdx;
     destIdx.Fill(0);
     VnlVectorImageType::RegionType destRegion;
     destRegion.SetSize(destSize);
     destRegion.SetIndex(destIdx);


     unsigned subIdx = 0;
     vnl_vector<float> timeSeries( srcSize[3], 0 );
     vnl_vector<float> zeroVec( srcSize[3], 0 );
     for (PathVec::const_iterator srcpathIt(sortedEntries.begin() ); srcpathIt != sortedEntries.end(); ++ srcpathIt) {

	  // Source.
	  srcReader->SetFileName( (*obspathIt).string() );
	  srcReader->Update();
	  srcPtr = srcReader->GetOutput();

	  // Destination.
	  fmriVec[subIdx] = VnlVectorImageType::New();
	  fmriVec[subIdx]->SetRegions( destRegion );
	  fmriVec[subIdx]->SetNumberOfComponentsPerPixel(srcSize[3]);
	  fmriVec[subIdx]->SetVectorLength(srcSize[3]);
	  fmriVec[subIdx]->Allocate();
	  fmriVec[subIdx]->FillBuffer ( zeroVec );

	  // transfer data.
	  IteratorType4DF destIt(fmriVec[subIdx], fmriVec[subIdx]->GetLargestPossibleRegion() );

	  for (maskIt.GoToBegin(), destIt.GoToBegin(); !destIt.IsAtEnd(); ++ maskIt, ++ destIt) {
	       if (maskIt.Get() > 0) {
		    destdx = fmriIt.GetIndex();
		    srcIdx[0] = destIdx[0];
		    srcIdx[1] = destIdx[1];
		    srcIdx[2] = destIdx[2];
	       
		    for (srcIdx[3] = 0; srcIdx[3] < srcSize[3]; srcIdx[3] ++) {
			 timeSeries[srcIdx[3]] = srcPtr->GetPixel( srcIdx );
		    }

		    destIt.Set ( timeSeries );
	       } // mask > 0
	  } // maskIt
	  subIdx ++;
     }
     return numSubs;
}
