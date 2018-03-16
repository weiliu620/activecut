
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "all.h"

int main(int argc, char *argv[])
{

	const char *  inFile = "/home/sci/weiliu/fmri/restdata_nii/s029a001.nii";
	const char *  outFile = "toy.nii";

	typedef  short PixelType;
	const unsigned int Dimension = 4;
	typedef itk::Image< PixelType, Dimension > ImageType;
	typedef itk::ImageFileWriter< ImageType > WriterType;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	WriterType::Pointer writer = WriterType::New();

	reader->SetFileName( inFile);
	writer->SetFileName(outFile);

	reader->Update();
	ImageType::Pointer image = reader->GetOutput();
// 	ImageType::Pointer image = ImageType::New();
	writer->SetInput(image);	//Setinput need a smart pointer as argument. see itk mannual p136.
	//itk::NiftiImageIO::Pointer niftiIO = itk::NiftiImageIO::New(); 
	//writer->SetImageIO(niftiIO);

// 	ImageType::IndexType start;
// 	start[0] = 0; // first index on X
// 	start[1] = 0; // first index on Y
// 	start[2] = 0; // first index on Z
// 	start[3] = 0;
// 	ImageType::SizeType size;
// 	size[0] = colNum; // size along X
// 	size[1] = rowNum; // size along Y
// 	size[2] = sliceNum; // size along Z
// 	size[3] = volNum;
// 
// 	ImageType::RegionType region;
// 	region.SetSize( size );
// 	region.SetIndex( start );
// 	image->SetRegions( region );
// 	image->Allocate();
	ImageType::IndexType pixelIndex;
	ImageType::SizeType imageSize = image->GetLargestPossibleRegion().GetSize();

	const int sliceNum = imageSize[2];
	const int rowNum = imageSize[1];
	const int colNum = imageSize[0];
	const int volNum = imageSize[3];

	
	int i, j, k, t;
	short pixelValue;
	const short toyMean = 100;
	double * ranNumber = NULL;
	
	/** Generate random number in normal distribution*/
	short * ranArray = NULL;
	long int ranSize = volNum * sliceNum * ceil(colNum/2) * ceil(rowNum/2);
	ranArray = new short[ranSize];
	GetNormalSample(ranArray, ranSize);
	short * currRanNum = ranArray;
	for (t = 0; t < volNum; t++){
		pixelIndex[3] = t;
		for (k = 0; k < sliceNum; k++){
			pixelIndex[2] = k;
			for (i = 0; i < rowNum/2; i++){
				pixelIndex[1] = i;
				for (j = 0; j < colNum/2; j++){
					pixelIndex[0] = j;
					pixelValue = toyMean + sin(M_PI*t/5)*(toyMean/5);
					image->SetPixel( pixelIndex,pixelValue);
				}
				for (j = colNum/2; j < colNum; j++){
					pixelIndex[0] = j;
					pixelValue = toyMean - sin(M_PI*t/5)*(toyMean/5);
					image->SetPixel( pixelIndex,pixelValue);
				}
			}
			for (i = rowNum/2; i < rowNum; i++){
				pixelIndex[1] = i;
				for (j = 0; j < colNum/2; j++){
					pixelIndex[0] = j;
					//pixelValue = toyMean + sin(M_PI*t/5)*(toyMean/5) + GetNormalSample(10);
					pixelValue = toyMean + sin(M_PI*t/5)*(toyMean/5) + *(currRanNum++) * 4;
					image->SetPixel( pixelIndex,pixelValue);
				}
				for (j = colNum/2; j < colNum; j++){
					pixelIndex[0] = j;
					pixelValue = toyMean;
					image->SetPixel( pixelIndex,pixelValue);
				}
			}
		}
	if (t%10 == 0) 	std::cout << "generating volume " << t << " of " << volNum << std::endl;
	}
	delete ranArray;




	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}


}
