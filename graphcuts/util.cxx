#include <common.h>

int label2hidden(MatrixXf & image)
{
     for (int i = 0; i < image.rows(); i++) {
	  for (int j = 0; j < image.cols(); j ++) {
	       image(i,j) = (image(i,j) == 0)?(-1):1;
	  }
     }
     return 0;
}

int hidden2label(MatrixXf & image)
{
     for (int i = 0; i < image.rows(); i++) {
	  for (int j = 0; j < image.cols(); j ++) {
	       image(i,j) = (image(i,j) == -1)?0:1;
	  }
     }
     return 0;
}

