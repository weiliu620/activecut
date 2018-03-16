#include <Eigen/Core>

#include <stdio.h>
#include <new>
#include <iostream>
#include <string>
#include <fstream>
#include <Eigen/Core>

#include <boost/math/distributions/normal.hpp> // for normal_distribution.
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>

typedef boost::mt19937 twister_base_gen_type;
using namespace std;
// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

int main(int, char *[])
{
//  Matrix3f m3;
     Matrix<float, 3, 3> m3;
  m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  Matrix4f m4 = Matrix4f::Identity();
  Vector4i v4(1, 2, 3, 4);

  std::cout << "m3\n" << m3 << "\nm4:\n"
    << m4 << "\nv4:\n" << v4 << std::endl;
}
