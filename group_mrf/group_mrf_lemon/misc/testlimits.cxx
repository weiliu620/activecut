// numeric_limits example
#include <iostream>
#include <limits>
#include <cmath>
using namespace std;

int main () {
  cout << boolalpha;
  cout << "Minimum value for double: " << numeric_limits<double>::min() << endl;
  cout << "Maximum value for double: " << numeric_limits<double>::max() << endl;
  cout << "double is signed: " << numeric_limits<double>::is_signed << endl;
  cout << "Non-sign bits in double: " << numeric_limits<double>::digits << endl;
  cout << "double has infinity: " << numeric_limits<double>::has_infinity << endl;

  double a = exp(-3090);
  cout << a << endl;
  return 0;
}
