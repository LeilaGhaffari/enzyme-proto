// Tested function: 2.*x*x*x
// First derivative: 2.*3.*x*x
// Second derivative: 2.*3.*2.*x

#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>
#include <adolc/adolc.h>

using namespace std;

int main () {

  double x1 = 3.;
  adouble ax1;
  double y1;
  adouble ay1;

  trace_on(1); // not sure what the input should be!
  ax1 <<= x1;

  ay1 = 2.*ax1*ax1*ax1;

  ay1 >>= y1;
  trace_off();

  double* yprim;
  yprim = myalloc1(1);
  yprim[0] = 2.*x1*x1*x1;


  double*** yDerivative;
  yDerivative = myalloc3(1, 3, 2);
  yDerivative[0][0][0] = 2.*3.*x1*x1;
  yDerivative[0][0][1] = 2.*3.*x1*x1 + 0.5*(2.*3.*2.*x1); // And what is this?
  yDerivative[0][1][0] = 2.*2.*3.*x1*x1;
  yDerivative[0][1][1] = 2.*2.*3.*x1*x1 + 0.5*(2.*3.*2.*x1)*2.*2.;
  yDerivative[0][2][0] = 3.*2.*3.*x1*x1;
  yDerivative[0][2][1] = 3.*2.*3.*x1*x1 + 0.5*(2.*3.*2.*x1)*3.*3.;

  double* x;
  x = myalloc1(1);
  x[0] = 3.;

  double*** X;
  X = myalloc3(1, 3, 2);
  X[0][0][0] = 1.;
  X[0][1][0] = 2.;
  X[0][2][0] = 3.;
  X[0][0][1] = 1.;
  X[0][1][1] = 2.;
  X[0][2][1] = 3.;

  double* y;
  y = myalloc1(1);

  double*** Y;
  Y = myalloc3(1, 3, 2);

  hov_forward(1, 1, 1, 2, 3, x, X, y, Y);

  for (int i=0; i<3; i++) for (int j=0; j<2; j++) cout << "Y[0, " << i << ", " << j << "] = " << Y[0][i][j] << "\t" << yDerivative[0][i][j] << endl;

  myfree1(yprim);
  myfree3(yDerivative);
  myfree1(x);
  myfree1(y);
  myfree3(X);
  myfree3(Y);
}
