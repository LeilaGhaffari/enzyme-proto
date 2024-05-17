// Tested function: 2.*x*x*x
// First derivative: 2.*3.*x*x
// Second derivative: 2.*3.*2.*x

#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>
#include <adolc/adolc.h>

using namespace std;

adouble my_func(adouble x[3]) {
  adouble u = x[0], v= x[1], w = x[2];
  return 5*u*u*v*w + 3*u*v*v*w + 6*w*w*v + 1.;
}

double my_func(double x[3]) {
  double u = x[0], v= x[1], w = x[2];
  return 5*u*u*u*u + 3*v*v*v + 6*w*w + 1.;
}

int main () {

  int m = 1, n = 3, d = 2, p =3, tag = 1; // not sure about tag

  // Initialize passive variables
  auto xp = new double[n];    // Independent vector
  auto yp = new double[m];    // Dependent vector
  xp[0] = 1., xp[1] = 2., xp[2] = 3.;

  // Initialize active variables
  auto xa = new adouble[n];
  auto ya = new adouble[m];

  trace_on(tag);
  for (int i=0; i<n; i++) xa[i] <<= xp[i];

  ya[0] = my_func(xa);
  ya[0] >>= yp[0];

  trace_off();

  double ***X = myalloc3(n, p, d);
  for (int i=0; i<n; i++)
    for (int j=0; j<p; j++)
      if (i == j)
        X[i][j][0] = 1.;
      else
        X[i][j][0] = 0.;  // Last index [0], activate only the x1 direction

  double ***Y;
  Y = myalloc3(m, p, d);

  hov_forward(tag, m, n, d, p, xp, X, yp, Y);

  for (int i=0; i<m; i++)
    for (int j=0; j<p; j++)
      for (int k=0; k<d; k++)
        cout << "Y[" << i << ", " << j << ", " << k << "] = " << Y[i][j][k] << endl;

  //myfree1(yprim);
  //myfree3(yDerivative);
  //myfree1(x);
  //myfree1(y);
  myfree3(X);
  myfree3(Y);
}
