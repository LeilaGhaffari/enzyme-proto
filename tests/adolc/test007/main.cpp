// Refactored vector example (vector-input-scalar-output)

#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>
#include <adolc/adolc.h>

using namespace std;

static int n = 3;
static int m = 1;

adouble func(adouble x[3]) {
  adouble u = x[0], v= x[1], w = x[2];
  return 5*u*u*v*w + 3*u*v*v*w + 6*w*w*v + 1.;
}

double func(double x[3]) {
  double u = x[0], v= x[1], w = x[2];
  return 5*u*u*v*w + 3*u*v*v*w + 6*w*w*v + 1.;
}

double *dfunc(double xp[3]) {
  int tag = 1, keep = 1;
  auto xa = new adouble[n];
  auto ya = new adouble[m];
  auto yp = new double[m];

  trace_on(tag);
  for (int i=0; i<n; i++) xa[i] <<= xp[i];
  ya[0] = func(xa);
  ya[0] >>= yp[0];
  trace_off();

  auto Y = new double[m];
  auto X = new double[n];
  auto dX = new double[n];
  for (int i=0; i<n; i++) X[i] = 0.;

  for (int i=0; i<n; i++) {
    X[i] = 1.;
    fos_forward(tag, m, n, keep, xp, X, yp, Y);
    X[i] = 0.;
    dX[i] = Y[0];
  }
  return dX;
}

int main () {

  // Define and initialize independent variable
  auto xp = new double[n];
  xp[0] = 1., xp[1] = 2., xp[2] = 3.;

  // Function evaluation
  auto f = func(xp);

  // First derivative
  auto dX = dfunc(xp);

  // Print
  cout << "f = " << f << endl << endl;
  cout << "dX =" << endl;
  for (int i=0; i<n; i++) cout << dX[i] << endl;
  cout << endl;

  return 0;
}

/*
f = 175

dX =
96
105
94
*/
