// Standard Neo-Hookean formulation in initial configuration with ADOL-C

#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>
#include <adolc/adolc.h>

using namespace std;
static int m = 1, // number of dependent variables
           n = 6; // number of independent variables

#include "nh-common.h"
#include "nh-adolc.h"



int main() {
  const double mu = 1., lambda = 1.;

  // Green Lagrange Strain Tensor: E
  double E_sym[6];
  E_sym[0] = 0.119405838995;
  E_sym[1] = 0.022906026341;
  E_sym[2] = 0.038223518388;
  E_sym[3] = 0.037745222581;
  E_sym[4] = 0.185032644972;
  E_sym[5] = 0.068185405852;

  double dE_sym[6];
  dE_sym[0] = 0.360249242839;
  dE_sym[1] = 0.414569769236;
  dE_sym[2] = -0.02167824525;
  dE_sym[3] = 0.158595277365;
  dE_sym[4] = 0.060981222528;
  dE_sym[5] = 0.186598372502;

  // ------------------------------------------------------------------------
  // Automatic Differentiation
  // ------------------------------------------------------------------------
  // Initialize passive variables
  auto Ep = new double[n];
  for (int i=0; i<n; i++) Ep[i] = E_sym[i];

  // Compute the gradient of Psi (S = dPsi/dE)
  double S_sym[6] = {0.};
  ComputeGradPsi(S_sym, Ep, lambda, mu);

  // Compute the hessian of Psi (d2Psi/dE2)
  double hessPsi[6][6] = {{0.}};
  ComputeHessianPsi(hessPsi, Ep, lambda, mu);

  // Compute dS
  double dS[6] = {0.};
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) dS[i] += hessPsi[i][j] * dE_sym[j];

  // ------------------------------------------------------------------------
  // Print
  // ------------------------------------------------------------------------
  cout.precision(12);
  cout.setf(ios::fixed);
  cout << "\n Strain energy = " << StrainEnergy(E_sym, lambda, mu) << endl;
  cout << "\n  S (ADOL-C gradient):" << endl << endl;
  for (int i=0; i<n; i++) cout << "\t" << S_sym[i] << endl;
  cout << endl;
  cout << "\n dS (ADOL-C hessian) =" << endl << endl;
  for (int i=0; i<n; i++) cout << "\t" << dS[i] << endl;
  cout << endl;

  return 0;
}

/*
  Strain energy = 0.922656758772

  S (ADOL-C gradient):

        -2.138040990828
        -2.061707342723
        -0.286357864495
        -0.676191587965
        1.080438244201
        2.596395931432


 dS (ADOL-C hessian) =

        1.804379625553
        2.195610276263
        2.226208937974
        1.672286083272
        -1.514998884900
        -1.831931547889
*/
