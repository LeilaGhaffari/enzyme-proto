// Standard Neo-Hookean formulation in initial configuration with ADOL-C

#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>
#include <adolc/adolc.h>

using namespace std;
static int m = 1, // number of dependent variables
           n = 6; // number of independent variables

#include "../include/nh-common.h"
#include "../include/nh-adolc.h"

int main() {
  const double mu = 1., lambda = 1.;

  // Green Lagrange Strain Tensor: E
  double E_sym[6];
  E_sym[0] = 0.356722230881;
  E_sym[1] = 0.053880729108;
  E_sym[2] = 0.192505475424;
  E_sym[3] = 0.105608348553;
  E_sym[4] = 0.355594178128;
  E_sym[5] = 0.150573971552;

  double dE_sym[6];
  dE_sym[0] = 0.705072432429;
  dE_sym[1] = 0.479848109568;
  dE_sym[2] = 0.127687977309;
  dE_sym[3] = 0.263932661797;
  dE_sym[4] = 0.307099644410;
  dE_sym[5] = 0.338054016477;

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
  cout << "\n  S (ADOL-C gradient):" << endl;
  for (int i=0; i<n; i++) cout << "\t" << S_sym[i] << endl;
  cout << "\n dS (ADOL-C hessian) =" << endl;
  for (int i=0; i<n; i++) cout << "\t" << dS[i] << endl;
  cout << endl;

  return 0;
}

/*
 Strain energy = 0.338798323728

  S (ADOL-C gradient):
        0.603043653441
        0.502395201627
        0.518371210470
        0.039367933591
        0.192987413544
        0.071116775739

 dS (ADOL-C hessian) =
        1.349719621895
        1.751709039535
        1.340798936767
        -0.027601318916
        -0.553382183654
        -0.172707056259
*/
