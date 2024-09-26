// Standard Neo-Hookean formulation in initial configuration with ADOL-C

#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>
#include <adolc/adolc.h>

using namespace std;
static int m = 1, //number of dependent variables
           n = 6; //number of independent variables

#include "../include/ratel-nh-initial.hpp"
#include "../include/ratel-nh-initial-adolc.hpp"



int main() {
  const double mu = 1., lambda = 1.;

  double grad_u[3][3] = {
    {0.0702417, 0.4799115, 0.3991242},
    {0.6756593, 0.0633284, 0.0959267},
    {0.2241923, 0.0281781, 0.0917613}
  };

  double grad_delta_u[3][3] = {
    {0.1425560,  0.115120,  0.551640},
    {0.0591922,  0.123535,  0.166572},
    {0.1617210,  0.478828,  0.646217}
  };

  // Compute the Deformation Gradient : F = I + grad_u
  const double F[3][3] = {
    {grad_u[0][0] + 1, grad_u[0][1],     grad_u[0][2]    },
    {grad_u[1][0],     grad_u[1][1] + 1, grad_u[1][2]    },
    {grad_u[2][0],     grad_u[2][1],     grad_u[2][2] + 1}
  };

  const double temp_grad_u[3][3] = {
    {grad_u[0][0], grad_u[0][1], grad_u[0][2]},
    {grad_u[1][0], grad_u[1][1], grad_u[1][2]},
    {grad_u[2][0], grad_u[2][1], grad_u[2][2]}
  };

  // deltaE - Green-Lagrange strain tensor
  const int ind_j[n] = {0, 1, 2, 1, 0, 0}, ind_k[n] = {0, 1, 2, 2, 2, 1};
  double    dE_sym[n];
  for (int mm = 0; mm < n; mm++) {
    dE_sym[mm] = 0;
    for (int nn = 0; nn < 3; nn++) {
      dE_sym[mm] += (grad_delta_u[nn][ind_j[mm]] * F[nn][ind_k[mm]] + F[nn][ind_j[mm]] * grad_delta_u[nn][ind_k[mm]]) / 2.;
    }
  }

  // Green Lagrange Strain Tensor: E
  double E_sym[n];
  GreenLagrangeStrain(temp_grad_u, E_sym);

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
