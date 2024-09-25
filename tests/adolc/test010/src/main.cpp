// Standard Neo-Hookean formulation in initial configuration with ADOL-C

#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>
#include <adolc/adolc.h>

using namespace std;
static int m = 1, n = 6, d = 2, p = 6;

#include "../include/ratel-nh-initial.hpp"
#include "../include/ratel-nh-initial-adolc.hpp"



int main() {
  const double mu = 1., lambda = 1.;
  int size = binomi(p + d, d), temp = 1;

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
  double    delta_E_Voigt[n];
  for (int mm = 0; mm < n; mm++) {
    delta_E_Voigt[mm] = 0;
    for (int nn = 0; nn < 3; nn++) {
      delta_E_Voigt[mm] += (grad_delta_u[nn][ind_j[mm]] * F[nn][ind_k[mm]] + F[nn][ind_j[mm]] * grad_delta_u[nn][ind_k[mm]]) / 2.;
    }
  }

  // Green Lagrange Strain Tensor: E (voigt)
  double E_Voigt[n];
  GreenLagrangeStrain(temp_grad_u, E_Voigt);

  // Second Piola-Kirchhoff: S
  double S_an[n];
  SecondPiolaKirchhoffStress_NeoHookean_Analytical(E_Voigt, S_an, lambda, mu); // Analytical

  // ------------------------------------------------------------------------
  // Automatic Differentiation
  // ------------------------------------------------------------------------
  // Initialize passive variables
  auto Ep = new double[n];
  for (int i=0; i<n; i++) Ep[i] = E_Voigt[i];

  // First derivative (forward vector mode)
  auto S_fwd = Stress(Ep, lambda, mu);

  // High order tensor
  // -----------------
  double F_tensor;
  double S_tensor[n] = {0.};
  double dS_tensor[n][n] ={{0.}};
  double dS[n] = {0.};

  // Tensor-AD evaluation
  auto tensor = dStress(Ep, lambda, mu);

  // Strain energy
  F_tensor = tensor[0];

  for (int i=0; i<n; i++) {
    S_tensor[i] = tensor[temp]; // Populate stress (1st derivative)
    if (i>2) S_tensor[i] /= 2.;
    for (int j=0; j<i+1; j++) {
      dS_tensor[i][j] = tensor[temp+j+1]; // Populate delta_S (2nd derivative)
      if (i != j) dS_tensor[j][i] = dS_tensor[i][j];
    }
    temp += i + 2;
  }
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) if (i > 2) dS_tensor[i][j] /= 2.;

  // Compute dS
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) dS[i] += dS_tensor[i][j] * delta_E_Voigt[j];
  }

  // ------------------------------------------------------------------------
  // Print
  // ------------------------------------------------------------------------
  cout.precision(12);
  cout.setf(ios::fixed);
  cout << "\n size = " << size << endl;
  cout << "\n Strain energy = " << F_tensor << endl;
  cout << "\n Stress =" << endl;
  cout << "\n   Forward Vector Mode:" << endl << endl;
  for (int i=0; i<n; i++) cout << "\t" << S_fwd[i] << endl;
  cout << endl;
  cout << "\n   Higher Order Tensor:" << endl << endl;
  for (int i=0; i<n; i++) cout << "\t" << S_tensor[i] << endl;
  cout << endl;
  cout << "\n dS =" << endl << endl;
  for (int i=0; i<n; i++) cout << "\t" << dS[i] << endl;
  cout << endl;

  return 0;
}

/*
 size = 28

 Strain energy = 0.922656758772

 Stress =

   Forward Vector Mode:

        -2.243659385920
        -2.164756543395
        -0.329653364318
        -0.698950459026
        1.116803018811
        2.683783945834


   Higher Order Tensor:

        -2.243659385920
        -2.164756543395
        -0.329653364318
        -0.698950459026
        1.116803018811
        2.683783945834


 dS =

        2.533390355923
        2.921532744535
        2.575081303544
        1.872573286101
        -1.796080925733
        -2.446520400374
*/
