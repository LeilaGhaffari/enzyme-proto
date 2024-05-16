// Standard Neo-Hookean formulation in initial configuration with ADOL-C

#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>
#include <adolc/adolc.h>
#include "ratel-nh-initial.hpp"
#include "ratel-nh-initial-adolc.hpp"

using namespace std;

int main() {
  int m = 1, n = 6;
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
  double    delta_E_Voigt[n];
  for (int mm = 0; mm < n; mm++) {
    delta_E_Voigt[mm] = 0;
    for (int nn = 0; nn < 3; nn++) {
      delta_E_Voigt[mm] += (grad_delta_u[nn][ind_j[mm]] * F[nn][ind_k[mm]] + F[nn][ind_j[mm]] * grad_delta_u[nn][ind_k[mm]]) / 2.;
    }
  }

  // J-1
  const double Jm1 = RatelMatDetAM1(temp_grad_u);

  // Green Lagrange Strain Tensor: E (voigt)
  double E_Voigt[n];
  GreenLagrangeStrain(temp_grad_u, E_Voigt);

  cout.precision(12);
  cout.setf(ios::fixed);

  // Strain energy
  cout << "\nStrain Energy = " << StrainEnergy(E_Voigt, lambda, mu) << endl;
  cout << endl;

  // Second Piola-Kirchhoff: S
  double S_an[n];
  SecondPiolaKirchhoffStress_NeoHookean_Analytical(E_Voigt, S_an, lambda, mu); // Analytical
  cout << "S: analytical  =" <<endl;
  for (int i=0; i<n; i++) cout << S_an[i] << endl;
  cout << endl;

  // ------------------------------------------------------------------------
  // Active section for automatic differentiation
  // ------------------------------------------------------------------------
  // Initialize passive variables
  auto Ep = new double[n];    // Independent vector
  auto Fp = new double[m];    // Dependent vector
  for (int i=0; i<n; i++) Ep[i] = E_Voigt[i];

  // Initialize active variables
  auto Ea = new adouble[n];
  auto Fa = new adouble[m];

  // Set the tag for the Automatic Differentiation trace
  int tag = 0;

  // Start tracing floating point operations
  trace_on(tag);  // Start of the active section

  // Assign independent variables
  for (int i=0; i<n; i++) Ea[i] <<= Ep[i];

  // Evaluate the body of the differentiated code
  Fa[0] = StrainEnergy(Ea, lambda, mu);

  // Assign dependent variables
  Fa[0] >>= Fp[0];
  trace_off();    // End of the active section

  // ------------------------------------------------------------------------
  // Compute the first derivative (forward scalar mode)
  // ------------------------------------------------------------------------
  // Declare the tangent vector
  auto E_fwd = new double[n];
  for (int i=0; i<n; i++) E_fwd[i] = 0.;

  // Declare the vector of first derivatives
  auto S_fwd = new double[m];

  // Define a flag to prepare for a reverse automatic differentiation
  int keep = 1;

  // Compute the derivative
  cout << "S: forward mode AD  =" << endl;
  for (int i=0; i<n; i++) {
    E_fwd[i] = 1.;
    fos_forward(tag, m, n, keep, Ep, E_fwd, Fp, S_fwd);
    if (i > 2) S_fwd[0] /= 2.;
    cout << S_fwd[0] << endl;
    E_fwd[i] = 0.;
  }
  cout << endl;

  // ------------------------------------------------------------------------
  // Compute the first derivative of the Strain Energy (reverse scalar mode)
  // ------------------------------------------------------------------------
  // Declare the weight vector
  auto f = new double[m];
  // Declare the vector of first derivatives (adjoint vector)
  auto S_rev = new double[n];
  // Compute the derivatives of f(x,y,z)
  f[0] = 1;
  fos_reverse(tag, m, n, f, S_rev);
  for (int i = 3; i < n; i++) S_rev[i] /= 2.;

  cout << "S: reverse mode AD  =" << endl;
  for (int i=0; i<n; i++) cout << S_rev[i] << endl;
  cout << endl;

  return 0;
}

/*
Strain Energy = 0.922656758772

S =
-2.243659409307
-2.164756566213
-0.329653373905
-0.698950464066
1.116803026863
2.683783965185

dS =
2.533391359245
2.921533726694
2.575081725908
1.872573511865
-1.796081277841
-2.446521233043
*/
