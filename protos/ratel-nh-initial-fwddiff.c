// Standard Neo-Hookean formulation in initial configuration with Enzyme-AD

#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ratel-nh-initial.h"

// -----------------------------------------------------------------------------
// Enzyme-AD
// -----------------------------------------------------------------------------
void __enzyme_fwddiff(void *, ...);
extern int enzyme_tape, enzyme_const, enzyme_dup, enzyme_nofree, enzyme_allocated;

void *__enzyme_function_like[2] = {(void *)RatelLog1pSeries, "log1p"};

void grad_S(double *dS, double *dE, double *S, double *E, const double lambda, const double mu) {
  __enzyme_fwddiff((void *)SecondPiolaKirchhoffStress_NeoHookean_Analytical,
                    enzyme_const, lambda,
                    enzyme_const, mu,
                    E, dE,
                    S, dS);
};

// -----------------------------------------------------------------------------
//  Main
// -----------------------------------------------------------------------------
int main() {
  double mu = 1.;
  double lambda = 1.;

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
  const int ind_j[6] = {0, 1, 2, 1, 0, 0}, ind_k[6] = {0, 1, 2, 2, 2, 1};
  double    delta_E_Voigt[6];
  for (int m = 0; m < 6; m++) {
    delta_E_Voigt[m] = 0;
    for (int n = 0; n < 3; n++) {
      delta_E_Voigt[m] += (grad_delta_u[n][ind_j[m]] * F[n][ind_k[m]] + F[n][ind_j[m]] * grad_delta_u[n][ind_k[m]]) / 2.;
    }
  }

  // J-1
  const double Jm1 = RatelMatDetAM1(temp_grad_u);

  // Green Lagrange Strain Tensor: E (voigt)
  double E_Voigt[6];
  GreenLagrangeStrain(temp_grad_u, E_Voigt);

  // Compute delta_S_Voigt with Enzyme-AD
  double S_Voigt[6], delta_S_Voigt[6];
  grad_S(delta_S_Voigt, delta_E_Voigt, S_Voigt, E_Voigt, lambda, mu);

  printf("\n\nS_ad         =\n\n");
  for (int i=0; i<6; i++) printf("\t\t%.12lf", S_Voigt[i]);
  printf("\n\n");

  printf("\n\ndS_ad         =\n\n");
  for (int i=0; i<6; i++) printf("\t\t%.12lf", delta_S_Voigt[i]);
  printf("\n\n");

  return 0;
}

/* Output:

S_ad         =

                -2.243659409307         -2.164756566213         -0.329653373905         -0.698950464066         1.116803026863          2.683783965185

dS_ad         =

                2.533391374087          2.921533741175          2.575081731992          1.872573515063          -1.796081282951         -2.446521245324
*/
