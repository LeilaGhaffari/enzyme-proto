// Standard Neo-Hookean formulation in initial configuration with Enzyme-AD

#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ratel-nh-initial.h"

// ----------------------------------------------------------------------------
// Strain Energy Function
// ----------------------------------------------------------------------------
double StrainEnergy(double E_Voigt[6], double lambda, double mu) {
  // Calculate 2*E 
  double E2_Voigt[6];
  for(int i = 0; i<6; i++) E2_Voigt[i] = E_Voigt[i] * 2; 

  // log(J)
  double detCm1 = RatelVoigtDetAM1(E2_Voigt);
  //double J      = sqrt(detCm1 + 1);
  double logJ   = RatelLog1pSeries(detCm1) / 2.;

  // trace(E)
  double traceE = RatelVoigtTrace(E_Voigt);

  //return lambda*(J*J - 1)/4 - lambda*logJ/2  + mu * (-logJ + traceE);
  return lambda*logJ*logJ/2  + mu * (-logJ + traceE);
};

// -----------------------------------------------------------------------------
// Enzyme-AD
// -----------------------------------------------------------------------------
void __enzyme_autodiff(void *, ...);
void __enzyme_augmentfwd(void *, ...);
void __enzyme_fwdsplit(void *, ...);
int  __enzyme_augmentsize(void *, ...);
extern int enzyme_tape, enzyme_const, enzyme_dup, enzyme_nofree, enzyme_allocated;

//  Compute Second Piola Kirchhoff stress:
void SecondPiolaKirchhoffStress_NeoHookean_AD(const double lambda, const double mu, double * __restrict__ E_Voigt, double * __restrict__ S_Voigt) {
  for (int i = 0; i < 6; i++) S_Voigt[i] = 0.;
  __enzyme_autodiff((void *)StrainEnergy, E_Voigt, S_Voigt, enzyme_const, lambda, enzyme_const, mu);
  for (int i = 3; i < 6; i++) S_Voigt[i] /= 2.;
};

void S_fwd(double *S, double *E, const double lambda, const double mu, double *tape) {
  int tape_bytes = __enzyme_augmentsize((void *)SecondPiolaKirchhoffStress_NeoHookean_AD, 
                                        enzyme_const, enzyme_const, enzyme_dup, enzyme_dup);

  __enzyme_augmentfwd((void *)SecondPiolaKirchhoffStress_NeoHookean_AD, 
                      enzyme_allocated, tape_bytes, enzyme_tape, tape, enzyme_nofree, 
                      enzyme_const, lambda, 
                      enzyme_const, mu, 
                      E, (double *)NULL, 
                      S, (double *)NULL);
};

void grad_S(double *dS, double *dE, const double lambda, const double mu, const double *tape) {
  int tape_bytes = __enzyme_augmentsize((void *)SecondPiolaKirchhoffStress_NeoHookean_AD, 
                                         enzyme_const, enzyme_const, enzyme_dup, enzyme_dup);
  for (int i=0; i<6; i++) dS[i] = 0;             
  __enzyme_fwdsplit((void *)SecondPiolaKirchhoffStress_NeoHookean_AD, 
                    enzyme_allocated, tape_bytes, enzyme_tape, tape, 
                    enzyme_const, lambda, 
                    enzyme_const, mu, 
                    (double *)NULL, dE, 
                    (double *)NULL, dS);
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

  // Strain energy
  double strain_energy = StrainEnergy(E_Voigt, mu, lambda);

  // Second Piola-Kirchhoff: S
  double S_ad[6], S_an[6];
  SecondPiolaKirchhoffStress_NeoHookean_AD(lambda, mu, E_Voigt, S_ad); // AD
  SecondPiolaKirchhoffStress_NeoHookean_Analytical(lambda, mu, E_Voigt, S_an); // Analytical

  double S_Voigt[6], tape[6];
  S_fwd(S_Voigt, E_Voigt, lambda, mu, tape);                           

  // Compute delta_S_Voigt with Enzyme-AD
  double delta_S_Voigt[6];
  grad_S(delta_S_Voigt, delta_E_Voigt, lambda, mu, tape);

  printf("\n\nStrain Energy = ");
  printf("\t   %.6lf", strain_energy);

  printf("\n\nS_ad from psi =\n\n");
  for (int i=0; i<6; i++) printf("\t\t%.12lf", S_ad[i]);
  printf("\n\n");

  printf("\n\nS_ad from fwd =\n\n");
  for (int i=0; i<6; i++) printf("\t\t%.12lf", S_Voigt[i]);
  printf("\n\n");

  printf("\n\nS_analytical  =\n\n");
  for (int i=0; i<6; i++) printf("\t\t%.12lf", S_an[i]);
  printf("\n\n");

  printf("\n\ndS_ad         =\n\n");
  for (int i=0; i<6; i++) printf("\t\t%.12lf", delta_S_Voigt[i]);
  printf("\n\n");

  return 0;
}

/* Output:

Strain Energy =            0.922657

S_ad from psi =

                -2.243659385920         -2.164756543395         -0.329653364318         -0.698950459026         1.116803018811          2.683783945834



S_ad from fwd =

                -2.243659385920         -2.164756543395         -0.329653364318         -0.698950459026         1.116803018811          2.683783945834



S_analytical  =

                -2.243659409307         -2.164756566213         -0.329653373905         -0.698950464066         1.116803026863          2.683783965185



dS_ad         =

                2.533390355923          2.921532744535          2.575081303544          1.872573286101          -1.796080925733         -2.446520400374

*/
