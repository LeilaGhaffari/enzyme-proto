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
  double J      = sqrt(detCm1 + 1);
  double logJ   = RatelLog1pSeries(detCm1) / 2.;

  // trace(E)
  double traceE = RatelVoigtTrace(E_Voigt);

  //return lambda*(J*J - 1)/4 - lambda*logJ/2  + mu * (-logJ + traceE);
  return lambda*logJ*logJ/2  + mu * (-logJ + traceE);
};

// -----------------------------------------------------------------------------
//  Compute Second Piola Kirchhoff stress:
// -----------------------------------------------------------------------------
void __enzyme_autodiff(void *, ...);
extern int enzyme_const;

void SecondPiolaKirchhoffStress_NeoHookean_AD(const double lambda, const double mu, double E_Voigt[6], double S_Voigt[6]) {
  for (int i=0; i<FSInitialNH_AD_TAPE_SIZE; i++) S_Voigt[i] = 0.;
  __enzyme_autodiff((void *)StrainEnergy, 
                     E_Voigt, S_Voigt,
                     enzyme_const, mu,
                     enzyme_const, lambda);
  for (int i=FSInitialNH_AD_TAPE_SIZE/2; i<FSInitialNH_AD_TAPE_SIZE; i++) S_Voigt[i] /= 2.;
}

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
   
  // J-1
  const double Jm1 = RatelMatDetAM1(temp_grad_u);

  // Green Lagrange Strain Tensor: E (voigt)
  double E_Voigt[6];
  GreenLagrangeStrain(temp_grad_u, E_Voigt);

  // Strain energy
  double strain_energy = StrainEnergy(E_Voigt, mu, lambda);

  // Second Piola-Kirchhoff: S
  double S_Voigt[6], S_an[6], S[3][3];
  SecondPiolaKirchhoffStress_NeoHookean_AD(lambda, mu, E_Voigt, S_Voigt); // AD
  S_analytical(S_an, E_Voigt, mu, lambda); // Analytical
  RatelVoigtUnpack(S_Voigt, S); // Unpack Voigt S

  // First Piola-Kirchhoff: P = F*S
  double P[3][3];
  RatelMatMatMult(1.0, F, S, P);

  // Print Results
  printf("\nStrain Energy = ");
  printf("%.6lf", strain_energy);
  printf("\nE_Voigt       =\n");
  for (int i=0; i<FSInitialNH_AD_TAPE_SIZE; i++) printf("\t\t%.12lf", E_Voigt[i]);
  printf("\n");
  printf("\nS_autodiff    =\n");
  for (int i=0; i<FSInitialNH_AD_TAPE_SIZE; i++) printf("\t\t%.12lf", S_Voigt[i]);
  printf("\n");
  printf("\nS_analytical  =\n");
  for (int i=0; i<FSInitialNH_AD_TAPE_SIZE; i++) printf("\t\t%.12lf", S_an[i]);
  printf("\n");
  return 0;
}

/* Output:

Strain Energy = 0.922657
E_Voigt       =
                0.326097486737          0.180888169699          0.180222397488          0.162154818512          0.368368803095          0.619193167536

S_autodiff    =
                -2.243659385920         -2.164756543395         -0.329653364318         -0.698950459026         1.116803018811          2.683783945834

S_analytical  =
                -2.243659409307         -2.164756566213         -0.329653373905         -0.698950464066         1.116803026863          2.683783965185
*/
