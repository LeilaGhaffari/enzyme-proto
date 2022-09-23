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

  double E_Voigt[6];
  E_Voigt[0] = 0.5895232828911128/2;
  E_Voigt[1] = 0.2362491738162759/2;
  E_Voigt[2] = 0.9793730522395296/2;
  E_Voigt[3] = 0.2190993957421843/2;
  E_Voigt[4] = 0.0126503210747925/2;
  E_Voigt[5] = 0.6570956167695403/2;

  // Compute S with Enzyme-AD Forward mode
  double strain_energy = StrainEnergy(E_Voigt, mu, lambda);
  double S_Voigt[6];
  SecondPiolaKirchhoffStress_NeoHookean_AD(lambda, mu, E_Voigt, S_Voigt);

  // Compute analytical S
  double S_an[6];
  S_analytical(S_an, E_Voigt, mu, lambda);

  printf("\n\nStrain Energy   = ");
  printf("\t   %.6lf", strain_energy);

  printf("\n\nS_autodiff      =\n\n");
  for (int i=0; i<FSInitialNH_AD_TAPE_SIZE; i++) printf("\t\t%.12lf", S_Voigt[i]);
  printf("\n\n");

  printf("\n\nS_analytical    =\n\n");
  for (int i=0; i<FSInitialNH_AD_TAPE_SIZE; i++) printf("\t\t%.12lf", S_an[i]);
  printf("\n\n");

  return 0;
}

/* Output:

Strain Energy   =          0.507029

S_autodiff      =

                0.629897288210          0.514638184498          0.763455742111          0.052445668610          -0.019798047937         0.200227118654



S_analytical    =

                0.629816767260          0.514532587340          0.763404278645          0.052457078889          -0.019802355275         0.200270680826
*/
