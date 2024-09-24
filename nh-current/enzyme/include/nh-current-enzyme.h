#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/nh-current.h"

double  StrainEnergy_NeoHookeanCurrentAD_Enzyme(double e_sym[6], double lambda, double mu) {
  double e2_sym[6];

  // J and log(J)
  for (int i = 0; i < 6; i++) e2_sym[i] = 2 * e_sym[i];
  const double detbm1 = MatDetAM1Symmetric(e2_sym);
  const double J      = sqrt(detbm1 + 1);
  const double logJ   = Log1pSeries(detbm1) / 2.;

  // trace(e)
  const double trace_e = MatTraceSymmetric(e_sym);

  return lambda * (J * J - 1) / 4 - lambda * logJ / 2 + mu * (-logJ + trace_e);
}

/// @}

// -----------------------------------------------------------------------------
//  Compute Kirchhoff stress
// -----------------------------------------------------------------------------
// -- Enzyme-AD
void       __enzyme_autodiff(void *, ...);
void       __enzyme_fwddiff(void *, ...);
extern int enzyme_const;

void *__enzyme_function_like[2] = {(void *)Log1pSeries, (void *)"log1p"};

/// @addtogroup Materials
/// @{

/**
  @brief Compute Kirchoff tau for neo-Hookean hyperelasticity in current configuration with Enzyme

  @param[in]   lambda    Lamé parameter
  @param[in]   mu        Shear modulus
  @param[in]   e_sym     Green Euler strain, in symmetric representation
  @param[out]  tau_sym   Kirchoff tau, in symmetric representation
**/
 void Kirchhofftau_sym_NeoHookean_AD(const double lambda, const double mu, double e_sym[6],
                                                               double tau_sym[6]) {
  double dPsi_sym[6] = {0.}, b_sym[6], dPsi[3][3], b[3][3], tau[3][3];

  // dPsi / de
  __enzyme_autodiff((void *)StrainEnergy_NeoHookeanCurrentAD_Enzyme, e_sym, dPsi_sym, enzyme_const, lambda, enzyme_const, mu);
  for (int i = 3; i < 6; i++) dPsi_sym[i] /= 2.;

  // b = 2 e + I
  for (int j = 0; j < 6; j++) b_sym[j] = 2 * e_sym[j] + (j < 3);

  // tau = (dPsi / de) b
  SymmetricMatUnpack(dPsi_sym, dPsi);
  SymmetricMatUnpack(b_sym, b);
  MatMatMult(1., dPsi, b, tau);
  SymmetricMatPack(tau, tau_sym);
}

void dtau_fwd(const double lambda, const double mu, double e_sym[6], double de_sym[6],
                                         double tau_sym[6], double dtau_sym[6]) {
  __enzyme_fwddiff((void *)Kirchhofftau_sym_NeoHookean_AD, enzyme_const, lambda, enzyme_const, mu, e_sym, de_sym, tau_sym, dtau_sym);
}
