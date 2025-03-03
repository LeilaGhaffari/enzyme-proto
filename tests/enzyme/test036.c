#include "test036.h"

// Strain Energy
double RatelStrainEnergy_IsochoricOgdenCurrentAD_Enzyme(const double bulk, const int N, const double *m, const double *alpha, double e_sym[6]) {
  double e2_sym[6], e_vals[3], pr_str[3], pr_str_bar[3], V, strain_energy;

  // J
  for (int i = 0; i < 6; i++) e2_sym[i] = 2 * e_sym[i];
  const double detbm1 = RatelMatDetAM1Symmetric(e2_sym);
  const double J      = sqrt(detbm1 + 1);

  RatelMatComputeEigensystemSymmetric(e_sym, e_vals, NULL);
  RatelPrincipalStretch(e_vals, pr_str);
  const double J_pow = pow(J, -1. / 3.);

  RatelScalarVecMult(J_pow, pr_str, pr_str_bar);
  VolumetricFunctionAndDerivatives(J - 1., &V, NULL, NULL);
  RatelStrainEnergy_IsochoricOgden(V, bulk, N, m, alpha, pr_str_bar, &strain_energy);

  return strain_energy;
}

// -- Enzyme-AD
extern int enzyme_const;
void       __enzyme_autodiff(void *, ...);
void       __enzyme_fwddiff(void *, ...);
void      *__enzyme_function_like[2] = {(void *)RatelLog1pSeries, (void *)"log1p"};

// Compute \tau
void RatelKirchhofftau_sym_IsochoricOgden_AD(const double bulk, const int N, const double *m, const double *alpha, double e_sym[6],
                                             double tau_sym[6]) {
  double dPsi_sym[6] = {0.}, b_sym[6], dPsi[3][3], b[3][3], tau[3][3];

  // dPsi / de
  __enzyme_autodiff((void *)RatelStrainEnergy_IsochoricOgdenCurrentAD_Enzyme, enzyme_const, bulk, enzyme_const, N, enzyme_const, m, enzyme_const,
                    alpha, e_sym, dPsi_sym);
  for (int i = 3; i < 6; i++) dPsi_sym[i] /= 2.;

  // b = 2 e + I
  for (int j = 0; j < 6; j++) b_sym[j] = 2 * e_sym[j] + (j < 3);

  // tau = (dPsi / de) b
  RatelSymmetricMatUnpack(dPsi_sym, dPsi);
  RatelSymmetricMatUnpack(b_sym, b);
  RatelMatMatMult(1., dPsi, b, tau);
  RatelSymmetricMatPack(tau, tau_sym);
}

// Compute d \tau
void Rateldtau_fwd(const double bulk, const double N, const double *m, const double *alpha, double e_sym[6], double de_sym[6], double tau_sym[6],
                   double dtau_sym[6]) {
  __enzyme_fwddiff((void *)RatelKirchhofftau_sym_IsochoricOgden_AD, enzyme_const, bulk, enzyme_const, N, enzyme_const, m, enzyme_const, alpha, e_sym,
                   de_sym, tau_sym, dtau_sym);
}

int main() {
  // Constants
  double bulk = 1.0, m[3] = {.1, .2, .3}, alpha[3] = {1., 2., 3.};
  int    N = 3;

  double e_sym[6]  = {0.0702417, 0.4799115, 0.3991242, 0.6756593, 0.0633284, 0.0959267};
  double de_sym[6] = {0.1425560, 0.115120, 0.551640, 0.0591922, 0.123535, 0.166572};

  // Strain energy
  double strain_energy = RatelStrainEnergy_IsochoricOgdenCurrentAD_Enzyme(bulk, N, m, alpha, e_sym);
  printf("\n\nStrain Energy = %f \n\n", strain_energy);

  // tau
  double tau_sym[6];
  RatelKirchhofftau_sym_IsochoricOgden_AD(bulk, N, m, alpha, e_sym, tau_sym);
  for (int j = 0; j < 6; j++) printf("%f\t", tau_sym[j]);
  printf("\n\n");

  // dtau
  double dtau_sym[6];
  Rateldtau_fwd(bulk, N, m, alpha, e_sym, de_sym, tau_sym, dtau_sym);
  for (int j = 0; j < 6; j++) printf("%f\t", dtau_sym[j]);
  printf("\n\n");

  return 0;
}