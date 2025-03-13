#include "test036.h"

// Strain Energy
double RatelStrainEnergy_IsochoricOgdenCurrentAD_Enzyme(const double bulk, const int N, const double *m, const double *alpha, double e_sym[6]) {
  double e2_sym[6], e_vals[3], pr_str[3], pr_str_bar[3], V;

  // J
  for (int i = 0; i < 6; i++) e2_sym[i] = 2 * e_sym[i];
  const double detbm1 = RatelMatDetAM1Symmetric(e2_sym);
  const double J      = sqrt(detbm1 + 1);

  RatelMatComputeEigensystemSymmetric(e_sym, e_vals, NULL);
  RatelPrincipalStretch(e_vals, pr_str);
  const double J_pow = pow(J, -1. / 3.);

  RatelScalarVecMult(J_pow, pr_str, pr_str_bar);
  VolumetricFunctionAndDerivatives(J - 1., &V, NULL, NULL);

  // RatelStrainEnergy_IsochoricOgden(V, bulk, N, m, alpha, pr_str_bar, &strain_energy);
  // It worked by copying over the content of RatelStrainEnergy_IsochoricOgden()!!!
  double omega[3], energy_iso = 0.;
  // energy_iso = sum_{j=1:3} omega[j] = sum_{j=1:3} (sum_{k=1:N} m_k/alpha_k (pr_bar_j^alpha_k - 1))
  for (int j = 0; j < 3; j++) {
    omega[j] = 0;
    for (int k = 0; k < N; k++) {
      omega[j] += (m[k] / alpha[k]) * (pow(pr_str_bar[j], alpha[k]) - 1.);
    }
    energy_iso += omega[j];
  }

  return bulk * V + energy_iso;
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
void Rateldtau_fwd(const double bulk, const int N, const double *m, const double *alpha, double e_sym[6], double de_sym[6], double tau_sym[6],
                   double dtau_sym[6]) {
  __enzyme_fwddiff((void *)RatelKirchhofftau_sym_IsochoricOgden_AD, enzyme_const, bulk, enzyme_const, N, enzyme_const, m, enzyme_const, alpha, e_sym,
                   de_sym, tau_sym, dtau_sym);
}

int main() {
  double bulk = 1.;
  int    N    = 3;

  // Get N from the user
  printf("Enter N: ");
  scanf("%d", &N);

  // Allocate memory for m and alpha based on N
  double *m     = (double *)malloc(N * sizeof(double));
  double *alpha = (double *)malloc(N * sizeof(double));

  for (int i = 0; i < N; i++) {
    m[i]     = .1 * (i + 1.);
    alpha[i] = i + 1.;
  }

  // Constants for strain
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

  // Free allocated memory
  free(m);
  free(alpha);

  return 0;
}

/*
# N =3
Strain Energy = 0.421338

0.090270        0.685211        0.581716        0.856592        0.089502        0.124111

1.529920        1.113393        1.719482        -0.392328       0.119887        0.139680
*/
