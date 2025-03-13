#include <math.h>
#include <stdio.h>

extern int enzyme_const;
void       __enzyme_autodiff(void *, ...);
void       __enzyme_fwddiff(void *, ...);

double     func(const double *m, const double *alpha, const int N, double x[6]) {
  double omega[3], energy_iso = 0.;
  for (int j = 0; j < 3; j++) {
    omega[j] = 0;
    for (int k = 0; k < N; k++) {
      omega[j] += (m[k] / alpha[k]) * (pow(x[j], alpha[k]) - 1.);
    }
    energy_iso += omega[j];
  }
  return energy_iso;
}

// dfunc
void dfunc(const double *m, const double *alpha, const int N, double x[6], double df[6]) {
  __enzyme_autodiff((void *)func, enzyme_const, m, enzyme_const, alpha, enzyme_const, N, x, df);
  for (int i = 3; i < 6; i++) df[i] /= 2.;
}

// d2func
void d2func_fwd(const double *m, const double *alpha, const int N, double x[6], double dx[6], double df[6], double d2f[6]) {
  __enzyme_fwddiff((void *)dfunc, enzyme_const, m, enzyme_const, alpha, enzyme_const, N, x, dx, df, d2f);
}

int main() {
  // Constants
  int    N = 3;
  double m[N], alpha[N];
  for (int i = 0; i < N; i++) {
    m[i]     = .1 * (i + 1.);
    alpha[i] = (i + 1.);
  }

  double x[6]  = {0.0702417, 0.4799115, 0.3991242, 0.6756593, 0.0633284, 0.0959267};
  double dx[6] = {0.1425560, 0.115120, 0.551640, 0.0591922, 0.123535, 0.166572};

  // Strain energy
  double strain_energy = func(m, alpha, N, x);
  printf("\n\nStrain Energy = %f \n\n", strain_energy);

  // tau
  double df[6] = {0.};
  dfunc(m, alpha, N, x, df);
  for (int j = 0; j < 6; j++) printf("%f\t", df[j]);
  printf("\n\n");

  // dtau
  double d2f[6] = {0.};
  d2func_fwd(m, alpha, N, x, dx, df, d2f);
  for (int j = 0; j < 6; j++) printf("%f\t", d2f[j]);
  printf("\n\n");

  return 0;
}
