#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define RATEL_EPSILON_DOUBLE 1e-16
#define RatelMax(a, b) ((a) > (b) ? (a) : (b))

int RatelScalarVecMult(double alpha, const double a[3], double b[3]) {
  b[0] = alpha * a[0];
  b[1] = alpha * a[1];
  b[2] = alpha * a[2];
  return 0;
}

double RatelDot3(const double a[3], const double b[3]) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

double RatelMatDetAM1Symmetric(const double A_sym[6]) {
  return A_sym[0] * (A_sym[1] * A_sym[2] - A_sym[3] * A_sym[3]) + A_sym[5] * (A_sym[3] * A_sym[4] - A_sym[5] * A_sym[2]) +
         A_sym[4] * (A_sym[5] * A_sym[3] - A_sym[4] * A_sym[1]) + A_sym[0] + A_sym[1] + A_sym[2] + A_sym[0] * A_sym[1] + A_sym[0] * A_sym[2] +
         A_sym[1] * A_sym[2] - A_sym[5] * A_sym[5] - A_sym[4] * A_sym[4] - A_sym[3] * A_sym[3];
}

int RatelSymmetricMatUnpack(const double sym[6], double full[3][3]) {
  full[0][0] = sym[0];
  full[0][1] = sym[5];
  full[0][2] = sym[4];
  full[1][0] = sym[5];
  full[1][1] = sym[1];
  full[1][2] = sym[3];
  full[2][0] = sym[4];
  full[2][1] = sym[3];
  full[2][2] = sym[2];
  return 0;
}

double RatelMatTraceSymmetric(const double A_sym[6]) { return A_sym[0] + A_sym[1] + A_sym[2]; }

int    RatelMatDeviatoricSymmetric(double trace_A, const double A_sym[6], double A_sym_dev[6]) {
  for (int j = 0; j < 6; j++) {
    A_sym_dev[j] = A_sym[j] - (j < 3) * (trace_A / 3.);
  }
  return 0;
}

int RatelMatMatMultSymmetric(double alpha, const double A_sym[6], const double B_sym[6], double C_sym[6]) {
  double    A[3][3], B[3][3];
  const int ind_j[6] = {0, 1, 2, 1, 0, 0}, ind_k[6] = {0, 1, 2, 2, 2, 1};

  RatelSymmetricMatUnpack(A_sym, A);
  RatelSymmetricMatUnpack(B_sym, B);
  for (int m = 0; m < 6; m++) {
    C_sym[m] = 0.;
    for (int n = 0; n < 3; n++) {
      C_sym[m] += A[ind_j[m]][n] * B[n][ind_k[m]];
    }
    C_sym[m] *= alpha;
  }
  return 0;
}

int RatelMatMatAddSymmetric(double alpha, const double A_sym[6], double beta, const double B_sym[6], double C_sym[6]) {
  for (int j = 0; j < 6; j++) {
    C_sym[j] = alpha * A_sym[j] + beta * B_sym[j];
  }
  return 0;
}

double RatelMatNorm(const double A[3][3]) {
  double sum_abs = 0.;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      sum_abs += A[i][j] * A[i][j];
    }
  }
  return sqrt(sum_abs);
}

int RatelSign(double x) { return (x > 0) - (x < 0); }

int RatelVecVecCross(const double a[3], const double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
  return 0;
}
double RatelAtanSeries(double x) {
  // Series approximation of `atan()` with max error 1.95815585968262e-14 on [-1,1]

  const double a_1  = 0.9999999999993525;
  const double a_3  = -0.33333333321533365;
  const double a_5  = 0.19999999355144213;
  const double a_7  = -0.14285697629655103;
  const double a_9  = 0.11110863734873448;
  const double a_11 = -0.09088554595195618;
  const double a_13 = 0.07676939640679648;
  const double a_15 = -0.06594663814817052;
  const double a_17 = 0.05632170111970379;
  const double a_19 = -0.046021855121716235;
  const double a_21 = 0.03405390225731899;
  const double a_23 = -0.021369751660910184;
  const double a_25 = 0.010577355252604226;
  const double a_27 = -0.003785491301089091;
  const double a_29 = 0.0008586183815558662;
  const double a_31 = -9.184922434344886e-05;
  const double x_sq = x * x;

  // Evaluate polynomial using Horner's method
  // clang-format off
    return x * (a_1 + x_sq * (a_3 + x_sq * (a_5 + x_sq * (a_7 + x_sq * (a_9 + x_sq * (a_11 + x_sq * (a_13 + x_sq * (a_15 + x_sq * (a_17 +
      x_sq * (a_19 + x_sq * (a_21 + x_sq * (a_23 + x_sq * (a_25 + x_sq * (a_27 + x_sq * (a_29 + x_sq * a_31)))))))))))))));
  // clang-format on
}

int RatelComputeEigenvector0(const double A_sym[6], double e_val_0, double e_vec_0[3]) {
  int    i_max = 0;
  double d_max;
  // rows of A - e_val_0 * I
  double row_0[3] = {A_sym[0] - e_val_0, A_sym[5], A_sym[4]};
  double row_1[3] = {A_sym[5], A_sym[1] - e_val_0, A_sym[3]};
  double row_2[3] = {A_sym[4], A_sym[3], A_sym[2] - e_val_0};
  double r0_x_r1[3], r0_x_r2[3], r1_x_r2[3], d_0, d_1, d_2;

  RatelVecVecCross(row_0, row_1, r0_x_r1);
  RatelVecVecCross(row_0, row_2, r0_x_r2);
  RatelVecVecCross(row_1, row_2, r1_x_r2);

  d_0   = RatelDot3(r0_x_r1, r0_x_r1);
  d_1   = RatelDot3(r0_x_r2, r0_x_r2);
  d_2   = RatelDot3(r1_x_r2, r1_x_r2);
  d_max = d_0;

  if (d_1 > d_max) {
    d_max = d_1;
    i_max = 1;
  }
  if (d_2 > d_max) {
    i_max = 2;
  }

  if (i_max == 0) RatelScalarVecMult(1.0 / sqrt(d_0), r0_x_r1, e_vec_0);
  else if (i_max == 1) RatelScalarVecMult(1.0 / sqrt(d_1), r0_x_r2, e_vec_0);
  else RatelScalarVecMult(1.0 / sqrt(d_2), r1_x_r2, e_vec_0);
  return 0;
}
#define FLOPS_ComputeEigenvector0 (4 + 3 * FLOPS_VecVecCross + 3 * FLOPS_Dot3 + FLOPS_ScalarVecMult)

int RatelOrthogonalComplement(const double w[3], double u[3], double v[3]) {
  double inv_length;

  // The vector w is guaranteed to be unit-length, so there is no need to worry about division by zero
  if (fabs(w[0]) > fabs(w[1])) {
    inv_length = 1. / sqrt(w[0] * w[0] + w[2] * w[2]);
    u[0]       = -w[2] * inv_length;
    u[1]       = 0.0;
    u[2]       = w[0] * inv_length;
  } else {
    inv_length = 1. / sqrt(w[1] * w[1] + w[2] * w[2]);
    u[0]       = 0.0;
    u[1]       = w[2] * inv_length;
    u[2]       = -w[1] * inv_length;
  }

  RatelVecVecCross(w, u, v);
  return 0;
}
int RatelMatVecMult(double alpha, const double A[3][3], const double x[3], double b[3]) {
  for (int j = 0; j < 3; j++) {
    b[j] = 0;
    for (int m = 0; m < 3; m++) {
      b[j] += alpha * A[j][m] * x[m];
    }
  }
  return 0;
}

int RatelVecVecSubtract(const double a[3], const double b[3], double c[3]) {
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
  return 0;
}

int RatelComputeEigenvector1(const double A_sym[6], const double e_vec_0[3], double e_val_1, double e_vec_1[3]) {
  double u[3], v[3], A[3][3], Au[3], Av[3];
  double d_00, d_11, m_00, m_01, m_11, max_abs_comp, mu[3], mv[3];

  RatelOrthogonalComplement(e_vec_0, u, v);
  RatelSymmetricMatUnpack(A_sym, A);
  RatelMatVecMult(1.0, A, u, Au);
  RatelMatVecMult(1.0, A, v, Av);

  d_00 = RatelDot3(u, Au);
  m_01 = RatelDot3(u, Av);
  d_11 = RatelDot3(v, Av);
  m_00 = d_00 - e_val_1;
  m_11 = d_11 - e_val_1;

  const double abs_m_00 = fabs(m_00), abs_m_01 = fabs(m_01), abs_m_11 = fabs(m_11);

  if (abs_m_00 >= abs_m_11) {
    max_abs_comp = RatelMax(abs_m_00, abs_m_01);
    if (max_abs_comp > 0.) {
      if (abs_m_00 >= abs_m_01) {
        m_01 /= m_00;
        m_00 = 1. / sqrt(1 + m_01 * m_01);
        m_01 *= m_00;
      } else {
        m_00 /= m_01;
        m_01 = 1. / sqrt(1 + m_00 * m_00);
        m_00 *= m_01;
      }

      RatelScalarVecMult(m_01, u, mu);
      RatelScalarVecMult(m_00, v, mv);
      RatelVecVecSubtract(mu, mv, e_vec_1);
    } else {
      e_vec_1[0] = u[0];
      e_vec_1[1] = u[1];
      e_vec_1[2] = u[2];
    }
  } else {
    max_abs_comp = RatelMax(abs_m_11, abs_m_01);
    if (max_abs_comp > 0.) {
      if (abs_m_11 >= abs_m_01) {
        m_01 /= m_11;
        m_11 = 1. / sqrt(1 + m_01 * m_01);
        m_01 *= m_11;
      } else {
        m_11 /= m_01;
        m_01 = 1. / sqrt(1 + m_11 * m_11);
        m_11 *= m_01;
      }

      RatelScalarVecMult(m_11, u, mu);
      RatelScalarVecMult(m_01, v, mv);
      RatelVecVecSubtract(mu, mv, e_vec_1);
    } else {
      e_vec_1[0] = u[0];
      e_vec_1[1] = u[1];
      e_vec_1[2] = u[2];
    }
  }
  return 0;
}

int RatelMatComputeEigensystemSymmetric(const double A_sym[6], double e_vals[3], double e_vecs[3][3]) {
  double I1, J2, s, b_00, b_11, b_22, off_diag;

  I1       = RatelMatTraceSymmetric(A_sym);
  b_00     = (A_sym[0] - A_sym[1]) * (A_sym[0] - A_sym[1]);
  b_11     = (A_sym[1] - A_sym[2]) * (A_sym[1] - A_sym[2]);
  b_22     = (A_sym[2] - A_sym[0]) * (A_sym[2] - A_sym[0]);
  J2       = (b_00 + b_11 + b_22) / 6. + (A_sym[3] * A_sym[3] + A_sym[4] * A_sym[4] + A_sym[5] * A_sym[5]);
  s        = sqrt(J2 / 3.);
  off_diag = A_sym[5] * A_sym[5] + A_sym[4] * A_sym[4] + A_sym[3] * A_sym[3];

  // zero matrix
  if (s < RATEL_EPSILON_DOUBLE && I1 == 0.) {
    e_vals[0] = 0.;
    e_vals[1] = 0.;
    e_vals[2] = 0.;
    if (e_vecs) {
      e_vecs[0][0] = 1., e_vecs[0][1] = 0., e_vecs[0][2] = 0.;
      e_vecs[1][0] = 0., e_vecs[1][1] = 1., e_vecs[1][2] = 0.;
      e_vecs[2][0] = 0., e_vecs[2][1] = 0., e_vecs[2][2] = 1.;
    }
    return 0;
  }

  // non-zero matrix
  if (off_diag > 0.) {
    double B[3][3], C[3][3], S_sym[6], S2_sym[6], T_sym[6], B_sym[6], C_sym[6];
    double d, sign_j;

    // S = A - I1/3 * I
    RatelMatDeviatoricSymmetric(I1, A_sym, S_sym);
    RatelMatMatMultSymmetric(1.0, S_sym, S_sym, S2_sym);
    // T = S*S - (2*J2/3) * I
    RatelMatDeviatoricSymmetric(2 * J2, S2_sym, T_sym);
    // B = T - s*S
    RatelMatMatAddSymmetric(1.0, T_sym, -s, S_sym, B_sym);
    // C = T + s*S
    RatelMatMatAddSymmetric(1.0, T_sym, s, S_sym, C_sym);
    RatelSymmetricMatUnpack(B_sym, B);
    RatelSymmetricMatUnpack(C_sym, C);
    const double norm_B = RatelMatNorm(B);
    const double norm_C = RatelMatNorm(C);

    d      = norm_B / norm_C;
    sign_j = RatelSign(1. - d);
    if (sign_j * (1. - d) < RATEL_EPSILON_DOUBLE) {
      e_vals[0] = -sqrt(J2) + I1 / 3.;
      e_vals[1] = 0. + I1 / 3.;
      e_vals[2] = sqrt(J2) + I1 / 3.;
    } else {
      double theta, cd, sd;

      theta     = 2 * RatelAtanSeries(sign_j > 0 ? d : (1.0 / d)) / 3.;
      cd        = sign_j * s * cos(theta);
      sd        = sqrt(J2) * sin(theta);
      e_vals[0] = -cd - sd + I1 / 3.;
      e_vals[1] = -cd + sd + I1 / 3.;
      e_vals[2] = 2 * cd + I1 / 3.;
    }
    if (e_vecs) {
      RatelComputeEigenvector0(A_sym, e_vals[2], e_vecs[2]);
      RatelComputeEigenvector1(A_sym, e_vecs[2], e_vals[1], e_vecs[1]);
      RatelVecVecCross(e_vecs[1], e_vecs[2], e_vecs[0]);
    }
  } else {  // A_sym is diagonal
    e_vals[0] = A_sym[0];
    e_vals[1] = A_sym[1];
    e_vals[2] = A_sym[2];
    if (e_vecs) {
      e_vecs[0][0] = 1., e_vecs[0][1] = 0., e_vecs[0][2] = 0.;
      e_vecs[1][0] = 0., e_vecs[1][1] = 1., e_vecs[1][2] = 0.;
      e_vecs[2][0] = 0., e_vecs[2][1] = 0., e_vecs[2][2] = 1.;
    }
  }
  return 0;
}

int RatelPrincipalStretch(const double e_vals[3], double pr_str[3]) {
  for (int i = 0; i < 3; i++) {
    pr_str[i] = sqrt(1 + 2 * e_vals[i]);
  }
  return 0;
}

double RatelLog1pSeries(double x) {
  double sum = 0;
  double y   = x / (2. + x);
  double y2  = y * y;

  sum += y;
  for (int i = 0; i < 5; i++) {
    y *= y2;
    sum += y / (2 * i + 3);
  }
  return 2 * sum;
}

int VolumetricFunctionAndDerivatives(double Jm1, double *V, double *J_dVdJ, double *J2_d2VdJ2) {
  const double Jp1  = Jm1 + 2.;
  const double logJ = RatelLog1pSeries(Jm1);
  const double A    = Jm1 * Jp1 - 2. * logJ;

  if (V) *V = A * 0.25;
  if (J_dVdJ) *J_dVdJ = Jm1 * Jp1 * 0.5;
  if (J2_d2VdJ2) {
    const double J = Jm1 + 1.;

    *J2_d2VdJ2 = (J * J + 1) * 0.5;
  }
  return 0;
}

int RatelStrainEnergy_IsochoricOgden(double V, double bulk, int N, const double *m, const double *alpha, const double pr_str_bar[3],
                                     double *strain_energy) {
  double omega[3], energy_iso = 0.;

  // energy_iso = sum_{j=1:3} omega[j] = sum_{j=1:3} (sum_{k=1:N} m_k/alpha_k (pr_bar_j^alpha_k - 1))
  for (int j = 0; j < 3; j++) {
    omega[j] = 0;
    for (int k = 0; k < N; k++) {
      omega[j] += (m[k] / alpha[k]) * (pow(pr_str_bar[j], alpha[k]) - 1.);
    }
    energy_iso += omega[j];
  }

  // Strain energy psi(e) for Ogden
  *strain_energy = bulk * V + energy_iso;
  return 0;
}

int RatelMatMatMult(double alpha, const double A[3][3], const double B[3][3], double C[3][3]) {
  for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {
      C[j][k] = 0;
      for (int m = 0; m < 3; m++) {
        C[j][k] += alpha * A[j][m] * B[m][k];
      }
    }
  }
  return 0;
}

int RatelSymmetricMatPack(const double full[3][3], double sym[6]) {
  sym[0] = full[0][0];
  sym[1] = full[1][1];
  sym[2] = full[2][2];
  sym[3] = full[1][2];
  sym[4] = full[0][2];
  sym[5] = full[0][1];
  return 0;
}
