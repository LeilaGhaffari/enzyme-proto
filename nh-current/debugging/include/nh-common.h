
int SymmetricMatUnpack(const double sym[6], double full[3][3]) {
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

int SymmetricMatPack(const double full[3][3], double sym[6]) {
  sym[0] = full[0][0];
  sym[1] = full[1][1];
  sym[2] = full[2][2];
  sym[3] = full[1][2];
  sym[4] = full[0][2];
  sym[5] = full[0][1];
  return 0;
}

int MatInverse(const double A[3][3], double det_A, double A_inv[3][3]) {
  double B[3][3] = {
      {A[1][1] * A[2][2] - A[1][2] * A[2][1], A[0][2] * A[2][1] - A[0][1] * A[2][2], A[0][1] * A[1][2] - A[0][2] * A[1][1]},
      {A[1][2] * A[2][0] - A[1][0] * A[2][2], A[0][0] * A[2][2] - A[0][2] * A[2][0], A[0][2] * A[1][0] - A[0][0] * A[1][2]},
      {A[1][0] * A[2][1] - A[1][1] * A[2][0], A[0][1] * A[2][0] - A[0][0] * A[2][1], A[0][0] * A[1][1] - A[0][1] * A[1][0]},
  };

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A_inv[i][j] = B[i][j] / (det_A);
    }
  }
  return 0;
}

int MatMatMult(double alpha, const double A[3][3], const double B[3][3], double C[3][3]) {
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

int MatComputeInverseSymmetric(const double A[3][3], const double det_A, double A_inv[6]) {
  // Compute A^(-1) : A-Inverse
  double B[6] = {
      A[1][1] * A[2][2] - A[1][2] * A[2][1],
      A[0][0] * A[2][2] - A[0][2] * A[2][0],
      A[0][0] * A[1][1] - A[0][1] * A[1][0],
      A[0][2] * A[1][0] - A[0][0] * A[1][2],
      A[0][1] * A[1][2] - A[0][2] * A[1][1],
      A[0][2] * A[2][1] - A[0][1] * A[2][2]
  };
  for (int m = 0; m < 6; m++) {
    A_inv[m] = B[m] / (det_A);
  }
  return 0;
};

double MatDetAM1Symmetric(const double A_sym[6]) {
  return A_sym[0] * (A_sym[1] * A_sym[2] - A_sym[3] * A_sym[3]) +          /* *NOPAD* */
         A_sym[5] * (A_sym[3] * A_sym[4] - A_sym[5] * A_sym[2]) +          /* *NOPAD* */
         A_sym[4] * (A_sym[5] * A_sym[3] - A_sym[4] * A_sym[1]) +          /* *NOPAD* */
         A_sym[0] + A_sym[1] + A_sym[2] +                                  /* *NOPAD* */
         A_sym[0] * A_sym[1] + A_sym[0] * A_sym[2] + A_sym[1] * A_sym[2] - /* *NOPAD* */
         A_sym[5] * A_sym[5] - A_sym[4] * A_sym[4] - A_sym[3] * A_sym[3];  /* *NOPAD* */
};

double MatTraceSymmetric(const double A_sym[6]) { return A_sym[0] + A_sym[1] + A_sym[2]; }


double MatDetAM1(const double A[3][3]) {
    return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) +
           A[0][1] * (A[1][2] * A[2][0] - A[1][0] * A[2][2]) +
           A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]) +
           A[0][0] + A[1][1] + A[2][2] +
           A[0][0] * A[1][1] + A[0][0] * A[2][2] + A[1][1] * A[2][2] -
           A[0][1] * A[1][0] - A[0][2] * A[2][0] - A[1][2] * A[2][1];
};

double Log1pSeries(double x) {
    double sum = 0;
    double y = x / (2. + x);
    double y2 = y*y;
    sum += y;
    for (int i=0; i<5; i++) {
      y *= y2;
      sum += y / (2*i + 3);
    }
    return 2 * sum;
};

void SecondPiolaKirchhoffStress_NeoHookean_Analytical(double E_sym[6], double S_sym[6], const double lambda, const double mu) {
  double E2_sym[6];
  for(int i = 0; i<6; i++) E2_sym[i] = 2*E_sym[i];

  double E2[3][3];
  SymmetricMatUnpack(E2_sym, E2);

  // Right Cauchy-Green tensor
  double C[3][3] = {{1 + E2[0][0], E2[0][1], E2[0][2]},
                     {E2[0][1], 1 + E2[1][1], E2[1][2]},
                     {E2[0][2], E2[1][2], 1 + E2[2][2]}
                    };

  // Compute C^(-1) : C-Inverse
  double Cinv_sym[6];
  double C_inv[3][3];
  double detCm1 = MatDetAM1(E2);
  MatComputeInverseSymmetric(C, detCm1+1, Cinv_sym);
  SymmetricMatUnpack(Cinv_sym, C_inv);

  // Compute the Second Piola-Kirchhoff (S)
  int indj[6] = {0, 1, 2, 1, 0, 0}, indk[6] = {0, 1, 2, 2, 2, 1};
  double logJ = Log1pSeries(detCm1) / 2.;
  for ( int m = 0; m < 6; m++) {
      S_sym[m] = lambda*logJ*Cinv_sym[m];
      for ( int n = 0; n < 3; n++)
          S_sym[m] += mu*C_inv[indj[m]][n]*E2[n][indk[m]];
  }
};

int LinearStrain(const double grad_u[3][3], double e_sym[6]) {
  e_sym[0] = grad_u[0][0];
  e_sym[1] = grad_u[1][1];
  e_sym[2] = grad_u[2][2];
  e_sym[3] = (grad_u[1][2] + grad_u[2][1]) / 2.;
  e_sym[4] = (grad_u[0][2] + grad_u[2][0]) / 2.;
  e_sym[5] = (grad_u[0][1] + grad_u[1][0]) / 2.;
  return 0;
}

int GreenEulerStrain(const double grad_u[3][3], double e_sym[6]) {
  const int ind_j[6] = {0, 1, 2, 1, 0, 0}, ind_k[6] = {0, 1, 2, 2, 2, 1};

  LinearStrain(grad_u, e_sym);
  // Add (grad_u * grad_u^T)/2 term to the linear part of e_sym
  for (int m = 0; m < 6; m++) {
    for (int n = 0; n < 3; n++) {
      e_sym[m] += (grad_u[ind_j[m]][n] * grad_u[ind_k[m]][n]) * 0.5;
    }
  }
  return 0;
}

int GreenEulerStrain_fwd(const double grad_du[3][3], const double b[3][3], double de_sym[6]) {
  const int ind_j[6] = {0, 1, 2, 1, 0, 0}, ind_k[6] = {0, 1, 2, 2, 2, 1};

  for (int m = 0; m < 6; m++) {
    de_sym[m] = 0;
    for (int n = 0; n < 3; n++) {
      de_sym[m] += (grad_du[ind_j[m]][n] * b[n][ind_k[m]] + b[ind_j[m]][n] * grad_du[ind_k[m]][n]) / 2.;
    }
  }
  return 0;
}

double StrainEnergy(double e_sym[6], const double lambda, const double mu) {
  double e2_sym[6];

  // J and log(J)
  for (int i = 0; i < 6; i++) e2_sym[i] = 2 * e_sym[i];
  const double detbm1 = MatDetAM1Symmetric(e2_sym);
  const double J      = sqrt(detbm1 + 1);
  const double logJ   = Log1pSeries(detbm1) / 2.;

  // trace(e)
  const double trace_e = MatTraceSymmetric(e_sym);

  return lambda * (J * J - 1) / 4 - lambda * logJ / 2 + mu * (-logJ + trace_e);
};

int MatMatTransposeMult(const double A[3][3], const double B[3][3], double C[3][3]) {
  for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {
      C[j][k] = 0;
      for (int m = 0; m < 3; m++) {
        C[j][k] += A[j][m] * B[k][m];
      }
    }
  }
  return 0;
};

void PullBack_symmetric(double Grad_u[3][3], double a_sym[6], double A_sym[6]) {
  // F = I + Grad_u
  const double F[3][3] = {
    {Grad_u[0][0] + 1, Grad_u[0][1],     Grad_u[0][2]    },
    {Grad_u[1][0],     Grad_u[1][1] + 1, Grad_u[1][2]    },
    {Grad_u[2][0],     Grad_u[2][1],     Grad_u[2][2] + 1}
  };
  // F^{-1}
  double F_inv[3][3];
  const double Jm1  = MatDetAM1(Grad_u);
  const double detF = Jm1 + 1.;
  MatInverse(F, detF, F_inv);
  // A = F_inv * a * F_inv_T (pull-back)
  double F_inv_a[3][3], a[3][3], A[3][3];
  SymmetricMatUnpack(a_sym, a);
  MatMatMult(1., F_inv, a, F_inv_a);
  MatMatTransposeMult(F_inv_a, F_inv, A);
  SymmetricMatPack(A, A_sym);
};

void dPullBack_symmetric(double F_inv[3][3], double dF[3][3], double a_sym[6], double da_sym[6], double dA_sym[6]) {
  // F = I + Grad_u => dF = Grad_du
  // F_inv * F = I => dF_inv * F + F_inv * dF = 0 => dF_inv = -F_inv * dF * F_inv
  // A = F_inv * a * F_inv^T => da = F_inv da F_inv^T + dF_inv a F_inv^T + F_inv a dF_inv^T
  double a[3][3], da[3][3];
  SymmetricMatUnpack(da_sym, da);
  SymmetricMatUnpack(a_sym, a);

  // dF_inv = -F_inv * dF * F_inv
  double F_inv_dF[3][3], dF_inv[3][3];
  MatMatMult(-1., F_inv, dF, F_inv_dF);
  MatMatMult(1., F_inv_dF, F_inv, dF_inv);

  // F_inv da F_inv^T
  double F_inv_da[3][3], F_inv_da_F_inv_T[3][3];
  MatMatMult(1., F_inv, da, F_inv_da);
  MatMatTransposeMult(F_inv_da, F_inv, F_inv_da_F_inv_T);

  // dF_inv a F_inv^T
  double dF_inv_a[3][3], dF_inv_a_F_inv_T[3][3];
  MatMatMult(1., dF_inv, a, dF_inv_a);
  MatMatTransposeMult(dF_inv_a, F_inv, dF_inv_a_F_inv_T);

  // F_inv a dF_inv^T
  double F_inv_a[3][3], F_inv_a_dF_inv_T[3][3];
  MatMatMult(1., F_inv, a, F_inv_a);
  MatMatTransposeMult(F_inv_a, dF_inv, F_inv_a_dF_inv_T);

  // da = F_inv da F_inv^T + dF_inv a F_inv^T + F_inv a dF_inv^T
  double dA[3][3] = {{0.}};
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      dA[i][j] = F_inv_da_F_inv_T[i][j] + dF_inv_a_F_inv_T[i][j] + F_inv_a_dF_inv_T[i][j];

  SymmetricMatPack(dA, dA_sym);
};

void PushForward_symmetric(double Grad_u[3][3], double A_sym[6], double a_sym[6]) {
  // F = I + Grad_u
  const double F[3][3] = {
    {Grad_u[0][0] + 1, Grad_u[0][1],     Grad_u[0][2]    },
    {Grad_u[1][0],     Grad_u[1][1] + 1, Grad_u[1][2]    },
    {Grad_u[2][0],     Grad_u[2][1],     Grad_u[2][2] + 1}
  };
  // a = F * A * F^T (push-forward)
  double F_A[3][3], a[3][3], A[3][3];
  SymmetricMatUnpack(A_sym, A);
  MatMatMult(1., F, A, F_A);
  MatMatTransposeMult(F_A, F, a);
  SymmetricMatPack(a, a_sym);
};

void dPushForward_symmetric(const double F[3][3], double dF[3][3], double A_sym[6], double dA_sym[6], double da_sym[6]) {
  // F = I + Grad_u => dF = Grad_du
  // a = F * A * F^T => da = F dA F^T + dF A F^T + F A dF^T
  double A[3][3], dA[3][3];
  SymmetricMatUnpack(dA_sym, dA);
  SymmetricMatUnpack(A_sym, A);

  // F dA F^T
  double F_dA[3][3], F_dA_FT[3][3];
  MatMatMult(1., F, dA, F_dA);
  MatMatTransposeMult(F_dA, F, F_dA_FT);

  // dF A F^T
  double dF_A[3][3], dF_A_FT[3][3];
  MatMatMult(1., dF, A, dF_A);
  MatMatTransposeMult(dF_A, F, dF_A_FT);

  // F A dF^T
  double F_A[3][3], F_A_dFT[3][3];
  MatMatMult(1., F, A, F_A);
  MatMatTransposeMult(F_A, dF, F_A_dFT);

  // da = F dA F^T + dF A F^T + F A dF^T
  double da[3][3] = {{0.}};
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      da[i][j] = F_dA_FT[i][j] + dF_A_FT[i][j] + F_A_dFT[i][j];

  SymmetricMatPack(da, da_sym);
};

void MatTransposeMatMult(double alpha, const double A[3][3], const double B[3][3], double C[3][3]) {
  for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {
      C[j][k] = 0;
      for (int m = 0; m < 3; m++) {
        C[j][k] += alpha * A[m][j] * B[m][k];
      }
    }
  }
}

// E = .5(F^t F - I) not the pull-back operator
void Compute_E_symmetric(double Grad_u[3][3], double E_sym[6]) {
  // F = I + Grad_u
  const double F[3][3] = {
    {Grad_u[0][0] + 1, Grad_u[0][1],     Grad_u[0][2]    },
    {Grad_u[1][0],     Grad_u[1][1] + 1, Grad_u[1][2]    },
    {Grad_u[2][0],     Grad_u[2][1],     Grad_u[2][2] + 1}
  };
  double E[3][3];
  MatTransposeMatMult(1., F, F, E);
  for (int i=0; i<3; i++) E[i][i] -= 1.;
  for (int i=0; i<3; i++) for (int j=0; j<3; j++) E[i][j] /= 2.;
  SymmetricMatPack(E, E_sym);
}

void RatelGreenLagrangeStrain_fwd(const double grad_du[3][3], const double F[3][3], double dE_sym[6]) {
  const int ind_j[6] = {0, 1, 2, 1, 0, 0}, ind_k[6] = {0, 1, 2, 2, 2, 1};

  for (int m = 0; m < 6; m++) {
    dE_sym[m] = 0;
    for (int n = 0; n < 3; n++) {
      dE_sym[m] += (grad_du[n][ind_j[m]] * F[n][ind_k[m]] + F[n][ind_j[m]] * grad_du[n][ind_k[m]]) / 2.;
    }
  }
};

// `C = alpha A + beta B` for 3x3 matrices
int MatMatAdd(double alpha, const double A[3][3], double beta, const double B[3][3], double C[3][3]) {
  for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {
      C[j][k] = alpha * A[j][k] + beta * B[j][k];
    }
  }
  return 0;
}
