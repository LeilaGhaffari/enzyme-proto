
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

void VoigtUnpack(const double sym[6], double full[3][3]) {
  full[0][0] = sym[0];
  full[0][1] = sym[5];
  full[0][2] = sym[4];
  full[1][0] = sym[5];
  full[1][1] = sym[1];
  full[1][2] = sym[3];
  full[2][0] = sym[4];
  full[2][1] = sym[3];
  full[2][2] = sym[2];
};

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

double MatDetAM1(const double A[3][3]) {
    return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) +
           A[0][1] * (A[1][2] * A[2][0] - A[1][0] * A[2][2]) +
           A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]) +
           A[0][0] + A[1][1] + A[2][2] +
           A[0][0] * A[1][1] + A[0][0] * A[2][2] + A[1][1] * A[2][2] -
           A[0][1] * A[1][0] - A[0][2] * A[2][0] - A[1][2] * A[2][1];
};

double VoigtDetAM1(const double V[6]) {
    return V[0] * (V[1] * V[2] - V[3] * V[3]) +
           V[5] * (V[3] * V[4] - V[5] * V[2]) +
           V[4] * (V[5] * V[3] - V[4] * V[1]) +
           V[0] + V[1] + V[2] +
           V[0] * V[1] + V[0] * V[2] + V[1] * V[2] -
           V[5] * V[5] - V[4] * V[4] - V[3] * V[3];
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

double VoigtTrace(double V[6]) { return V[0] + V[1] + V[2]; };

void SecondPiolaKirchhoffStress_NeoHookean_Analytical(double E_Voigt[6], double S_Voigt[6], const double lambda, const double mu) {
  double E2_Voigt[6];
  for(int i = 0; i<6; i++) E2_Voigt[i] = 2*E_Voigt[i];

  double E2[3][3];
  VoigtUnpack(E2_Voigt, E2);

  // Right Cauchy-Green tensor
  double C[3][3] = {{1 + E2[0][0], E2[0][1], E2[0][2]},
                     {E2[0][1], 1 + E2[1][1], E2[1][2]},
                     {E2[0][2], E2[1][2], 1 + E2[2][2]}
                    };

  // Compute C^(-1) : C-Inverse
  double Cinv_Voigt[6];
  double C_inv[3][3];
  double detCm1 = MatDetAM1(E2);
  MatComputeInverseSymmetric(C, detCm1+1, Cinv_Voigt);
  VoigtUnpack(Cinv_Voigt, C_inv);

  // Compute the Second Piola-Kirchhoff (S)
  int indj[6] = {0, 1, 2, 1, 0, 0}, indk[6] = {0, 1, 2, 2, 2, 1};
  double logJ = Log1pSeries(detCm1) / 2.;
  for ( int m = 0; m < 6; m++) {
      S_Voigt[m] = lambda*logJ*Cinv_Voigt[m];
      for ( int n = 0; n < 3; n++)
          S_Voigt[m] += mu*C_inv[indj[m]][n]*E2[n][indk[m]];
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

double *tau_from_gradPsi(double gradPsi_sym[6], double e_p[6]) {

  auto tau_sym = new double[n];
  double b_sym[6], gradPsi[3][3], b[3][3], tau[3][3];
  // b = 2 e + I
  for (int j = 0; j < 6; j++) b_sym[j] = 2. * e_p[j] + (j < 3);
  // tau = (dPsi / de) b
  SymmetricMatUnpack(gradPsi_sym, gradPsi);
  SymmetricMatUnpack(b_sym, b);
  MatMatMult(1., gradPsi, b, tau);
  SymmetricMatPack(tau, tau_sym);
  return tau_sym;
}
