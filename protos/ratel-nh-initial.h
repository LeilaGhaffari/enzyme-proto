static const int FSInitialNH_AD_TAPE_SIZE = 6;

void RatelVoigtUnpack(const double sym[6], double full[3][3]) {
  full[0][0] = sym[0];
  full[0][1] = sym[5];
  full[0][2] = sym[4];
  full[1][0] = sym[5];
  full[1][1] = sym[1];
  full[1][2] = sym[3];
  full[2][0] = sym[4];
  full[2][1] = sym[3];
  full[2][2] = sym[2];
}

int RatelMatComputeInverseSymmetric(const double A[3][3], const double det_A, double A_inv[6]) {
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
}

double RatelMatDetAM1(const double A[3][3]) {
    return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) +  
           A[0][1] * (A[1][2] * A[2][0] - A[1][0] * A[2][2]) +        
           A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]) +        
           A[0][0] + A[1][1] + A[2][2] +                              
           A[0][0] * A[1][1] + A[0][0] * A[2][2] + A[1][1] * A[2][2] -
           A[0][1] * A[1][0] - A[0][2] * A[2][0] - A[1][2] * A[2][1]; 
}

double RatelVoigtDetAM1(const double V[6]) {
    return V[0] * (V[1] * V[2] - V[3] * V[3]) +     
           V[5] * (V[3] * V[4] - V[5] * V[2]) +     
           V[4] * (V[5] * V[3] - V[4] * V[1]) +     
           V[0] + V[1] + V[2] +
           V[0] * V[1] + V[0] * V[2] + V[1] * V[2] -
           V[5] * V[5] - V[4] * V[4] - V[3] * V[3]; 
}

double RatelLog1pSeries(double x) {
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

double RatelVoigtTrace(double V[6]) { return V[0] + V[1] + V[2]; }

void S_analytical(double S_an[6], double E_Voigt[6], double mu, double lambda) {
  double E2_Voigt[6];
  for(int i = 0; i<6; i++) E2_Voigt[i] = 2*E_Voigt[i];

  double E2[3][3];
  RatelVoigtUnpack(E2_Voigt, E2);
   
  // C : right Cauchy-Green tensor
  double C[3][3] = {{1 + E2[0][0], E2[0][1], E2[0][2]},
                    {E2[0][1], 1 + E2[1][1], E2[1][2]},
                    {E2[0][2], E2[1][2], 1 + E2[2][2]}
                   };

  // Compute C^(-1) : C-Inverse
  double Cinv_Voigt[6];
  double C_inv[3][3];
  double detCm1 = RatelMatDetAM1(E2);
  double J = sqrt(detCm1 + 1);
  RatelMatComputeInverseSymmetric(C, detCm1+1, Cinv_Voigt);
  RatelVoigtUnpack(Cinv_Voigt, C_inv);

  // Compute the Second Piola-Kirchhoff (S)
  int indj[6] = {0, 1, 2, 1, 0, 0}, indk[6] = {0, 1, 2, 2, 2, 1};
  double logJ = RatelLog1pSeries(detCm1) / 2.;
  for ( int m = 0; m < FSInitialNH_AD_TAPE_SIZE; m++) {
      S_an[m] = lambda*logJ*Cinv_Voigt[m];
      for ( int n = 0; n < 3; n++)
          S_an[m] += mu*C_inv[indj[m]][n]*E2[n][indk[m]];
  }
};

void GreenLagrangeStrain(const double grad_u[3][3], double E[6]) {
  const int ind_j[6] = {0, 1, 2, 1, 0, 0}, ind_k[6] = {0, 1, 2, 2, 2, 1};
  for (int m = 0; m < 6; m++) {
    E[m] = (grad_u[ind_j[m]][ind_k[m]] + grad_u[ind_k[m]][ind_j[m]]) * 0.5;
    for (int n = 0; n < 3; n++) {
      E[m] += (grad_u[n][ind_j[m]] * grad_u[n][ind_k[m]]) * 0.5;
    }
  }
};

void ComputeCinverse(double E_Voigt[6], const double Jm1, double C_inv_Voigt[6]) {
  double E[3][3];
  RatelVoigtUnpack(E_Voigt, E);

  // C : right Cauchy-Green tensor
  // C = I + 2E
  const double C[3][3] = {
      {2 * E[0][0] + 1, 2 * E[0][1],     2 * E[0][2]    },
      {2 * E[0][1],     2 * E[1][1] + 1, 2 * E[1][2]    },
      {2 * E[0][2],     2 * E[1][2],     2 * E[2][2] + 1}
  };

  // Compute C^(-1) : C-Inverse
  const double detC = (Jm1 + 1.) * (Jm1 + 1.);
  RatelMatComputeInverseSymmetric(C, detC, C_inv_Voigt);
}

int RatelMatMatMult(const double alpha, const double A[3][3], const double B[3][3], double C[3][3]) {
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
