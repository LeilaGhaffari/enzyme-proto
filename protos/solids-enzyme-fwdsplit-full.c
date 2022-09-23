// Standard Neo-Hookean formulation in initial configuration with Enzyme-AD
// S = dStrainEnergy/dE

#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>

#define VECSIZE 6

double __enzyme_autodiff(void*, ...);

int enzyme_dup;
int enzyme_const;

double StrainEnergy(double E_Voigt[VECSIZE], double mu, double lambda);

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

double log1p_series(double x) {
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

void S_analytical(double S_an[VECSIZE], double E_Voigt[VECSIZE], double mu, double lambda) {
  double E2_Voigt[VECSIZE];
  for(int i = 0; i<6; i++) E2_Voigt[i] = 2*E_Voigt[i];

  double E2[3][3];
  RatelVoigtUnpack(E2_Voigt, E2);
   
  // C : right Cauchy-Green tensor
  double C[3][3] = {{1 + E2[0][0], E2[0][1], E2[0][2]},
                    {E2[0][1], 1 + E2[1][1], E2[1][2]},
                    {E2[0][2], E2[1][2], 1 + E2[2][2]}
                   };

  // Compute C^(-1) : C-Inverse
  double Cinv_Voigt[VECSIZE];
  double C_inv[3][3];
  double detCm1 = RatelMatDetAM1(E2);
  double J = sqrt(detCm1 + 1);
  RatelMatComputeInverseSymmetric(C, detCm1+1, Cinv_Voigt);
  RatelVoigtUnpack(Cinv_Voigt, C_inv);

  // Compute the Second Piola-Kirchhoff (S)
  int indj[VECSIZE] = {0, 1, 2, 1, 0, 0}, indk[VECSIZE] = {0, 1, 2, 2, 2, 1};
  double logJ = log1p_series(detCm1) / 2.;
  for ( int m = 0; m < VECSIZE; m++) {
      S_an[m] = lambda*logJ*Cinv_Voigt[m];
      for ( int n = 0; n < 3; n++)
          S_an[m] += mu*C_inv[indj[m]][n]*E2[n][indk[m]];
  }
};

double StrainEnergy(double E_Voigt[VECSIZE], double mu, double lambda) {
  // Calculate 2*E 
  double E2_Voigt[VECSIZE];
  for(int i = 0; i<6; i++) E2_Voigt[i] = E_Voigt[i]*2; 

  // log(J)
  double E2[3][3];
  RatelVoigtUnpack(E2_Voigt, E2);
  double detCm1 = RatelMatDetAM1(E2);
  double J = sqrt(detCm1 + 1);
  double logJ = log1p_series(detCm1) / 2.;

  // trace(E)
  double traceE = (E_Voigt[0] + E_Voigt[1] + E_Voigt[2]);

  //return lambda*(J*J - 1)/4 - lambda*logJ/2  + mu * (-logJ + traceE);
  return lambda*logJ*logJ/2  + mu * (-logJ + traceE);
};

void S_autodiff(double *S_ad, double *E_Voigt, double mu, double lambda) {
  for (int i=0; i<VECSIZE; i++) S_ad[i] = 0.;
  __enzyme_autodiff((void *)StrainEnergy, 
                     E_Voigt, S_ad,
                     enzyme_const, mu,
                     enzyme_const, lambda);
  for (int i=VECSIZE/2; i<VECSIZE; i++) S_ad[i] /= 2.;
}

int main() {
  double mu = 1.;
  double lambda = 1.;

  double E_Voigt[VECSIZE];
  E_Voigt[0] = 0.5895232828911128/2;
  E_Voigt[1] = 0.2362491738162759/2;
  E_Voigt[2] = 0.9793730522395296/2;
  E_Voigt[3] = 0.2190993957421843/2;
  E_Voigt[4] = 0.0126503210747925/2;
  E_Voigt[5] = 0.6570956167695403/2;

  // Compute S with Enzyme-AD Forward mode
  double strain_energy = StrainEnergy(E_Voigt, mu, lambda);
  double S_ad[VECSIZE];
  S_autodiff(S_ad, E_Voigt, mu, lambda);

  // Compute analytical S
  double S_an[VECSIZE];
  S_analytical(S_an, E_Voigt, mu, lambda);

  printf("\n\nStrain Energy   = ");
  printf("\t   %.6lf", strain_energy);

  printf("\n\nS_autodiff      =\n\n");
  for (int i=0; i<VECSIZE; i++) printf("\t\t%.12lf", S_ad[i]);
  printf("\n\n");

  printf("\n\nS_analytical    =\n\n");
  for (int i=0; i<VECSIZE; i++) printf("\t\t%.12lf", S_an[i]);
  printf("\n\n");

  return 0;
}

/* Output:

Strain Energy   =          0.507029

S_autodiff      =

                0.629897288210          0.514638184498          0.763455742111          0.052445668610          -0.019798047937         0.200227118654



S_analytical    =

                1.389589119252          1.510916770389          1.248998632432          -0.055207003871         0.020840441891          -0.210769346783
*/
