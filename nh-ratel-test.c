// Standard Neo-Hookean formulation in initial configuration with Enzyme-AD
// S = dStrainEnergy/dE

#include <stdio.h>
#include <math.h>
#include <string.h>

#define VECSIZE 6

double __enzyme_autodiff(void*, ...);

int enzyme_dup;
int enzyme_const;

double StrainEnergy(double E_Voigt[VECSIZE], double mu, double lambda);

double log1p_series_shifted(double x) {
  double left = sqrt(2.)/2 - 1;
  double right = sqrt(2.) - 1;
  double sum = 0;
  if (x < left) {
      sum -= log(2.) / 2;
      x = 1 + 2 * x;
  } else if (right < x) {
      sum += log(2.) / 2;
      x = (x - 1) / 2;
  }
  double y = x / (2. + x);
  double y2 = y*y;
  sum += y;
  y *= y2;
  sum += y / 3;
  y *= y2;
  sum += y / 5;
  y *= y2;
  sum += y / 7;
  return 2 * sum;
};

double computeDetCM1(double *E2_Voigt) {
  return E2_Voigt[0]*(E2_Voigt[1]*E2_Voigt[2]-E2_Voigt[3]*E2_Voigt[3]) +
         E2_Voigt[5]*(E2_Voigt[4]*E2_Voigt[3]-E2_Voigt[5]*E2_Voigt[2]) +
         E2_Voigt[4]*(E2_Voigt[5]*E2_Voigt[3]-E2_Voigt[4]*E2_Voigt[1]) +
         E2_Voigt[0] + E2_Voigt[1] + E2_Voigt[2] +
         E2_Voigt[0]*E2_Voigt[1] + E2_Voigt[0]*E2_Voigt[2] +
         E2_Voigt[1]*E2_Voigt[2] - E2_Voigt[5]*E2_Voigt[5] -
         E2_Voigt[4]*E2_Voigt[4] - E2_Voigt[3]*E2_Voigt[3];
};

int computeMatinvSym(double A[3][3], double detA, double *Ainv) {
  // Compute A^(-1) : A-Inverse
  double B[VECSIZE] = {A[1][1]*A[2][2] - A[1][2]*A[2][1], 
                 A[0][0]*A[2][2] - A[0][2]*A[2][0], 
                 A[0][0]*A[1][1] - A[0][1]*A[1][0], 
                 A[0][2]*A[1][0] - A[0][0]*A[1][2], 
                 A[0][1]*A[1][2] - A[0][2]*A[1][1], 
                 A[0][2]*A[2][1] - A[0][1]*A[2][2]  
                };               
  for ( int m = 0; m < VECSIZE; m++) Ainv[m] = B[m] / (detA);
  return 0;
};

void S_analytical(double S_an[VECSIZE], double E_Voigt[VECSIZE], double mu, double lambda) {
  double E2_Voigt[VECSIZE];
  for(int i = 0; i<6; i++) E2_Voigt[i] = 2*E_Voigt[i];

  double E2[3][3] = {{E2_Voigt[0], E2_Voigt[5], E2_Voigt[4]},
                     {E2_Voigt[5], E2_Voigt[1], E2_Voigt[3]},
                     {E2_Voigt[4], E2_Voigt[3], E2_Voigt[2]}
                    };
   
  // C : right Cauchy-Green tensor
  double C[3][3] = {{1 + E2[0][0], E2[0][1], E2[0][2]},
                    {E2[0][1], 1 + E2[1][1], E2[1][2]},
                    {E2[0][2], E2[1][2], 1 + E2[2][2]}
                   };

  // Compute C^(-1) : C-Inverse
  double Cinvwork[VECSIZE];
  double detCm1 = computeDetCM1(E2_Voigt);
  computeMatinvSym(C, detCm1+1, Cinvwork);

  double C_inv[3][3] = {{Cinvwork[0], Cinvwork[5], Cinvwork[4]},
                        {Cinvwork[5], Cinvwork[1], Cinvwork[3]},
                        {Cinvwork[4], Cinvwork[3], Cinvwork[2]}
                       };

  // Compute the Second Piola-Kirchhoff (S)
  int indj[VECSIZE] = {0, 1, 2, 1, 0, 0}, indk[VECSIZE] = {0, 1, 2, 2, 2, 1};
  double logJ = log1p_series_shifted(detCm1) / 2.;
  for ( int m = 0; m < VECSIZE; m++) {
      S_an[m] = lambda*logJ*Cinvwork[m];;
      for ( int n = 0; n < 3; n++)
          S_an[m] += mu*C_inv[indj[m]][n]*E2[n][indk[m]];
  }
};

double StrainEnergy(double E_Voigt[VECSIZE], double mu, double lambda) {
  // Calculate 2*E 
  double E2_Voigt[VECSIZE];
  for(int i = 0; i<6; i++) E2_Voigt[i] = E_Voigt[i]*2; 

  // log(J)
  double detCm1 = computeDetCM1(E2_Voigt);
  double logJ = log1p_series_shifted(detCm1) / 2.;

  // trace(E)
  double traceE = (E_Voigt[0] + E_Voigt[1] + E_Voigt[2]);

  return lambda*logJ*logJ/2. + mu * (-logJ + traceE);
};

void S_autodiff(double *S_ad, double *E_Voigt, double mu, double lambda) {
  for (int i=0; i<VECSIZE; i++) S_ad[i] = 0.;
  __enzyme_autodiff((void *)StrainEnergy, 
                     E_Voigt, S_ad,
                     enzyme_const, mu,
                     enzyme_const, lambda);
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
