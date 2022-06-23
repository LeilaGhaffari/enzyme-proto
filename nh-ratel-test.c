// Standard Neo-Hookean formulation in initial configuration with Enzyme-AD
// S = dStrainEnergy/dE

#include <stdio.h>
#include <math.h>
#include <string.h>

#define VECSIZE 6

double __enzyme_autodiff(void*, ...);

int enzyme_dup;
int enzyme_const;

double StrainEnergy(double Ework[VECSIZE], double mu, double lambda);

double RatelLog1pSeriesShifted(double x) {
  const double left = sqrt(2.) / 2 - 1, right = sqrt(2.) - 1;
  double       sum = 0;
  // Shift first
  // Replace if with while for arbitrary range (may hurt vectorization)
  if (x < left) {
    sum -= log(2.) / 2;
    x = 1 + 2 * x;
  } else if (right < x) {
    sum += log(2.) / 2;
    x = (x - 1) / 2;
  }
  double       y  = x / (2. + x);
  const double y2 = y * y;
  sum += y;
  y *= y2;
  sum += y / 3;
  y *= y2;
  sum += y / 5;
  y *= y2;
  sum += y / 7;
  return 2 * sum;
}

double computeDetCM1(double *E2work) {
  return E2work[0]*(E2work[1]*E2work[2]-E2work[3]*E2work[3]) +
         E2work[5]*(E2work[4]*E2work[3]-E2work[5]*E2work[2]) +
         E2work[4]*(E2work[5]*E2work[3]-E2work[4]*E2work[1]) +
         E2work[0] + E2work[1] + E2work[2] +
         E2work[0]*E2work[1] + E2work[0]*E2work[2] +
         E2work[1]*E2work[2] - E2work[5]*E2work[5] -
         E2work[4]*E2work[4] - E2work[3]*E2work[3];
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

void S_analytical(double S_an[VECSIZE], double Ework[VECSIZE], double mu, double lambda) {
  
  double E[3][3] =  {{Ework[0], Ework[5], Ework[4]},
                     {Ework[5], Ework[1], Ework[3]},
                     {Ework[4], Ework[3], Ework[2]}
                    };
   
  double E2work[6];
  for(int i = 0; i<6; i++){
    E2work[i] = Ework[i]*2.;
  } 

  // C : right Cauchy-Green tensor
  double C[3][3]; 
  for(int i = 0; i<3; i++){
    for(int j = 0; j<3; j++){
        C[i][j] = 2*E[i][j];
        if(i == j) C[i][j] += 1;
    }
  }

  // Compute C^(-1) : C-Inverse
  double Cinvwork[VECSIZE];
  double detCm1 = computeDetCM1(E2work);
  computeMatinvSym(C, detCm1+1, Cinvwork);

  double C_inv[3][3] = {{Cinvwork[0], Cinvwork[5], Cinvwork[4]},
                        {Cinvwork[5], Cinvwork[1], Cinvwork[3]},
                        {Cinvwork[4], Cinvwork[3], Cinvwork[2]}
                       };

  // Compute the Second Piola-Kirchhoff (S)
  int indj[VECSIZE] = {0, 1, 2, 1, 0, 0}, indk[VECSIZE] = {0, 1, 2, 2, 2, 1};
  double logJ = RatelLog1pSeriesShifted(detCm1) / 2.;
  for ( int m = 0; m < VECSIZE; m++) {
      S_an[m] = lambda*logJ*Cinvwork[m];
      for ( int n = 0; n < 3; n++)
          S_an[m] += mu*C_inv[indj[m]][n]*E[n][indk[m]];
  }
};

double StrainEnergy(double Ework[VECSIZE], double mu, double lambda) {
   
  // log(J)
  double E2work[6];
  for(int i = 0; i<6; i++){
    E2work[i] = Ework[i]*2.;
  } 

  double detCm1 = computeDetCM1(E2work);
  double logJ = RatelLog1pSeriesShifted(detCm1)/2.;

  // trace(E)
  double traceE = (Ework[0] + Ework[1] + Ework[2]);

  return lambda*logJ*logJ/2. + mu * (-logJ + traceE);
};

void S_autodiff(double *S_ad, double *Ework, double mu, double lambda) {
  for (int i=0; i<VECSIZE; i++) S_ad[i] = 0.;
  __enzyme_autodiff((void *)StrainEnergy, 
                     Ework, S_ad,
                     enzyme_const, mu,
                     enzyme_const, lambda);
  // The first 3 entries (diagonal in Voigt notation) worked with 2E here, but
  // we need the gradient with respect to E. The next three (off-diagonal
  // components) represent entries that appear twice in the matrix, thus they
  // already have the necessary scaling.
  for (int i=0; i<3; i++)
    S_ad[i] *= 2.;
}

int main() {
  double mu = 1.;
  double lambda = 1.;

  double Ework[VECSIZE];
  Ework[0] = 0.5895232828911128;
  Ework[1] = 0.2362491738162759;
  Ework[2] = 0.9793730522395296;
  Ework[3] = 0.2190993957421843;
  Ework[4] = 0.0126503210747925;
  Ework[5] = 0.6570956167695403;

  // Compute S with Enzyme-AD Forward mode
  double strain_energy = StrainEnergy(Ework, mu, lambda);
  double S_ad[VECSIZE];
  S_autodiff(S_ad, Ework, mu, lambda);
  
  // Compute analytical S
  double S_an[VECSIZE];
  S_analytical(S_an, Ework, mu, lambda);

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
