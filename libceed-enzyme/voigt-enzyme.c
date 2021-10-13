// Mocking  the function computeS() in 
// "libCEED/examples/solids/qfunctions/finite-strain-neo-hookean-initial-1.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define VECSIZE 6


int  __enzyme_augmentsize(void *, ...);
void __enzyme_augmentfwd(void *, ...);
void __enzyme_reverse(void *, ...);

int enzyme_dup;
int enzyme_tape;
int enzyme_const;
int enzyme_allocated;
void computeS(double *Swork, double *E2work, double lambda, double mu);
void grad_S(double *S, double *dS, double *E, double *dE, double lambda, double mu) {
    int size = __enzyme_augmentsize((void *)computeS, enzyme_dup, enzyme_dup, enzyme_const, enzyme_const); // We need to add Const as well
    void *data = malloc(size);
    __enzyme_augmentfwd((void *)computeS, enzyme_allocated, size, enzyme_tape, data, S, dS, E, dE, lambda, mu);
    __enzyme_reverse((void *)computeS,    enzyme_allocated, size, enzyme_tape, data, S, dS, E, dE, lambda, mu);
    free(data);
}


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

void computeS(double *Swork, double *E2work, double lambda, double mu) {
  
  double E2[3][3] = {{E2work[0], E2work[5], E2work[4]},
                     {E2work[5], E2work[1], E2work[3]},
                     {E2work[4], E2work[3], E2work[2]}
                    };
   
  // C : right Cauchy-Green tensor
  double C[3][3] = {{1 + E2[0][0], E2[0][1], E2[0][2]},
                    {E2[0][1], 1 + E2[1][1], E2[1][2]},
                    {E2[0][2], E2[1][2], 1 + E2[2][2]}
                   };

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
  double logJ = log1p_series_shifted(detCm1) / 2.;
  for ( int m = 0; m < VECSIZE; m++) {
      Swork[m] = lambda*logJ*Cinvwork[m];
      for ( int n = 0; n < 3; n++)
          Swork[m] += mu*C_inv[indj[m]][n]*E2[n][indk[m]];
  }
};

int main() {
  double E = .3;
  double nu = .3;
  double TwoMu = E / (1 + nu);
  double mu = TwoMu / 2;
  double Kbulk = E / (3*(1 - 2*nu));
  double lambda = (3*Kbulk - TwoMu) / 3;

  double E2work[VECSIZE] = {0., 0., 0., 0., 0., 0.};
  E2work[0] = 0.5895232828911128;
  E2work[1] = 0.2362491738162759;
  E2work[2] = 0.9793730522395296;
  E2work[3] = 0.2190993957421843;
  E2work[4] = 0.0126503210747925;
  E2work[5] = 0.6570956167695403;

  double J[VECSIZE][VECSIZE];
  for (int i=0; i<VECSIZE; i++) for (int j=0; j<6; j++) J[i][j] = 0.;

  double Swork[VECSIZE];
  for (int i=0; i<VECSIZE; i++) {
      double dSwork[VECSIZE]  = {0., 0., 0., 0., 0., 0.}; dSwork[i] = 1.;
      grad_S(Swork, dSwork, E2work, J[i], lambda, mu);
  }

  double deltaEwork[VECSIZE] = {0., 0., 0., 0., 0., 0.};
  deltaEwork[0] = 0.9681576729097205;
  deltaEwork[1] = 0.7994338113484318;
  deltaEwork[2] = 0.2755183472001872;
  deltaEwork[3] = 0.6500440500146469;
  deltaEwork[4] = 0.0593948875992271;
  deltaEwork[5] = 0.6002528007029311;

  double deltaS_[VECSIZE];
  for (int i=0; i<VECSIZE; i++) {
    deltaS_[i] = 0;
    for (int j=0; j<VECSIZE; j++)
      deltaS_[i] += J[i][j] * 2. * deltaEwork[j]; // we need deltaE2work = 2 .* deltaEwork 
  }

  double deltaS[3][3] = {{deltaS_[0], deltaS_[5], deltaS_[4]},
                         {deltaS_[5], deltaS_[1], deltaS_[3]},
                         {deltaS_[4], deltaS_[3], deltaS_[2]}
                        };

  printf("\n\nSwork       = ");
  for (int i=0; i<VECSIZE; i++) printf("\t   %.6lf", Swork[i]);

  printf("\n\ndeltaSwork  = ");
  for (int i=0; i<VECSIZE; i++) printf("\t   %.6lf", deltaS_[i]);

  printf("\n\ndeltaS      =\n\n");
  for (int i=0; i<3; i++) printf("\t\t%.12lf", deltaS[0][i]);
  printf("\n\n");
  for (int i=0; i<3; i++) printf("\t\t%.12lf", deltaS[1][i]);
  printf("\n\n");
  for (int i=0; i<3; i++) printf("\t\t%.12lf", deltaS[2][i]);
  printf("\n\n");

  return 0;
}

/*
Compile:
    clang voigt-enzyme.c -Xclang -load -Xclang /home/linuxbrew/.linuxbrew/Cellar/enzyme/HEAD-6e45ead/lib/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops

Output:

Swork       = 	   0.098041	   0.092640	   0.104300	   0.002458	   -0.000928	   0.009383

deltaSwork  = 	   0.169671	   0.219554	   0.099307	   -0.010067	   0.003474	   -0.085117

deltaS      =

		0.169670848078		-0.085117173016		0.003473687038

		-0.085117173016		0.219553855108		-0.010067073475

		0.003473687038		-0.010067073475		0.099307342055

*/
