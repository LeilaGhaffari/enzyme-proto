// Standard Neo-Hookean formulation in initial configuration with Enzyme-AD
// S = dPhi/dE

#include <stdio.h>
#include <math.h>
#include <string.h>

#define VECSIZE 6

double __enzyme_fwddiff(void*, ...);

int enzyme_dup;
int enzyme_const;

void computePhi(double *Phi, double E2work[VECSIZE], double mu);

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

void S_analytical(double S_an[VECSIZE], double E2work[VECSIZE], double mu) {
  
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
      S_an[m] = 0.;
      for ( int n = 0; n < 3; n++)
          S_an[m] += mu*C_inv[indj[m]][n]*E2[n][indk[m]];
  }
};

void computePhi(double *Phi, double E2work[VECSIZE], double mu) {
   
  // log(J)
  double detCm1 = computeDetCM1(E2work);
  double logJ = log1p_series_shifted(detCm1) / 2.;

  // trace(E)
  double traceE = (E2work[0] + E2work[1] + E2work[2]) / 2.;

  *Phi = mu * (-logJ + traceE);
};

// Kelvin-Mandel notation
void KMStrainRate(const double grad_u[3][3], double strain_rate[6]) {
  const double weight = 1 / sqrt(2.);
  strain_rate[0] = grad_u[0][0];
  strain_rate[1] = grad_u[1][1];
  strain_rate[2] = grad_u[2][2];
  strain_rate[3] = weight * (grad_u[2][1] + grad_u[1][2]);
  strain_rate[4] = weight * (grad_u[2][0] + grad_u[0][2]);
  strain_rate[5] = weight * (grad_u[1][0] + grad_u[0][1]);
}

void KMUnpack(const double v[6], double A[3][3]) {
  const double weight = 1 / sqrt(2.);
  A[0][0] = v[0];
  A[1][1] = v[1];
  A[2][2] = v[2];
  A[2][1] = A[1][2] = weight * v[3];
  A[2][0] = A[0][2] = weight * v[4];
  A[1][0] = A[0][1] = weight * v[5];
}

int main() {
  double mu = 1.;

  double grad_u[3][3] = {
    {
      0.5895232828911128,
      0.2362491738162759,
      0.9793730522395296
    },
    {
      0.2190993957421843,
      0.0126503210747925,
      0.6570956167695403
    },
    {
      0.4984241659585802,
      0.6163324783748984,
      0.9352322238627951
    }
  };

  double strain_rate[VECSIZE];
  KMStrainRate(grad_u, strain_rate);

  // Compute S with Enzyme-AD Forward mode
  double Phi;
  double S_ad[VECSIZE];
  for (int i=0; i<VECSIZE; i++) { 
    double dE[VECSIZE] = {0.}; dE[i] = 1.;
    __enzyme_fwddiff((void *)computePhi, 
                     &Phi, &S_ad[i], 
                     strain_rate, dE,
                     enzyme_const, mu);
  }
  // The first 3 entries (diagonal in Voigt notation) worked with 2E here, but
  // we need the gradient with respect to E. The next three (off-diagonal
  // components) represent entries that appear twice in the matrix, thus they
  // already have the necessary scaling.
  //for (int i=0; i<3; i++)
  //  S_ad[i] *= sqrt(2.);

  // Compute analytical S
  double S_an[VECSIZE];
  S_analytical(S_an, strain_rate, mu);

  printf("\n\nPhi       = ");
  printf("\t   %.6lf", Phi);

  printf("\n\nS_autodiff      =\n\n");
  for (int i=0; i<VECSIZE; i++) printf("\t\t%.12lf", S_ad[i]);
  printf("\n\n");

  printf("\n\nS_analytical    =\n\n");
  for (int i=0; i<VECSIZE; i++) printf("\t\t%.12lf", S_an[i]);
  printf("\n\n");

  return 0;
}

/*

Phi       =        0.359631

S_autodiff      =

                0.190093858759          -0.062130869557         0.482363568176          0.114768327105          -0.043324623403         0.438162617768



S_analytical    =

                0.190092241607          -0.062132990331         0.482362534603          0.114768556264          -0.043324709910         0.438163492654

*/
