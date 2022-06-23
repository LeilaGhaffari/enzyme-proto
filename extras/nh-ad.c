#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define VECSIZE 6

int __enzyme_autodiff(void *, ...);
int enzyme_dup; 
int enzyme_const;
// Computes J - volume change
double get_Jm1(double E[VECSIZE]){
    // J = |F| = sqrt(|C|) = sqrt(|I_3 + 2*E|)
    // Calculate C = [a, d, f, e, c, b]
    double C[6] = {2.*E[0]+1., 2.*E[1]+1., 2.*E[2]+1., 2.*E[3], 2.*E[4], 2.*E[5]};

    // Calculate determinate of C = adf - (ae^2 + dc^2 + fb^2) + 2bce
    double det_C = C[0]*C[1]*C[3] - (C[0]*C[3]*C[3] + C[1]*C[4]*C[4] + C[2]*C[5]*C[5]) + 2.*C[5]*C[4]*C[3];

    // Calculate and return J
    double Jm1 = sqrtf(fabs(det_C)) - 1.;
    return Jm1;
    
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

/* Function that defines the Strain energy. Currently unsing the neo-hookean model, but would like
user input for this function*/
void psi_func(double *psi, double E[VECSIZE], double lambda, double mu){

/* 
Inputs
---------------------------------------------------------------------
E      - 
lambda - Lame's first parameter 
mu     - Shear modulus
J      - Determintate of Deformation (Volume Change)

Output
---------------------------------------------------------------------
psi - Elastic Stain Energy
*/

double Jm1 = get_Jm1(E);

// psi = lambda/2*log(J)^2 - mu*log(J) + mu*trace(E)
*psi = lambda/2.*log1p_series_shifted(Jm1)*log1p_series_shifted(Jm1) - mu*log1p_series_shifted(Jm1) + mu*(E[0] + E[1] + E[2]); 
}

// Computes S via eq 15(St. Venant-Kirchoff model) to check S computed by enzyme
void test_compute_S(double Swork[VECSIZE], double E2work[VECSIZE], double lambda, double mu){
/*double trace_E = E[0] + E[1] + E[2];

for(int i = 0; i<6; i++){
    S[i] = 2.*mu*E[i];
    if(i<3){
        S[i] += lambda*trace_E;
    }
}*/

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
}


// Function that prints 3x3 a symetric Matrix from vector
void print_symmatrix( double vec[6]) {
  double M[3][3] = {{vec[0], vec[5], vec[4]},
                    {vec[5], vec[1], vec[3]},
                    {vec[4], vec[3], vec[2]}};
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      printf("%f ", M[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

// Computes the second Piola-Kirchoff tensor using enzyme
void compute_S(double S[VECSIZE], double E[VECSIZE], double lambda, double mu){

    // Enzyme call
    double psi;
    for(int i=0; i<VECSIZE; i++){
    double dE[VECSIZE] = {0.}; dE[i] = 1.;
    __enzyme_autodiff((void *)psi_func,
                        enzyme_dup, &psi, &S[i],
                        enzyme_dup, E, dE,
                        enzyme_const, lambda,
                        enzyme_const, mu);

    }
}

int main(){
    // Define E
    double E[VECSIZE] = {.2, .5, .4, .8, .6, 1};

    // Define constants 
    double lambda = 1.; 
    double mu = 1.; 
    double Jm1 = get_Jm1(E);

    // Compute S's 
    double S_ad[VECSIZE];
    double S_an[VECSIZE];

    compute_S(S_ad, E, lambda, mu);
    test_compute_S(S_an, E, lambda, mu);

    print_symmatrix(S_ad);
    print_symmatrix(S_an);

}