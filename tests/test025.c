// Compute d2q/dx2 - x[3] passed in

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

void __enzyme_autodiff(void *, ...);
int enzyme_const;

int ExactSolution(double q[], double *time, double X[]) {
  // -- Time
  double t = time[0];
  // -- Coordinates
  double x = X[0];
  double y = X[1];
  double z = X[2];

  // Exact solutions
  q[0] = 1.*x*x+x + 6. *y*y+y + 11.*z*z+z;
  q[1] = 2.*x*x+x + 7. *y*y+y + 12.*z*z+z;
  q[2] = 3.*x*x+x + 8. *y*y+y + 13.*z*z+z;
  q[3] = 4.*x*x+x + 9. *y*y+y + 14.*z*z+z;
  q[4] = 5.*x*x+x + 10.*y*y+y + 15.*z*z+z;

  return 0;
}

// -- grad_q
void computeGrad_q(double grad_q[5][3], double *t, double x[]) {
  double q[5];
  for (int i=0; i<5; i++) {
    double q_[5] = {0.}; q_[i]= 1.;
    __enzyme_autodiff((void *)ExactSolution,
                      q, q_,
                      enzyme_const, t,
                      x, grad_q[i]);          
  }
}

void compute_dq_x(double dq[5], double *t, double x[]) {
  double grad_q[5][3];
  computeGrad_q(grad_q, t, x);
  for (int i=0; i<5; i++) dq[i] = grad_q[i][0];
}

void compute_d2q_x(double dq[5], double *t, double x[]) {
  double dq_x[5];
  for (int i=0; i<5; i++) {
    double q_[5] = {0.}; q_[i]= 1.;
    double dq_[3] = {0.};
    __enzyme_autodiff((void *)compute_dq_x,
                      dq_x, q_,
                      enzyme_const, t,
                      x, dq_);
    dq[i] = dq_[0];
  }
}

// ***************************************************************************
// Main Function
// ***************************************************************************
int main() {
  // Declarations
  double X[] = {.5, 2.5, 5.};
  //double x[1] = {X[0]}, y[1] = {X[1]}, z[1] = {X[2]};
  double time[1] = {.2};
// -------------------------------------------------------------------------
  // grad_q
  // -------------------------------------------------------------------------
  double grad_q[5][3] = {{0.}};
  computeGrad_q(grad_q, time, X);
  for (int i=0; i<3; i++) {
    printf("\nDerivative in direction %d:\n", i);
    for (int j=0; j<5; j++) printf("%f\t", grad_q[j][i]);
    printf("\n");
  }
  // -------------------------------------------------------------------------
  // d2q/dx2
  // -------------------------------------------------------------------------
  double d2q_dx2[5] = {0.}; // Must initialize
  compute_d2q_x(d2q_dx2, time, X);
  // Print output
  printf("\nd2q/dx2:\n");
  for (int j=0; j<5; j++) printf("%f\t", d2q_dx2[j]);

  printf("\n");
  return 0;
}
// ***************************************************************************
/*
clang test025.c -Xclang -load -Xclang /home/leila/Enzyme/enzyme/build12DHB/Enzyme/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops

Derivative in direction 0:
2.000000	3.000000	4.000000	5.000000	6.000000	

Derivative in direction 1:
31.000000	36.000000	41.000000	46.000000	51.000000	

Derivative in direction 2:
111.000000	121.000000	131.000000	141.000000	151.000000	

d2q/dx2:
2.000000	4.000000	6.000000	8.000000	10.000000

*/