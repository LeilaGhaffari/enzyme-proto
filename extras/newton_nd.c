#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int enzyme_const;
int enzyme_dup;

// Define a function pointer?
typedef void(f_ptr)(double *, double *);

extern void __enzyme_autodiff(f_ptr, int, double *, double *, int, double *, double *);

void        func(double *F, double *x) {
  // Define equation and size
  F[0] = 3 * x[0] * x[0] - x[1] * x[1];
  F[1] = 3 * x[0] * x[1] - x[0] * x[0] * x[0] - 1;
}

void get_Jacobian(double *x, double J[2][2]) {
  // Pre define
  double F[2] = {0., 0.};

  // Calculate Jacobian
  for (int i = 0; i < 2; i++) {
    double dF[2] = {0., 0.};
    dF[i]        = 1;
    __enzyme_autodiff(func, enzyme_dup, F, dF, enzyme_dup, x, J[i]);
  }
}

// Does Gauss Elim to solve Ax = b for a 2x2 A matrix
void gauss(double A[2][2], double x[2], double b[2]) {
  double ratio = A[1][0] / A[0][0];
  A[1][0]      = 0;
  A[1][1]      = A[1][1] - ratio * A[0][1];
  b[1]         = b[1] - ratio * b[0];
  x[1]         = b[1] / A[1][1];
  x[0]         = (b[0] - A[0][1] * x[1]) / A[0][1];
}

// Function for newtons method in 2d
void newton_nd(int n_max, double target_tol, double guess[2], double xn[2]) {
  double tol = 1;
  int    n   = 0;

  printf("|  n  |  x  |  f(x)  |\n|  %i  |  {%f,%f}  |    |\n", n, xn[0], xn[0]);
  while (tol > target_tol && n < n_max) {
    // Solve J*diff = F for diff
    double A[2][2] = {0};
    double diff[2] = {0};
    double F[2];
    func(F, xn);
    F[0] = -F[0];
    F[1] = -F[1];
    get_Jacobian(xn, A);
    gauss(A, diff, F);
    for (int i = 0; i < 2; i++) {
      xn[i] = diff[i] + xn[i];
    }
    tol = sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
    n++;
    printf("|  n  |  x  |  f(x)  |\n|  %i  |  {%f,%f}  |  {%f,%f}  |\n", n, xn[0], xn[0], -F[0], -F[1]);
  }
}

int main() {
  double guess[2] = {2, 2};
  double xn[2]    = {0};
  newton_nd(100, 1e10, guess, xn);

  printf("The solution is [%f,%f]\n", xn[0], xn[1]);
  return 0;
}