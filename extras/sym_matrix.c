// Takes the derivative of  a symetric matrix wrt another symetric matirx

#include <stdio.h>
#include <stdlib.h>

#define VECSIZE 6

int  __enzyme_autodiff(void *, ...);
int  enzyme_dup;

void func(double *V, double *x) {
  V[0] = x[2] * x[1];
  V[1] = x[1] * x[1] + 3 * x[5];
  V[2] = x[0] + 4 * x[3];
  V[3] = x[4] * x[4] + 5 * x[5];
  V[4] = x[1] * x[4] * x[5] + 2 * x[0];
  V[5] = x[2];
}

int main() {
  double x[6] = {2., 6., 11., 5., 7., 4.};
  double V[6] = {0., 0., 0., 0., 0., 0.};

  double dx[6][6] = {0};

  for (int i = 0; i < 6; i++) {
    double dV[6] = {0., 0., 0., 0., 0., 0.};
    dV[i]        = 1;
    __enzyme_autodiff((void *)func, enzyme_dup, V, dV, enzyme_dup, x, dx[i]);
  }

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      printf("%f ", dx[i][j]);
    }
    printf("\n");
  }

  return 0;
}