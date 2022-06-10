#include <stdio.h>
#include <stdlib.h>

void __enzyme_augmentfwd(void *, ...);
void __enzyme_reverse(void *, ...);
int  __enzyme_augmentsize(void *, ...);

int  enzyme_dup;
int  enzyme_tape;
int  enzyme_allocated;

void func(double *V, double *x) {
  V[0] = x[0] * x[0] + 3 * x[1] + 5;
  V[1] = x[0] * x[1] + 3 * x[2] * x[1];
  V[2] = x[2] * x[2] * x[2];
}

void get_grad(double *V, double *dv, double *x, double *dx) {
  int   size = __enzyme_augmentsize((void *)func, enzyme_dup, enzyme_dup);
  void *data = malloc(size);
  __enzyme_augmentfwd((void *)func, enzyme_allocated, size, enzyme_tape, data, V, dv, x, dx);
  __enzyme_reverse((void *)func, enzyme_allocated, size, enzyme_tape, data, V, dv, x, dx);
  free(data);
}

// Function that prints a Matrix
void print_matrix(int row, int col, double M[row][col]) {
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      printf("%f ", M[i][j]);
    }
    printf("\n");
  }
}

int main() {
  double x[3]     = {2, 3, 5};
  double V[3]     = {0};
  double dx[3][3] = {0};

  for (int i = 0; i < 3; i++) {
    double dV[3] = {0.};
    dV[i]        = 1.;
    get_grad(V, dV, x, dx[i]);
  }

  print_matrix(3, 3, dx);
  return 0;
}