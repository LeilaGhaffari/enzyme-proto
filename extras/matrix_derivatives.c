// Derivative af a matrix wrt a scalar, a vector, and a matrix
#include <stdio.h>
#include <stdlib.h>

int __enzyme_autodiff(void *, ...);
int enzyme_dup;

// Function that prints a Matrix
void print_matrix(int row, int col, double M[row][col]) {
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      printf("%f ", M[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

// 2x2 Matrix with a 2x2 Matrix input - done as 2 4x1 vectors
void mat_mat(double *M, double *x) {
  M[0] = x[2] * x[3] + x[0] * x[0];
  M[1] = 4 * x[1] + 15;
  M[2] = x[3] * x[3] * x[3] + 2 * x[0] * x[2];
  M[3] = 3 * x[0] + 4 * x[3];
}

// Goes from gradient matrix to tensor
void mat2tensor(int row1, int col1, int row2, int col2, double Mat[row1 + row2][col1 + col2], double tensor[row1][col1][row2][col2]) {
  int m = 0;
  // Outer row
  for (int i = 0; i < row1; i++) {
    // outer column
    for (int j = 0; j < col1; j++) {
      int n = 0;
      // inner row
      for (int k = 0; k < row2; k++) {
        // inner column
        for (int l = 0; l < col2; l++) {
          tensor[i][j][k][l] = Mat[m][n];
          n++;
        }
      }
      m++;
    }
  }
}

void get_grad(double *x, double grad[4][4]) {
  // Pre define
  double F[4] = {0.};

  // Calculate Jacobian
  for (int i = 0; i < 4; i++) {
    double dF[4] = {0.};
    dF[i]        = 1.;
    __enzyme_autodiff(mat_mat, enzyme_dup, F, dF, enzyme_dup, x, grad[i]);
  }
}

int main() {
  // Matrix wrt matrix - m x n x k x l tensor
  double z[4]     = {3., 5., 1., 7.};
  double dz[4][4] = {0.};

  get_grad(z, dz);

  print_matrix(4, 4, dz);

  double tens[2][2][2][2] = {0};
  mat2tensor(2, 2, 2, 2, dz, tens);

  print_matrix(2, 2, tens[0][0]);
  print_matrix(2, 2, tens[0][1]);
  print_matrix(2, 2, tens[1][0]);
  print_matrix(2, 2, tens[1][1]);
}