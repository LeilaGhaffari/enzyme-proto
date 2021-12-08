// Compute div(U)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

void __enzyme_autodiff(void *, ...);
int enzyme_const;

int ExactSolution(double q[], double *time, double *X, double *Y, double *Z) {

  // -- Time
  double t = time[0];
  // -- Coordinates
  double x = X[0];
  double y = Y[0];
  double z = Z[0];

  // Exact solutions
  q[0] = 2*(t*t) + (1*x*x) + (6 *y*y) +  (20*z*z);
  q[1] = 3*(t*t) + (2*x*x) + (7 *y*y) +  (30*z*z);
  q[2] = 4*(t*t) + (3*x*x) + (8 *y*y) +  (40*z*z);
  q[3] = 5*(t*t) + (4*x*x) + (9 *y*y) +  (50*z*z);
  q[4] = 6*(t*t) + (5*x*x) + (10*y*y) +  (60*z*z);

  return 0;
}

void Grad_q(double grad_q[3][5], double *x, double *y, double *z, double *t) {
    double q[5];
    // Derivative wrt x
    for (int i=0; i<5; i++) {
        double q_[5] = {0.}; q_[i] = 1.;
      __enzyme_autodiff((void *)ExactSolution,
                    q, q_,
                    enzyme_const, t,
                    x, &grad_q[0][i],
                    enzyme_const, y,
                    enzyme_const, z);
    }
    // Derivative wrt y
    for (int i=0; i<5; i++) {
        double q_[5] = {0.}; q_[i] = 1.;
      __enzyme_autodiff((void *)ExactSolution,
                    q, q_,
                    enzyme_const, t,
                    enzyme_const, x,
                    y, &grad_q[1][i],
                    enzyme_const, z);
    }
    // Derivative wrt z
    for (int i=0; i<5; i++) {
        double q_[5] = {0.}; q_[i] = 1.;
      __enzyme_autodiff((void *)ExactSolution,
                    q, q_,
                    enzyme_const, t,
                    enzyme_const, x,
                    enzyme_const, y,
                    z, &grad_q[2][i]);
    }
}

//int ConvectiveFlux(double *F_conv_x, double *F_conv_y, double *F_conv_z,  
//                   double *X, double *Y, double *Z, double *time) {
//  double q[5];
//  ExactSolution(q, time, X, Y, Z);
//  double rho = q[0];
//  double u[3] = {q[1]/rho, q[2]/rho, q[3]/rho};
//  double E = q[4];
//}

void computeQdot(double *q_dot, double *t, double *x, double *y, double *z) {
    double q[5];
    for (int i=0; i<5; i++) {
        double q_[5] = {0.}; q_[i] = 1.;
        __enzyme_autodiff((void *)ExactSolution,
                          q, q_,
                          t, &q_dot[i],
                          enzyme_const, x,
                          enzyme_const, y,
                          enzyme_const, z);
    }
}

int main() {
    // Declarations
    double X[] = {.5, 2.5, 5.};
    double x[1] = {X[0]}, y[1] = {X[1]}, z[1] = {X[2]};
    double time[1] = {.2};

    // -------------------------------------------------------------------------
    // q and q_dot
    // -------------------------------------------------------------------------
    double q_dot[5] = {0.};
    computeQdot(q_dot, time, x, y, z);

    // Print output
    printf("\nqdot\n");
    for (int j=0; j<5; j++) printf("%g\n", q_dot[j]);

    // -------------------------------------------------------------------------
    // div(U)
    // -------------------------------------------------------------------------
    // grad_q
    double grad_q[3][5]; // TODO: everything looks correct except for grad[0][4] (adds 1 to it)
    Grad_q(grad_q, x, y, z, time);

    // Print output
    printf("\ngrad_q\n");
    for (int i=0; i<3; i++) {
        printf("\nDerivative in the %d direction:\n", i);
      for (int j=0; j<5; j++) printf("%f\t", grad_q[i][j]);
      printf("\n\n");
    }

    return 0;
}

/*

clang test019.c -Xclang -load -Xclang /home/leila/Enzyme/enzyme/build12DHB/Enzyme/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops

*/
