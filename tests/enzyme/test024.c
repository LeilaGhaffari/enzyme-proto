// Compute d2q/dx2

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
  q[0] = 1.*x*x+x + 6. *y*y+y + 11.*z*z+z;
  q[1] = 2.*x*x+x + 7. *y*y+y + 12.*z*z+z;
  q[2] = 3.*x*x+x + 8. *y*y+y + 13.*z*z+z;
  q[3] = 4.*x*x+x + 9. *y*y+y + 14.*z*z+z;
  q[4] = 5.*x*x+x + 10.*y*y+y + 15.*z*z+z;

  return 0;
}

// -- grad_q
void computeGrad_q(double grad_q[3][5], double *t, double *x, double *y, double *z) {
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

void compute_dq_x(double dq[5], double *t, double *x, double *y, double *z) {
  double grad_q[3][5];
  computeGrad_q(grad_q, t, x, y, z);
  for (int i=0; i<5; i++) dq[i] = grad_q[0][i];
}

void compute_d2q_x(double dq[5], double *t, double *x, double *y, double *z) {
  double dq_x[5];
  // Derivative wrt x
  for (int i=0; i<5; i++) {
    double q_[5] = {0.}; q_[i] = 1.;
    __enzyme_autodiff((void *)compute_dq_x,
                      dq_x, q_,
                      enzyme_const, t,
                      x, &dq[i],
                      enzyme_const, y,
                      enzyme_const, z);
  }
}

 // -- Temperature: T =  (E / rho - (u u)/2 - g z) / cv
void computeT(double *T, double *t, double *x, double *y, double *z, double g, double cv) {
  // Compute state variables                
  double q[5];
  ExactSolution(q, t, x, y, z);
  double rho = q[0];
  double u[3] = {q[1]/rho, q[2]/rho, q[3]/rho};
  double E = q[4];

  // Kinetic Energy
  double kinetic_energy = (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]) / 2.;

  // Compute T
  T[0] =  (E/rho - kinetic_energy - g*z[0]) / cv;

}

// -- Grad_T
void computeGrad_T(double dT[3], double *t, double *x, double *y, double *z, double g, double cv) {
  double T[0], T_ = 1.;
  // Derivative wrt x
  __enzyme_autodiff((void *)computeT,
                    T, &T_,
                    enzyme_const, t,
                    x, &dT[0],
                    enzyme_const, y,
                    enzyme_const, z,
                    enzyme_const, g,
                    enzyme_const, cv);
  // Derivative wrt y
  __enzyme_autodiff((void *)computeT,
                    T, &T_,
                    enzyme_const, t,
                    enzyme_const, x,
                    y, &dT[1],
                    enzyme_const, z,
                    enzyme_const, g,
                    enzyme_const, cv);
  // Derivative wrt z
  __enzyme_autodiff((void *)computeT,
                    T, &T_,
                    enzyme_const, t,
                    enzyme_const, x,
                    enzyme_const, y,
                    z, &dT[2],
                    enzyme_const, g,
                    enzyme_const, cv);
}

// ***************************************************************************
// Main Function
// ***************************************************************************
int main() {
  // Declarations
  double X[] = {.5, 2.5, 5.};
  double x[1] = {X[0]}, y[1] = {X[1]}, z[1] = {X[2]};
  double time[1] = {.2};
  // -------------------------------------------------------------------------
  // grad_q
  // -------------------------------------------------------------------------
  double grad_q[3][5] = {{0.}};
  computeGrad_q(grad_q, time, x, y, z);
  for (int i=0; i<3; i++) {
    printf("\nDerivative in direction %d:\n", i);
    for (int j=0; j<5; j++) printf("%f\t", grad_q[i][j]);
    printf("\n");
  }
  // -------------------------------------------------------------------------
  // d2q/dx2
  // -------------------------------------------------------------------------
  double d2q_dx2[5] = {0.}; // Must initialize
  compute_d2q_x(d2q_dx2, time, x, y, z);
  // Print output
  printf("\nd2q/dx2:\n");
  for (int j=0; j<5; j++) printf("%f\t", d2q_dx2[j]);

  // -------------------------------------------------------------------------
  // dT
  // -------------------------------------------------------------------------
  double g = 58., cv = 2649.;
  double dT[3] = {0.};
  computeGrad_T(dT, time, x, y, z, g, cv);
  // Print output
  printf("\n\ndT:\n");
  for (int j=0; j<3; j++) printf("%f\t", dT[j]);

  printf("\n");
  return 0;
}
// ***************************************************************************

/*

Derivative in direction 0:
2.000000	3.000000	4.000000	5.000000	6.000000	

Derivative in direction 1:
31.000000	36.000000	41.000000	46.000000	51.000000	

Derivative in direction 2:
111.000000	121.000000	131.000000	141.000000	151.000000	

d2q/dx2:
2.000000	4.000000	6.000000	8.000000	10.000000	

dT:
-0.000003	0.000000	0.000000

*/