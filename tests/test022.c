
/*const double Fu[6]     =  {mu*(dudx[0][0] * (2 + lambda) + 
                               lambda * (dudx[1][1] + dudx[2][2])),
                                   mu*(dudx[0][1] + dudx[1][0]), 
                                   mu*(dudx[0][2] + dudx[2][0]), 
                                   mu*(dudx[1][1] * (2 + lambda) + 
                                       lambda * (dudx[0][0] + dudx[2][2])),
                                   mu*(dudx[1][2] + dudx[2][1]), 
                                   mu*(dudx[2][2] * (2 + lambda) + 
                                       lambda * (dudx[0][0] + dudx[1][1]))
                                  };
const CeedInt Fuviscidx[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}}; // symmetric matrix indices
const CeedScalar Fe[3]     =  {u[0]*Fu[0] + u[1]*Fu[1] + u[2]*Fu[2] + 
                                   k*grad_T[0], 
                                   u[0]*Fu[1] + u[1]*Fu[3] + u[2]*Fu[4] + 
                                   k*grad_T[1], 
                                   u[0]*Fu[2] + u[1]*Fu[4] + u[2]*Fu[5] + 
                                   k*grad_T[2] 
                                  };
*/                                 

// Use autodiff for computing gradients
//   MMS in libCEED/examples/fluids follows this test

// Include C libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// ***************************************************************************
// Exact Solution
// ***************************************************************************
int ExactSolution(double q[], double *time, double *X, double *Y, double *Z) {

  // -- Time
  double t = time[0];
  // -- Coordinates
  double x = X[0];
  double y = Y[0];
  double z = Z[0];

  // Exact solutions
  q[0] = 10.*(t*t+t+1.) + (1.*x) + (6. *y) + (11.*z);
  q[1] = 20.*(t*t+t+1.) + (2.*x) + (7. *y) + (12.*z);
  q[2] = 30.*(t*t+t+1.) + (3.*x) + (8. *y) + (13.*z);
  q[3] = 40.*(t*t+t+1.) + (4.*x) + (9. *y) + (14.*z);
  q[4] = 50.*(t*t+t+1.) + (5.*x) + (10.*y) + (15.*z);

  return 0;
}
// ***************************************************************************

// ---------------------------------------------------------------------------
// Enzyme AD
// ---------------------------------------------------------------------------
// -- Enzyme functions and variables
void __enzyme_autodiff(void *, ...);
int enzyme_const;

// -- Q_dot
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

void computeF(double f[3][5], double *t, double *x, double *y, double *z, 
              double lambda, double mu, double k, double cv, double cp, double g) {
  // Compute state variables                
  double q[5];
  ExactSolution(q, t, x, y, z);
  double rho = q[0];
  double u[3] = {q[0]/rho, q[1]/rho, q[2]/rho};
  double E = q[4];

  // Density equation - Advective flux
  double F_adv_density[3] = {rho*u[0], rho*u[1], rho*u[2]};

  // Momentum equations 
  // -- Advective flux
  double gamma = cp/cv;
  double kinetic_energy = (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]) / 2.;
  double P = (E - kinetic_energy * rho - rho*g*z[0]) * (gamma - 1.);
  double F_adv_momentum = {{rho*u[0]*u[0] + P, rho*u[0]*u[1], rho*u[0]*u[2]}
                           {}, 
                           {}};


  


}

void computeGrad_F(double grad_f[3][5], double *t, double *x, double *y, double *z,
                   double lambda, double mu, double k, double cv, double cp, double g) { 
  double q[5];
  // Derivative wrt x
  for (int i=0; i<5; i++) {
    double q_[5] = {0.}; q_[i] = 1.;
    __enzyme_autodiff((void *)computeF,
                      q, q_,
                      enzyme_const, t,
                      x, &grad_f[0][i],
                      enzyme_const, y,
                      enzyme_const, z,
                      enzyme_const, lambda,
                      enzyme_const, mu,
                      enzyme_const, k,
                      enzyme_const, cv,
                      enzyme_const, cp,
                      enzyme_const, g);
  }
  // Derivative wrt y
  for (int i=0; i<5; i++) {
    double q_[5] = {0.}; q_[i] = 1.;
    __enzyme_autodiff((void *)computeF,
                      q, q_,
                      enzyme_const, t,
                      enzyme_const, x,
                      y, &grad_f[1][i],
                      enzyme_const, z,
                      enzyme_const, lambda,
                      enzyme_const, mu,
                      enzyme_const, k,
                      enzyme_const, cv,
                      enzyme_const, cp,
                      enzyme_const, g);                      
  }
  // Derivative wrt z
  for (int i=0; i<5; i++) {
    double q_[5] = {0.}; q_[i] = 1.;
    __enzyme_autodiff((void *)computeF,
                      q, q_,
                      enzyme_const, t,
                      enzyme_const, x,
                      enzyme_const, y,
                      z, &grad_f[2][i],
                      enzyme_const, lambda,
                      enzyme_const, mu,
                      enzyme_const, k,
                      enzyme_const, cv,
                      enzyme_const, cp,
                      enzyme_const, g);
  }
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
  // q_dot
  // -------------------------------------------------------------------------
  double q_dot[5] = {0.};       // Must initialize
  computeQdot(q_dot, time, x, y, z);
  // Print output
  printf("\nq_dot:\n");
  for (int j=0; j<5; j++) printf("%f\t", q_dot[j]);
  printf("\n");
  // -------------------------------------------------------------------------
  // grad_q
  // -------------------------------------------------------------------------
  double grad_q[3][5] = {{0.}}; // Must initialize
  computeGrad_q(grad_q, time, x, y, z);
  // Print output
  for (int i=0; i<3; i++) {
    printf("\nDerivative in direction %d:\n", i);
    for (int j=0; j<5; j++) printf("%f\t", grad_q[i][j]);
    printf("\n");
  }
  // -------------------------------------------------------------------------
  // Flux
  // -------------------------------------------------------------------------
  // -- Physical properties
  double lambda = -2./3.;
  double mu     = 75.;
  double k      = 0.02638;
  double cv     = 717.;
  double cp     = 1004.;
  double g      = 9.81;

  // -- Primary Units
  double meter    = 1e-2;  // 1 meter in scaled length units
  double kilogram = 1e-6;  // 1 kilogram in scaled mass units
  double second   = 1e-2;  // 1 second in scaled time units
  double Kelvin   = 1;     // 1 Kelvin in scaled temperature units

  // -- Secondary Units
  double W_per_m_K, Pascal, J_per_kg_K, m_per_squared_s;
  Pascal          = kilogram / (meter * PetscSqr(second));
  J_per_kg_K      =  PetscSqr(meter) / (PetscSqr(second) * Kelvin);
  m_per_squared_s = meter / PetscSqr(second);
  W_per_m_K       = kilogram * meter / (pow(second,3) * Kelvin);

  // -- Unit conversion
  cv *= J_per_kg_K;
  cp *= J_per_kg_K;
  g  *= m_per_squared_s;
  mu *= Pascal * second;
  k  *= W_per_m_K;
  double gamma  = cp / cv;
  
  // -- Compute div(Flux)
  double grad_F[3][5] = {{0.}}; // Must initialize
  computeGrad_F(grad_F, time, x, y, z, lambda, mu, k, cv, cp, g);
  // Print output
  for (int i=0; i<3; i++) {
    printf("\nDerivative in direction %d:\n", i);
    for (int j=0; j<5; j++) printf("%f\t", grad_F[i][j]);
    printf("\n");
  }


  return 0;
}
// ***************************************************************************

/*
clang test022.c -Xclang -load -Xclang /home/leila/Enzyme/enzyme/build12DHB/Enzyme/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops

*/
