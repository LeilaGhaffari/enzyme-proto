                             

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
  
  // Compute gradient of state variables
  //double dq[3][5] = {{1.,	2.,	3.,	4.,	5.,},
  //                   {6.,	7.,	8.,	9.,	10.},
  //                   {11., 12.,	13.,	14., 15.}
  //                  };
  double dq[3][5] = {{0.}};
  computeGrad_q(dq, t, x, y, z);

  double drho[3] =   {dq[0][0],
                      dq[1][0],
                      dq[2][0]
                     };
  double dU[3][3] = {{dq[0][1],
                      dq[1][1],
                      dq[2][1]},
                     {dq[0][2],
                      dq[1][2],
                      dq[2][2]},
                     {dq[0][3],
                      dq[1][3],
                      dq[2][3]}
                    };
  double dE[3] =     {dq[0][4],
                      dq[1][4],
                      dq[2][4]
                     };

  double du[3][3] = {{0.}};
  for (int j=0; j<3; j++)
    for (int k=0; k<3; k++)
      du[j][k] = (dU[j][k] - drho[k]*u[j]) / rho;                   

  // Density 
  // -- Advective flux
  //for (int j=0; j<3; j++) f[j][0] += rho*u[j];
  //for (int j=0; j<3; j++) {
  //  printf("\n");
  //  for (int k=0; k<5; k++) printf("force[%d][%d] = ")
  //}
  double F_adv_density[3] = {rho*u[0], rho*u[1], rho*u[2]};

  // -- No diffusive flux

  // Momentum 
  // -- Advective flux
  double gamma = cp/cv;
  double kinetic_energy = (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]) / 2.;
  double P = (E - kinetic_energy * rho - rho*g*z[0]) * (gamma - 1.);
  double F_adv_momentum[3][3] = {{rho*u[0]*u[0] + P, rho*u[0]*u[1],     rho*u[0]*u[2]},
                                 {rho*u[1]*u[0],     rho*u[1]*u[1] + P, rho*u[1]*u[2]}, 
                                 {rho*u[2]*u[0],     rho*u[2]*u[1],     rho*u[2]*u[02] + P}};
  // -- Diffusive Flux
  double Fu[6] = {mu*(du[0][0] * (2 + lambda) + lambda * (du[1][1] + du[2][2])),
                  mu*(du[0][1] + du[1][0]), 
                  mu*(du[0][2] + du[2][0]), 
                  mu*(du[1][1] * (2 + lambda) + lambda * (du[0][0] + du[2][2])),
                  mu*(du[1][2] + du[2][1]), 
                  mu*(du[2][2] * (2 + lambda) + lambda * (du[0][0] + du[1][1]))
                 };

  const int Fuviscidx[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}}; // symmetric matrix indices
  double F_dif_momentum[3][3];
  for (int j=0; j<3; j++) 
    for (int k=0; k<3; k++) 
      F_dif_momentum[j][k] = Fu[Fuviscidx[j][k]];

  // Total Energy
  // -- Advective flux
  double F_adv_energy[3] = {(E + P)*u[0], (E + P)*u[1], (E + P)*u[2]};

  // -- Diffusive Flux 
  double dT[3] = {(dE[0]/rho - E*drho[0]/(rho*rho) - (u[0]*du[0][0] + u[1]*du[1][0] + u[2]*du[2][0]))    /cv,
                  (dE[1]/rho - E*drho[1]/(rho*rho) - (u[0]*du[0][1] + u[1]*du[1][1] + u[2]*du[2][1]))    /cv,
                  (dE[2]/rho - E*drho[2]/(rho*rho) - (u[0]*du[0][2] + u[1]*du[1][2] + u[2]*du[2][2]) - g)/cv
                 };
  double F_dif_energy[3] = {u[0]*Fu[0] + u[1]*Fu[1] + u[2]*Fu[2] + k*dT[0], 
                            u[0]*Fu[1] + u[1]*Fu[3] + u[2]*Fu[4] + k*dT[1], 
                            u[0]*Fu[2] + u[1]*Fu[4] + u[2]*Fu[5] + k*dT[2] 
                           };
  
  // Populate Flux
  // -- Zero f
  for (int j=0; j<3; j++)  for (int k=0; k<5; k++) f[j][k] = 0.;
  
  // -- Density
  for (int i=0; i<3; i++) f[i][0] += F_adv_density[i];

  // -- Momentum
  for (int i=0; i<3; i++) 
    for (int j=0; j<3; j++) {
      f[i][j+1] += F_adv_momentum[j][i];
      f[i][j+1] -= F_dif_momentum[j][i];
    }

  // -- Energy
  for (int i=0; i<3; i++) {
      f[i][4] += F_adv_energy[i];
      f[i][4] -= F_dif_energy[i];
  } 
}

// -- dFlux[0]/dx
void computeF0(double f0[5], double *t, double *x, double *y, double *z, 
               double lambda, double mu, double k, double cv, double cp, double g) {              
  // Compute state variables                
  double q[5];
  ExactSolution(q, t, x, y, z);
  double rho = q[0];
  double u[3] = {q[0]/rho, q[1]/rho, q[2]/rho};
  double E = q[4];
  
  // Compute gradient of state variables
  double dq[3][5] = {{0.}};
  computeGrad_q(dq, t, x, y, z);

  double drho[3] =   {dq[0][0],
                      dq[1][0],
                      dq[2][0]
                     };
  double dU[3][3] = {{dq[0][1],
                      dq[1][1],
                      dq[2][1]},
                     {dq[0][2],
                      dq[1][2],
                      dq[2][2]},
                     {dq[0][3],
                      dq[1][3],
                      dq[2][3]}
                    };
  double dE[3] =     {dq[0][4],
                      dq[1][4],
                      dq[2][4]
                     };

  double du[3][3] = {{0.}};
  for (int j=0; j<3; j++)
    for (int k=0; k<3; k++)
      du[j][k] = (dU[j][k] - drho[k]*u[j]) / rho;                   

  // Density 
  // -- Advective flux
  double F_adv_density[3] = {rho*u[0], rho*u[1], rho*u[2]};

  // -- No diffusive flux

  // Momentum 
  // -- Advective flux
  double gamma = cp/cv;
  double kinetic_energy = (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]) / 2.;
  double P = (E - kinetic_energy * rho - rho*g*z[0]) * (gamma - 1.);
  double F_adv_momentum[3][3] = {{rho*u[0]*u[0] + P, rho*u[0]*u[1],     rho*u[0]*u[2]},
                                 {rho*u[1]*u[0],     rho*u[1]*u[1] + P, rho*u[1]*u[2]}, 
                                 {rho*u[2]*u[0],     rho*u[2]*u[1],     rho*u[2]*u[02] + P}};
  // -- Diffusive Flux
  double Fu[6] = {mu*(du[0][0] * (2 + lambda) + lambda * (du[1][1] + du[2][2])),
                  mu*(du[0][1] + du[1][0]), 
                  mu*(du[0][2] + du[2][0]), 
                  mu*(du[1][1] * (2 + lambda) + lambda * (du[0][0] + du[2][2])),
                  mu*(du[1][2] + du[2][1]), 
                  mu*(du[2][2] * (2 + lambda) + lambda * (du[0][0] + du[1][1]))
                 };

  const int Fuviscidx[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}}; // symmetric matrix indices
  double F_dif_momentum[3][3];
  for (int j=0; j<3; j++) 
    for (int k=0; k<3; k++) 
      F_dif_momentum[j][k] = Fu[Fuviscidx[j][k]];

  // Total Energy
  // -- Advective flux
  double F_adv_energy[3] = {(E + P)*u[0], (E + P)*u[1], (E + P)*u[2]};

  // -- Diffusive Flux  
  double dT[3] = {(dE[0]/rho - E*drho[0]/(rho*rho) - (u[0]*du[0][0] + u[1]*du[1][0] + u[2]*du[2][0]))    /cv,
                  (dE[1]/rho - E*drho[1]/(rho*rho) - (u[0]*du[0][1] + u[1]*du[1][1] + u[2]*du[2][1]))    /cv,
                  (dE[2]/rho - E*drho[2]/(rho*rho) - (u[0]*du[0][2] + u[1]*du[1][2] + u[2]*du[2][2]) - g)/cv
                 };
  double F_dif_energy[3] = {u[0]*Fu[0] + u[1]*Fu[1] + u[2]*Fu[2] + k*dT[0], 
                            u[0]*Fu[1] + u[1]*Fu[3] + u[2]*Fu[4] + k*dT[1], 
                            u[0]*Fu[2] + u[1]*Fu[4] + u[2]*Fu[5] + k*dT[2] 
                           };
  
  // Populate Flux
  double f[3][5] = {{0.}};
  // -- Density
  for (int i=0; i<3; i++) f[i][0] += F_adv_density[i];

  // -- Momentum
  for (int i=0; i<3; i++) 
    for (int j=0; j<3; j++) {
      f[i][j+1] += F_adv_momentum[j][i];
      f[i][j+1] -= F_dif_momentum[j][i];
    }

  // -- Energy
  for (int i=0; i<3; i++) {
      f[i][4] += F_adv_energy[i];
      f[i][4] -= F_dif_energy[i];
  }

  for (int i=0; i<5; i++) f0[i] = f[0][i];

}

void compute_dF0_dx(double df0_dx[5], double *t, double *x, double *y, double *z,
                    double lambda, double mu, double k, double cv, double cp, double g) {
  double f0[5];
  for (int i=0; i<5; i++) {
    double f0_[5] = {0.}; f0_[i] = 1.;
    __enzyme_autodiff((void *)computeF0,
                      f0, f0_,
                      enzyme_const, t,
                      x, &df0_dx[i],
                      enzyme_const, y,
                      enzyme_const, z,
                      enzyme_const, lambda,
                      enzyme_const, mu,
                      enzyme_const, k,
                      enzyme_const, cv,
                      enzyme_const, cp,
                      enzyme_const, g);
  }
}

// -- dFlux[1]/dy
void computeF1(double f1[5], double *t, double *x, double *y, double *z, 
               double lambda, double mu, double k, double cv, double cp, double g) {
  // Compute state variables                
  double q[5];
  ExactSolution(q, t, x, y, z);
  double rho = q[0];
  double u[3] = {q[0]/rho, q[1]/rho, q[2]/rho};
  double E = q[4];
  
  // Compute gradient of state variables
  double dq[3][5] = {{0.}};
  computeGrad_q(dq, t, x, y, z);

  double drho[3] =   {dq[0][0],
                      dq[1][0],
                      dq[2][0]
                     };
  double dU[3][3] = {{dq[0][1],
                      dq[1][1],
                      dq[2][1]},
                     {dq[0][2],
                      dq[1][2],
                      dq[2][2]},
                     {dq[0][3],
                      dq[1][3],
                      dq[2][3]}
                    };
  double dE[3] =     {dq[0][4],
                      dq[1][4],
                      dq[2][4]
                     };

  double du[3][3] = {{0.}};
  for (int j=0; j<3; j++)
    for (int k=0; k<3; k++)
      du[j][k] = (dU[j][k] - drho[k]*u[j]) / rho;                   

  // Density 
  // -- Advective flux
  double F_adv_density[3] = {rho*u[0], rho*u[1], rho*u[2]};

  // -- No diffusive flux

  // Momentum 
  // -- Advective flux
  double gamma = cp/cv;
  double kinetic_energy = (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]) / 2.;
  double P = (E - kinetic_energy * rho - rho*g*z[0]) * (gamma - 1.);
  double F_adv_momentum[3][3] = {{rho*u[0]*u[0] + P, rho*u[0]*u[1],     rho*u[0]*u[2]},
                                 {rho*u[1]*u[0],     rho*u[1]*u[1] + P, rho*u[1]*u[2]}, 
                                 {rho*u[2]*u[0],     rho*u[2]*u[1],     rho*u[2]*u[02] + P}};
  // -- Diffusive Flux
  double Fu[6] = {mu*(du[0][0] * (2 + lambda) + lambda * (du[1][1] + du[2][2])),
                  mu*(du[0][1] + du[1][0]), 
                  mu*(du[0][2] + du[2][0]), 
                  mu*(du[1][1] * (2 + lambda) + lambda * (du[0][0] + du[2][2])),
                  mu*(du[1][2] + du[2][1]), 
                  mu*(du[2][2] * (2 + lambda) + lambda * (du[0][0] + du[1][1]))
                 };

  const int Fuviscidx[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}}; // symmetric matrix indices
  double F_dif_momentum[3][3];
  for (int j=0; j<3; j++) 
    for (int k=0; k<3; k++) 
      F_dif_momentum[j][k] = Fu[Fuviscidx[j][k]];

  // Total Energy
  // -- Advective flux
  double F_adv_energy[3] = {(E + P)*u[0], (E + P)*u[1], (E + P)*u[2]};

  // -- Diffusive Flux  
  double dT[3] = {(dE[0]/rho - E*drho[0]/(rho*rho) - (u[0]*du[0][0] + u[1]*du[1][0] + u[2]*du[2][0]))    /cv,
                  (dE[1]/rho - E*drho[1]/(rho*rho) - (u[0]*du[0][1] + u[1]*du[1][1] + u[2]*du[2][1]))    /cv,
                  (dE[2]/rho - E*drho[2]/(rho*rho) - (u[0]*du[0][2] + u[1]*du[1][2] + u[2]*du[2][2]) - g)/cv
                 };
  double F_dif_energy[3] = {u[0]*Fu[0] + u[1]*Fu[1] + u[2]*Fu[2] + k*dT[0], 
                            u[0]*Fu[1] + u[1]*Fu[3] + u[2]*Fu[4] + k*dT[1], 
                            u[0]*Fu[2] + u[1]*Fu[4] + u[2]*Fu[5] + k*dT[2] 
                           };
  
  // Populate Flux
  double f[3][5] = {{0.}};
  // -- Density
  for (int i=0; i<3; i++) f[i][0] += F_adv_density[i];

  // -- Momentum
  for (int i=0; i<3; i++) 
    for (int j=0; j<3; j++) {
      f[i][j+1] += F_adv_momentum[j][i];
      f[i][j+1] -= F_dif_momentum[j][i];
    }

  // -- Energy
  for (int i=0; i<3; i++) {
      f[i][4] += F_adv_energy[i];
      f[i][4] -= F_dif_energy[i];
  }

  for (int i=0; i<5; i++) f1[i] = f[1][i];
}
void compute_dF1_dy(double df1_dy[5], double *t, double *x, double *y, double *z,
                   double lambda, double mu, double k, double cv, double cp, double g) {
  double f1[5];
  for (int i=0; i<5; i++) {
    double f1_[5] = {0.}; f1_[i] = 1.;
    __enzyme_autodiff((void *)computeF1,
                      f1, f1_,
                      enzyme_const, t,
                      enzyme_const, x,
                      y, &df1_dy[i],
                      enzyme_const, z,
                      enzyme_const, lambda,
                      enzyme_const, mu,
                      enzyme_const, k,
                      enzyme_const, cv,
                      enzyme_const, cp,
                      enzyme_const, g);
  }
}

// -- dFlux[2]/dz
void computeF2(double f2[5], double *t, double *x, double *y, double *z, 
              double lambda, double mu, double k, double cv, double cp, double g) {
  // Compute state variables                
  double q[5];
  ExactSolution(q, t, x, y, z);
  double rho = q[0];
  double u[3] = {q[0]/rho, q[1]/rho, q[2]/rho};
  double E = q[4];
  
  // Compute gradient of state variables
  double dq[3][5] = {{0.}};
  computeGrad_q(dq, t, x, y, z);

  double drho[3] =   {dq[0][0],
                      dq[1][0],
                      dq[2][0]
                     };
  double dU[3][3] = {{dq[0][1],
                      dq[1][1],
                      dq[2][1]},
                     {dq[0][2],
                      dq[1][2],
                      dq[2][2]},
                     {dq[0][3],
                      dq[1][3],
                      dq[2][3]}
                    };
  double dE[3] =     {dq[0][4],
                      dq[1][4],
                      dq[2][4]
                     };

  double du[3][3] = {{0.}};
  for (int j=0; j<3; j++)
    for (int k=0; k<3; k++)
      du[j][k] = (dU[j][k] - drho[k]*u[j]) / rho;                   

  // Density 
  // -- Advective flux
  double F_adv_density[3] = {rho*u[0], rho*u[1], rho*u[2]};

  // -- No diffusive flux

  // Momentum 
  // -- Advective flux
  double gamma = cp/cv;
  double kinetic_energy = (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]) / 2.;
  double P = (E - kinetic_energy * rho - rho*g*z[0]) * (gamma - 1.);
  double F_adv_momentum[3][3] = {{rho*u[0]*u[0] + P, rho*u[0]*u[1],     rho*u[0]*u[2]},
                                 {rho*u[1]*u[0],     rho*u[1]*u[1] + P, rho*u[1]*u[2]}, 
                                 {rho*u[2]*u[0],     rho*u[2]*u[1],     rho*u[2]*u[02] + P}};
  // -- Diffusive Flux
  double Fu[6] = {mu*(du[0][0] * (2 + lambda) + lambda * (du[1][1] + du[2][2])),
                  mu*(du[0][1] + du[1][0]), 
                  mu*(du[0][2] + du[2][0]), 
                  mu*(du[1][1] * (2 + lambda) + lambda * (du[0][0] + du[2][2])),
                  mu*(du[1][2] + du[2][1]), 
                  mu*(du[2][2] * (2 + lambda) + lambda * (du[0][0] + du[1][1]))
                 };

  const int Fuviscidx[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}}; // symmetric matrix indices
  double F_dif_momentum[3][3];
  for (int j=0; j<3; j++) 
    for (int k=0; k<3; k++) 
      F_dif_momentum[j][k] = Fu[Fuviscidx[j][k]];

  // Total Energy
  // -- Advective flux
  double F_adv_energy[3] = {(E + P)*u[0], (E + P)*u[1], (E + P)*u[2]};

  // -- Diffusive Flux  
  double dT[3] = {(dE[0]/rho - E*drho[0]/(rho*rho) - (u[0]*du[0][0] + u[1]*du[1][0] + u[2]*du[2][0]))    /cv,
                  (dE[1]/rho - E*drho[1]/(rho*rho) - (u[0]*du[0][1] + u[1]*du[1][1] + u[2]*du[2][1]))    /cv,
                  (dE[2]/rho - E*drho[2]/(rho*rho) - (u[0]*du[0][2] + u[1]*du[1][2] + u[2]*du[2][2]) - g)/cv
                 };
  double F_dif_energy[3] = {u[0]*Fu[0] + u[1]*Fu[1] + u[2]*Fu[2] + k*dT[0], 
                            u[0]*Fu[1] + u[1]*Fu[3] + u[2]*Fu[4] + k*dT[1], 
                            u[0]*Fu[2] + u[1]*Fu[4] + u[2]*Fu[5] + k*dT[2] 
                           };
  
  // Populate Flux
  double f[3][5] = {{0.}};
  // -- Density
  for (int i=0; i<3; i++) f[i][0] += F_adv_density[i];

  // -- Momentum
  for (int i=0; i<3; i++) 
    for (int j=0; j<3; j++) {
      f[i][j+1] += F_adv_momentum[j][i];
      f[i][j+1] -= F_dif_momentum[j][i];
    }

  // -- Energy
  for (int i=0; i<3; i++) {
      f[i][4] += F_adv_energy[i];
      f[i][4] -= F_dif_energy[i];
  }

  for (int i=0; i<5; i++) f2[i] = f[2][i];
}

void compute_dF2_dz(double df2_dz[5], double *t, double *x, double *y, double *z,
                   double lambda, double mu, double k, double cv, double cp, double g) {
  double f2[5];
  for (int i=0; i<5; i++) {
    double f2_[5] = {0.}; f2_[i] = 1.;
    __enzyme_autodiff((void *)computeF2,
                      f2, f2_,
                      enzyme_const, t,
                      enzyme_const, x,
                      enzyme_const, y,
                      z, &df2_dz[i],
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
  double wdetJ = 1.;
  double force[5];
  // Zero force so all future terms can safely sum into it
  for (int j=0; j<5; j++) force[j] = 0.;
  // -------------------------------------------------------------------------
  // q_dot
  // -------------------------------------------------------------------------
  double q_dot[5] = {0.};
  computeQdot(q_dot, time, x, y, z);

  // Add Qdot to the forcing term
  for (int j=0; j<5; j++) force[j] += wdetJ*q_dot[j];

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
  Pascal          = kilogram / (meter * second*second);
  J_per_kg_K      = (meter*meter) / (second*second * Kelvin);
  m_per_squared_s = meter / (second*second);
  W_per_m_K       = kilogram * meter / (pow(second,3) * Kelvin);

  // -- Unit conversion
  cv *= J_per_kg_K;
  cp *= J_per_kg_K;
  g  *= m_per_squared_s;
  mu *= Pascal * second;
  k  *= W_per_m_K;
  double gamma  = cp / cv;

  // -- Compute grad(flux)
  double grad_F[3][5] = {{0.}}; // Must initialize
  compute_dF0_dx(grad_F[0], time, x, y, z, lambda, mu, k, cv, cp, g);
  compute_dF1_dy(grad_F[1], time, x, y, z, lambda, mu, k, cv, cp, g);
  compute_dF2_dz(grad_F[2], time, x, y, z, lambda, mu, k, cv, cp, g);
  // Print output
  printf("\n-------------------------------------------------------------------------\n");
  printf("grad(flux):");
  printf("\n-------------------------------------------------------------------------");
  for (int i=0; i<3; i++) {
    printf("\nDerivative in direction %d:\n", i);
    for (int j=0; j<5; j++) printf("%f\t", grad_F[i][j]);
    printf("\n");
  }

  // -- Compute div(flux)
  double div_f[5] = {0.};
  for (int j=0; j<5; j++) for (int k=0; k<3; k++) div_f[j] += grad_F[k][j];

  // Add div(flux) to the forcing term
  for (int j=0; j<5; j++) force[j] += wdetJ*div_f[j];

  // -------------------------------------------------------------------------
  // Body force
  // -------------------------------------------------------------------------
  // Compute state variables                
  double q[5];
  ExactSolution(q, time, x, y, z);
  double rho = q[0];
  force[3] += wdetJ*rho*g;

  // -------------------------------------------------------------------------
  // Print Forcing Terms
  // -------------------------------------------------------------------------
  // Print output
  printf("\n-------------------------------------------------------------------------\n");
  printf("force:");
  printf("\n-------------------------------------------------------------------------\n");
  for (int j=0; j<5; j++) printf("%f\t", force[j]);
  printf("\n\n");

  return 0;
}
// ***************************************************************************

/*
clang test022.c -Xclang -load -Xclang /home/leila/Enzyme/enzyme/build12DHB/Enzyme/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops

Output:
-------------------------------------------------------------------------
Flux:
-------------------------------------------------------------------------
Flux in direction 0:
82.900000	-162693.796709	103.275085	123.668264	-162612.275110	

Flux in direction 1:
103.300000	103.275085	-162647.968064	154.146018	-202627.748570	

Flux in direction 2:
123.700000	123.668264	154.146018	-162592.054844	-242643.222031

-------------------------------------------------------------------------
grad(flux):
-------------------------------------------------------------------------
Derivative in direction 0:
1.000000	-1962.601219	2.000383	3.000547	-1958.600729	

Derivative in direction 1:
7.000000	7.001751	-11773.008432	9.256028	-13733.141883	

Derivative in direction 2:
13.000000	13.003838	13.647028	-54137.027758	-74083.800481	

-------------------------------------------------------------------------
force:
-------------------------------------------------------------------------
35.000000	-1914.595630	-11715.361021	27256.128818	-89705.543093

*/
