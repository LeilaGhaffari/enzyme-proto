// ***************************************************************************
// Exact Solution
// ***************************************************************************
int ExactSolution(double q[], double *time, double X[]) {
  // -- Time
  double t = time[0];
  // -- Coordinates
  double x = X[0];
  double y = X[1];
  double z = X[2];

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

// -- q_diff
void q_diff(double q[5], double q_dot[5], double dq[5][3], double *t, double x[]) {
  for (int i=0; i<5; i++) {
      double q_[5] = {0.}; q_[i] = 1.;
      __enzyme_autodiff((void *)ExactSolution,
                        q, q_,
                        t, &q_dot[i],
                        x, dq[i]);
  }
}

// -- dFlux[0]/dx
void computeF0(double f0[5], double *t, double x[], 
               double lambda, double mu, double k, double cv, double cp, double g) {              
  // Compute state variables and the derivatives
  double q[5];
  double q_dot[5] = {0.};
  double dq[5][3] = {{0.}};
  q_diff(q, q_dot, dq, t, x);

  // -- q
  double rho = q[0];
  double u[3] = {q[1]/rho, q[2]/rho, q[3]/rho};
  double E = q[4];
  // -- dq
  double drho[3]   = {dq[0][0], dq[0][1], dq[0][2]};
  double dU[3][3] = {{dq[1][0], dq[1][1], dq[1][2]},
                     {dq[2][0], dq[2][1], dq[2][2]},
                     {dq[3][0], dq[3][1], dq[3][2]}
                    };
  double dE[3] =     {dq[4][0], dq[4][1], dq[4][2]};
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
  double P = (E - kinetic_energy * rho - rho*g*x[2]) * (gamma - 1.);
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
  // -- Zero f0
  for (int j=0; j<5; j++) f0[j] = 0.;

  // -- Density
  f0[0] += F_adv_density[0];

  // -- Momentum
  for (int j=0; j<3; j++) {
    f0[j+1] += F_adv_momentum[j][0];
    f0[j+1] -= F_dif_momentum[j][0];
  }

  // -- Energy
  f0[4] += F_adv_energy[0];
  f0[4] -= F_dif_energy[0];

}

void compute_dF0_dx(double df0_dx[5], double *t, double x[],
                    double lambda, double mu, double k, double cv, double cp, double g) {
  double f[5];
  for (int i=0; i<5; i++) {
    double f_[5] = {0.}; f_[i] = 1.;
    double df[3] = {0.};
    __enzyme_autodiff((void *)computeF0,
                      f, f_,
                      enzyme_const, t,
                      x, df,
                      enzyme_const, lambda,
                      enzyme_const, mu,
                      enzyme_const, k,
                      enzyme_const, cv,
                      enzyme_const, cp,
                      enzyme_const, g);
    df0_dx[i] = df[0];
  }
}

// -- dFlux[1]/dy
void computeF1(double f1[5], double *t, double x[], 
               double lambda, double mu, double k, double cv, double cp, double g) {              
  // Compute state variables and the derivatives
  double q[5];
  double q_dot[5] = {0.};
  double dq[5][3] = {{0.}};
  q_diff(q, q_dot, dq, t, x);

  // -- q
  double rho = q[0];
  double u[3] = {q[1]/rho, q[2]/rho, q[3]/rho};
  double E = q[4];
  // -- dq
  double drho[3]   = {dq[0][0], dq[0][1], dq[0][2]};
  double dU[3][3] = {{dq[1][0], dq[1][1], dq[1][2]},
                     {dq[2][0], dq[2][1], dq[2][2]},
                     {dq[3][0], dq[3][1], dq[3][2]}
                    };
  double dE[3] =     {dq[4][0], dq[4][1], dq[4][2]};
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
  double P = (E - kinetic_energy * rho - rho*g*x[2]) * (gamma - 1.);
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
  // -- Zero f0
  for (int j=0; j<5; j++) f1[j] = 0.;

  // -- Density
  f1[0] += F_adv_density[1];

  // -- Momentum
  for (int j=0; j<3; j++) {
    f1[j+1] += F_adv_momentum[j][1];
    f1[j+1] -= F_dif_momentum[j][1];
  }

  // -- Energy
  f1[4] += F_adv_energy[1];
  f1[4] -= F_dif_energy[1];
}

void compute_dF1_dy(double df1_dy[5], double *t, double x[],
                   double lambda, double mu, double k, double cv, double cp, double g) {
  double f[5];
  for (int i=0; i<5; i++) {
    double f_[5] = {0.}; f_[i] = 1.;
    double df[3] = {0.};
    __enzyme_autodiff((void *)computeF1,
                      f, f_,
                      enzyme_const, t,
                      x, df,
                      enzyme_const, lambda,
                      enzyme_const, mu,
                      enzyme_const, k,
                      enzyme_const, cv,
                      enzyme_const, cp,
                      enzyme_const, g);
    df1_dy[i] = df[1];
  }
}

// -- dFlux[2]/dz
void computeF2(double f2[5], double *t, double x[], 
               double lambda, double mu, double k, double cv, double cp, double g) {              
  // Compute state variables and the derivatives
  double q[5];
  double q_dot[5] = {0.};
  double dq[5][3] = {{0.}};
  q_diff(q, q_dot, dq, t, x);

  // -- q
  double rho = q[0];
  double u[3] = {q[1]/rho, q[2]/rho, q[3]/rho};
  double E = q[4];
  // -- dq
  double drho[3]   = {dq[0][0], dq[0][1], dq[0][2]};
  double dU[3][3] = {{dq[1][0], dq[1][1], dq[1][2]},
                     {dq[2][0], dq[2][1], dq[2][2]},
                     {dq[3][0], dq[3][1], dq[3][2]}
                    };
  double dE[3] =     {dq[4][0], dq[4][1], dq[4][2]};
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
  double P = (E - kinetic_energy * rho - rho*g*x[2]) * (gamma - 1.);
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
  // -- Zero f0
  for (int j=0; j<5; j++) f2[j] = 0.;

  // -- Density
  f2[0] += F_adv_density[2];

  // -- Momentum
  for (int j=0; j<3; j++) {
    f2[j+1] += F_adv_momentum[j][2];
    f2[j+1] -= F_dif_momentum[j][2];
  }

  // -- Energy
  f2[4] += F_adv_energy[2];
  f2[4] -= F_dif_energy[2];
}

void compute_dF2_dz(double df2_dz[5], double *t, double x[],
                   double lambda, double mu, double k, double cv, double cp, double g) {
  double f[5];
  for (int i=0; i<5; i++) {
    double f_[5] = {0.}; f_[i] = 1.;
    double df[3] = {0.};
    __enzyme_autodiff((void *)computeF2,
                      f, f_,
                      enzyme_const, t,
                      x, df,
                      enzyme_const, lambda,
                      enzyme_const, mu,
                      enzyme_const, k,
                      enzyme_const, cv,
                      enzyme_const, cp,
                      enzyme_const, g);
    df2_dz[i] = df[2];
  }
}
