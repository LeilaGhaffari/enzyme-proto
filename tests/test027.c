                             

// Use autodiff for computing gradients
//   MMS in libCEED/examples/fluids follows this test

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "test026.h"


int main() {
  // Declarations
  double x[] = {.5, 2.5, 5.};
  double time[1] = {.2};
  double wdetJ = 1.;
  double force[5];

  // Zero force so all future terms can safely sum into it
  for (int j=0; j<5; j++) force[j] = 0.;

  // -------------------------------------------------------------------------
  // q_diff
  // -------------------------------------------------------------------------
  double q[5];
  double q_dot[5] = {0.};
  double dq[5][3] = {{0.}};
  q_diff(q, q_dot, dq, time, x);

  // Print q_dot
  printf("\n-------------------------------------------------------------------------\n");
  printf("q_dot:");
  printf("\n-------------------------------------------------------------------------\n");
  for (int j=0; j<5; j++) printf("%f\t", q_dot[j]);
  printf("\n");

  // Print dq
  printf("\n-------------------------------------------------------------------------\n");
  printf("dq:");
  printf("\n-------------------------------------------------------------------------");
  for (int i=0; i<3; i++) {
    printf("\nDerivative in direction %d:\n", i);
    for (int j=0; j<5; j++) printf("%f\t", dq[j][i]);
    printf("\n");
  }

  // Add Qdot to the forcing term
  for (int j=0; j<5; j++) force[j] += wdetJ*q_dot[j];

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

  // -------------------------------------------------------------------------
  // Flux
  // -------------------------------------------------------------------------
  double F[3][5] = {{0.}}; // Must initialize
  computeF0(F[0], time, x, lambda, mu, k, cv, cp, g);
  computeF1(F[1], time, x, lambda, mu, k, cv, cp, g);
  computeF2(F[2], time, x, lambda, mu, k, cv, cp, g);
  // Print output
  printf("\n-------------------------------------------------------------------------\n");
  printf("Flux:");
  printf("\n-------------------------------------------------------------------------");
  for (int i=0; i<3; i++) {
    printf("\nFlux in direction %d:\n", i);
    for (int j=0; j<5; j++) printf("%f\t", F[i][j]);
    printf("\n");
  }

  // -------------------------------------------------------------------------
  // Spacial derivative of Flux
  // -------------------------------------------------------------------------
  double grad_F[3][5] = {{0.}};
  compute_dF0_dx(grad_F[0], time, x, lambda, mu, k, cv, cp, g);
  compute_dF1_dy(grad_F[1], time, x, lambda, mu, k, cv, cp, g);
  compute_dF2_dz(grad_F[2], time, x, lambda, mu, k, cv, cp, g);
  // Print output
  printf("\n-------------------------------------------------------------------------\n");
  printf("grad(flux):");
  printf("\n-------------------------------------------------------------------------");
  for (int i=0; i<3; i++) {
    printf("\nDerivative in direction %d:\n", i);
    for (int j=0; j<5; j++) printf("%f\t", grad_F[i][j]);
    printf("\n");
  }

  // -------------------------------------------------------------------------
  // div(Flux)
  // -------------------------------------------------------------------------
  double div_f[5] = {0.};
  for (int j=0; j<5; j++) for (int k=0; k<3; k++) div_f[j] += grad_F[k][j];

  // Add div(flux) to the forcing term
  for (int j=0; j<5; j++) force[j] += wdetJ*div_f[j];

  // -------------------------------------------------------------------------
  // Body force
  // -------------------------------------------------------------------------                
  double rho = q[0];

  // Add body force to the forcing term
  force[3] += wdetJ*rho*g;

  // -------------------------------------------------------------------------
  // Print Forcing Terms
  // -------------------------------------------------------------------------
  printf("\n-------------------------------------------------------------------------\n");
  printf("force:");
  printf("\n-------------------------------------------------------------------------\n");
  for (int j=0; j<5; j++) printf("%f\t", force[j]);
  printf("\n\n");

  return 0;
}
// ***************************************************************************

/*

Output:
-------------------------------------------------------------------------
q_dot:
-------------------------------------------------------------------------
14.000000	28.000000	42.000000	56.000000	70.000000	

-------------------------------------------------------------------------
dq:
-------------------------------------------------------------------------
Derivative in direction 0:
1.000000	2.000000	3.000000	4.000000	5.000000	

Derivative in direction 1:
6.000000	7.000000	8.000000	9.000000	10.000000	

Derivative in direction 2:
11.000000	12.000000	13.000000	14.000000	15.000000	

-------------------------------------------------------------------------
Flux:
-------------------------------------------------------------------------
Flux in direction 0:
103.300000	-162681.538364	154.130717	179.555052	-202669.586479	

Flux in direction 1:
123.700000	154.130717	-162625.647407	215.063961	-242693.258650	

Flux in direction 2:
144.100000	179.555052	215.063961	-162559.671876	-282716.930461	

-------------------------------------------------------------------------
grad(flux):
-------------------------------------------------------------------------
Derivative in direction 0:
2.000000	-1962.148435	4.863979	6.296025	-3922.190909	

Derivative in direction 1:
8.000000	9.257916	-11772.054837	11.767718	-15696.907074	

Derivative in direction 2:
14.000000	14.476215	14.943049	-54136.789826	-84057.035506	

-------------------------------------------------------------------------
force:
-------------------------------------------------------------------------
38.000000	-1910.414304	-11710.247809	27262.173917	-103606.133490

*/
