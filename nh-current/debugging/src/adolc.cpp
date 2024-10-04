// Standard Neo-Hookean formulation in initial configuration with ADOL-C

#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>
#include <adolc/adolc.h>

using namespace std;
static int m = 1, // number of dependent variables
           n = 6; // number of independent variables

#include "../include/nh-common.h"
#include "../include/nh-adolc.h"

int main() {
  const double mu = 1., lambda = 1.;

  double dXdx_initial[3][3] = {
    {0.0702417, 0.4799115, 0.3991242},
    {0.6756593, 0.0633284, 0.0959267},
    {0.2241923, 0.0281781, 0.0917613}
  };

  double dudX[3][3] = {
    {0.1425560,  0.115120,  0.551640},
    {0.0591922,  0.123535,  0.166572},
    {0.1617210,  0.478828,  0.646217}
  };

  // Grad_u = du/dx_initial = du/dX * dX/dx_initial
  double Grad_u[3][3];
  MatMatMult(1.0, dudX, dXdx_initial, Grad_u);

  // Deformation Gradient : F = I + Grad_u
  const double F[3][3] = {
    {Grad_u[0][0] + 1, Grad_u[0][1],     Grad_u[0][2]    },
    {Grad_u[1][0],     Grad_u[1][1] + 1, Grad_u[1][2]    },
    {Grad_u[2][0],     Grad_u[2][1],     Grad_u[2][2] + 1}
  };

  // F^{-1} = dx_initial/dx
  double F_inv[3][3];
  const double Jm1  = MatDetAM1(Grad_u);
  const double detF = Jm1 + 1.;
  MatInverse(F, detF, F_inv);

  // dXdx = dX/dx_current = dX/dx_initial * F^{-1}
  double dXdx[3][3];
  MatMatMult(1.0, dXdx_initial, F_inv, dXdx);

  // Green Euler Strain tensor (e)
  double e_sym[6];
  GreenEulerStrain(Grad_u, e_sym);

  // ------------------------------------------------------------------------
  // Automatic Differentiation
  // ------------------------------------------------------------------------
  // Initialize passive variables
  auto e_p = new double[n];
  for (int i=0; i<n; i++) e_p[i] = e_sym[i];

  // gradPsi = dPsi/de with ADOLC
  double gradPsi_sym[6] = {0.}, gradPsi[3][3];
  ComputeGradPsi(gradPsi_sym, e_p, lambda, mu);
  SymmetricMatUnpack(gradPsi_sym, gradPsi);

  // b = 2 e + I
  double b_sym[6], b[3][3];
  for (int j = 0; j < 6; j++) b_sym[j] = 2. * e_p[j] + (j < 3);
  SymmetricMatUnpack(b_sym, b);

  // tau = dPsi/de * b
  double tau_sym[6], tau[3][3];
  MatMatMult(1.0, gradPsi, b, tau);
  SymmetricMatPack(tau, tau_sym);

  // grad_du = ddu/dX * dX/dx
  double grad_du[3][3];
  double ddudX[3][3] = {
      {0.23617975,  0.60250516,  0.1717169},
      {0.86615524,  0.3365063,   0.17371375},
      {0.0441905,   0.16762188,  0.45047923}
  };
  MatMatMult(1.0, ddudX, dXdx, grad_du);

  // de = db / 2 = (grad_du b + b (grad_du)^T) / 2
  double de_sym[6], de[3][3];
  GreenEulerStrain_fwd(grad_du, b, de_sym);
  SymmetricMatUnpack(de_sym, de);

  // The hessian of Psi (d2Psi/de2) with ADOLC
  double hessPsi_curr[6][6] = {{0.}};
  ComputeHessianPsi(hessPsi_curr, e_p, lambda, mu);

  // dGradPsi = hessPsi : de
  double dGradPsi[3][3], dGradPsi_sym[6] = {0.};
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) dGradPsi_sym[i] += hessPsi_curr[i][j] * de_sym[j];
  SymmetricMatUnpack(dGradPsi_sym, dGradPsi);

  // dtau_1 = dGradPsi b
  double dtau_1[3][3];
  MatMatMult(1.0, dGradPsi, b, dtau_1);

  // dtau_2 = 2 (gradPsi * de ) : de
  double gradPsi_de[3][3], dtau_2[3][3];
  MatMatMult(2., gradPsi, de, gradPsi_de);
  for (int i=0; i<3; i++) for (int j=0; j<3; j++) dtau_2[i][j] = gradPsi_de[i][j] * de[i][j];

  // dtau = dtau_1 + dtau_2 (TODO: all members are ~.3 smaller!)
  double dtau[3][3], dtau_sym[6];
  MatMatAdd(1., dtau_1, 1., dtau_2, dtau);
  SymmetricMatPack(dtau, dtau_sym);

  // ------------------------------------------------------------------------
  // More info for debugging
  // ------------------------------------------------------------------------
  printf("\n\ne =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", e_sym[i]);

  double E_sym[6];
  Compute_E_symmetric(Grad_u, E_sym);
  printf("\n\nE =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", E_sym[i]);

  printf("\n\nStrain Energy from e = ");
  printf(" %.12lf", StrainEnergy(e_sym, lambda, mu));

  printf("\nStrain Energy from E = ");
  printf(" %.12lf", StrainEnergy(E_sym, lambda, mu));

  // S = dPsi/dE with ADOLC
  double S_sym_ad[6];
  auto Ep = new double[n];
  for (int i=0; i<n; i++) Ep[i] = E_sym[i];
  ComputeGradPsi(S_sym_ad, Ep, lambda, mu);
  printf("\n\nS from AD =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", S_sym_ad[i]);

  double S_sym_pb[6];
  PullBack_symmetric(Grad_u, tau_sym, S_sym_pb);
  printf("\n\nS from pull-back =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", S_sym_pb[i]);

  printf("\n\ntau from AD =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", tau_sym[i]);

  double tau_sym_pf[6];
  PushForward_symmetric(Grad_u, S_sym_ad, tau_sym_pf);
  printf("\n\ntau from push-forward =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", tau_sym_pf[i]);

  printf("\n\nde =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", de_sym[i]);

  // Grad_du = ddu/dX * dX/dx_initial
  double Grad_du[3][3], dE_sym[6];
  MatMatMult(1.0, ddudX, dXdx_initial, Grad_du);
  RatelGreenLagrangeStrain_fwd(Grad_du, F, dE_sym);
  printf("\n\ndE =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", dE_sym[i]);

  // dS = dS/dE with ADOLC
  double dS_sym_ad[6];
  double hessPsi_init[6][6] = {{0.}};
  ComputeHessianPsi(hessPsi_init, Ep, lambda, mu);
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) dS_sym_ad[i] += hessPsi_init[i][j] * dE_sym[j];
  printf("\n\ndS from AD =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", dS_sym_ad[i]);

  double dS_sym_pb[6] = {0.};
  dPullBack_symmetric(F_inv, Grad_du, tau_sym, dtau_sym, dS_sym_pb);
  printf("\n\ndS from pull-back =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", dS_sym_pb[i]);

  printf("\n\ndtau from AD =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", dtau_sym[i]);

  double dtau_sym_pf[6];
  dPushForward_symmetric(F, Grad_du, S_sym_ad, dS_sym_ad, dtau_sym_pf);
  printf("\n\ndtau from push-forward =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", dtau_sym_pf[i]);

  printf("\n\n");
  return 0;
}
