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
  double F[3][3];
  DeformationGradient(Grad_u, F);

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
    double gradPsi_sym[6] = {0.}, gradPsi[3][3], hessPsi_curr[6][6] = {{0.}};
  // Initialize passive variables
  auto e_p = new double[n];
  for (int i=0; i<n; i++) e_p[i] = e_sym[i];

  // ADOL-C
  ComputeGradPsi(gradPsi_sym, e_p, lambda, mu); // gradient: dPsi/de
  ComputeHessianPsi(hessPsi_curr, e_p, lambda, mu); // hessian: d2Psi/de2
  // ------------------------------------------------------------------------
  // ------------------------------------------------------------------------

  // b = 2 e + I
  double b_sym[6], b[3][3];
  for (int j = 0; j < 6; j++) b_sym[j] = 2. * e_sym[j] + (j < 3);
  SymmetricMatUnpack(b_sym, b);

  // tau = dPsi/de * b
  double tau_sym[6], tau[3][3];
  SymmetricMatUnpack(gradPsi_sym, gradPsi);
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

  // dtau = (hessPsi : de) b + 2 gradPsi (I_4 : de)
  //      = dGradPsi b + 2 gradPsi de
  //      = dtau_1 + dtau_2

  // dGradPsi = hessPsi : de
  double dGradPsi[3][3], dGradPsi_sym[6] = {0.};
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) dGradPsi_sym[i] += hessPsi_curr[i][j] * de_sym[j];
  SymmetricMatUnpack(dGradPsi_sym, dGradPsi);

  // dtau_1 = dGradPsi b
  double dtau_1[3][3];
  MatMatMult(1.0, dGradPsi, b, dtau_1);

  // dtau_2 = 2 gradPsi de
  double dtau_2[3][3];
  MatMatMult(2., gradPsi, de, dtau_2);

  // dtau = dtau_1 + dtau_2
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
  PushForward_symmetric(F, S_sym_ad, tau_sym_pf);
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

/*
e =
        0.245019612050
        0.050858346002
        0.307230477361
        0.125320141114
        0.365707271890
        0.126198905768

E =
        0.356722230881
        0.053880729108
        0.192505475424
        0.105608348553
        0.355594178128
        0.150573971552

Strain Energy from e =  0.338798323727
Strain Energy from E =  0.338798323727

S from AD =
        0.603043653442
        0.502395201627
        0.518371210470
        0.039367933592
        0.192987413544
        0.071116775739

S from pull-back =
        0.603043653442
        0.502395201627
        0.518371210470
        0.039367933592
        0.192987413544
        0.071116775739

tau from AD =
        0.968543642678
        0.580221110582
        1.092965373300
        0.250640282229
        0.731414543779
        0.252397811536

tau from push-forward =
        0.968543642678
        0.580221110582
        1.092965373300
        0.250640282229
        0.731414543779
        0.252397811536

de =
        0.594073472562
        0.520846110579
        0.197688936164
        0.375403844308
        0.357021134303
        0.356184921415

dE =
        0.705072432429
        0.479848109568
        0.127687977309
        0.263932661797
        0.307099644410
        0.338054016477

dS from AD =
        1.349719621893
        1.751709039532
        1.340798936766
        -0.027601318917
        -0.553382183652
        -0.172707056257

dS from pull-back =
        1.349719621893
        1.751709039532
        1.340798936766
        -0.027601318917
        -0.553382183652
        -0.172707056257

dtau from AD =
        2.662096617697
        2.515641893730
        1.869327544899
        0.750807688616
        0.714042268606
        0.712369842831

dtau from push-forward =
        2.662096617697
        2.515641893730
        1.869327544899
        0.750807688616
        0.714042268606
        0.712369842831
*/