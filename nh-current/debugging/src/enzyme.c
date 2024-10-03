// Standard Neo-Hookean formulation  with Enzyme
#include "../include/nh-enzyme.h"

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

  // tau = dPsi b with Enzyme-AD
  double tau_sym[6];
  tau_sym_Enzyme(lambda, mu, e_sym, tau_sym);

  // grad_du = ddu/dX * dX/dx
  double grad_du[3][3];
  double ddudX[3][3] = {
      {0.23617975,  0.60250516,  0.1717169},
      {0.86615524,  0.3365063,   0.17371375},
      {0.0441905,   0.16762188,  0.45047923}
  };
  MatMatMult(1.0, ddudX, dXdx, grad_du);

  // b = 2 e + I
  double b_sym[6], b[3][3];
  for (int j = 0; j < 6; j++) b_sym[j] = 2 * e_sym[j] + (j < 3);
  SymmetricMatUnpack(b_sym, b);

  // de = db / 2 = (grad_du b + b (grad_du)^T) / 2
  double de_sym[6];
  GreenEulerStrain_fwd(grad_du, b, de_sym);

  // dtau with Enzyme-AD
  double dtau_sym[6];
  dtau_fwd_Enzyme(lambda, mu, e_sym, de_sym, tau_sym, dtau_sym);

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

  // S = dPsi/dE with Enzyme-AD
  double S_sym_ad[6];
  S_sym_Enzyme(lambda, mu, E_sym, S_sym_ad);
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

  // dS = dS/dE with Enzyme-AD
  double dS_sym_ad[6];
  dS_fwd_Enzyme(lambda, mu, E_sym, dE_sym, S_sym_ad, dS_sym_ad);
  printf("\n\ndS from AD =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", dS_sym_ad[i]);

  double dS_sym_pb[6] = {0.};
  dPullBack_symmetric(F_inv, Grad_du, tau_sym, dtau_sym, dS_sym_pb);
  printf("\n\ndS from pull-back =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", dS_sym_pb[i]);

  printf("\n\ndtau from AD =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", dtau_sym[i]);

  double dtau_sym_pf[6]; // note: we don't use S_sym_ad here since it changed with dS_fwd_Enzyme()
  dPushForward_symmetric(F, Grad_du, S_sym_pb, dS_sym_ad, dtau_sym_pf);
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
        0.603042145694
        0.502393311590
        0.518369381113
        0.039368083122
        0.192988146562
        0.071117045859

S from pull-back =
        0.603042145694
        0.502393311590
        0.518369381113
        0.039368083122
        0.192988146562
        0.071117045859

tau from AD =
        0.968541661897
        0.580219129800
        1.092963392519
        0.250640282229
        0.731414543779
        0.252397811536

tau from push-forward =
        0.968541661897
        0.580219129800
        1.092963392519
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
        1.349682793048
        1.751663099730
        1.340753124766
        -0.027597262846
        -0.553363888159
        -0.172700195711

dS from pull-back =
        1.349682793048
        1.751663099730
        1.340753124766
        -0.027597262846
        -0.553363888159
        -0.172700195711

dtau from AD =
        2.662047097805
        2.515592373839
        1.869278025008
        0.750807688616
        0.714042268606
        0.712369842831

dtau from push-forward =
        2.662047097805
        2.515592373839
        1.869278025008
        0.750807688616
        0.714042268606
        0.712369842831
*/
