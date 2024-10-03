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

  // Compute the Deformation Gradient : F = I + Grad_u
  const double F[3][3] = {
    {Grad_u[0][0] + 1, Grad_u[0][1],     Grad_u[0][2]    },
    {Grad_u[1][0],     Grad_u[1][1] + 1, Grad_u[1][2]    },
    {Grad_u[2][0],     Grad_u[2][1],     Grad_u[2][2] + 1}
  };

  // Compute F^{-1}
  double F_inv[3][3];
  const double Jm1  = MatDetAM1(Grad_u);
  const double detF = Jm1 + 1.;

  MatInverse(F, detF, F_inv);

  // x is current config coordinate system
  // dXdx = dX/dx_current = dX/dx_initial * F^{-1}
  // Note that F^{-1} = dx_initial/dx
  double dXdx[3][3];
  MatMatMult(1.0, dXdx_initial, F_inv, dXdx);

  // Compute Green Euler Strain tensor (e)
  double e_sym[6];
  GreenEulerStrain(Grad_u, e_sym);

  // ------------------------------------------------------------------------
  // Automatic Differentiation
  // ------------------------------------------------------------------------
  // tau = dPsi b
  double tau_sym[6];
  tau_sym_Enzyme(lambda, mu, e_sym, tau_sym);

  // Compute grad_du = ddu/dX * dX/dx
  double grad_du[3][3];
  double ddudX[3][3] = {
      {0.23617975,  0.60250516,  0.1717169},
      {0.86615524,  0.3365063,   0.17371375},
      {0.0441905,   0.16762188,  0.45047923}
  };
  MatMatMult(1.0, ddudX, dXdx, grad_du);

  // Compute b = 2 e + I
  double b_sym[6], b[3][3];
  for (int j = 0; j < 6; j++) b_sym[j] = 2 * e_sym[j] + (j < 3);
  SymmetricMatUnpack(b_sym, b);

  // Compute de = db / 2 = (grad_du b + b (grad_du)^T) / 2
  double de_sym[6];
  GreenEulerStrain_fwd(grad_du, b, de_sym);

  // Compute dtau with Enzyme-AD
  double dtau_sym[6];
  dtau_fwd_Enzyme(lambda, mu, e_sym, de_sym, tau_sym, dtau_sym);

  // ------------------------------------------------------------------------
  // Pull-back for debugging
  // ------------------------------------------------------------------------
  double E_sym[6];
  Compute_E_symmetric(Grad_u, E_sym);
  printf("\n\nE =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", E_sym[i]);

  printf("\n\nStrain Energy from e = ");
  printf(" %.12lf", StrainEnergy(e_sym, lambda, mu));

  printf("\nStrain Energy from E = ");
  printf(" %.12lf", StrainEnergy(E_sym, lambda, mu));

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

  // Compute Grad_du = ddu/dX * dX/dx_initial
  double Grad_du[3][3], dE_sym[6];
  MatMatMult(1.0, ddudX, dXdx_initial, Grad_du);
  RatelGreenLagrangeStrain_fwd(Grad_du, F, dE_sym);
  printf("\n\ndE =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", dE_sym[i]);

  double dS_sym_ad[6];
  dS_fwd_Enzyme(lambda, mu, E_sym, dE_sym, S_sym_ad, dS_sym_ad);
  printf("\n\ndS from AD =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", dS_sym_ad[i]);

  double dtau_sym_pf[6];
  dPushForward_symmetric(F, Grad_du, S_sym_pb, dS_sym_ad, dtau_sym_pf);
  printf("\n\ndtau from push-forward =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", dtau_sym_pf[i]);

  printf("\n\ndtau from AD =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", dtau_sym[i]);
  printf("\n\n");

  return 0;
}

/*

Strain Energy =  0.338798

tau =

        0.968541661897
        0.580219129800
        1.092963392519
        0.250640282229
        0.731414543779
        0.252397811536


dtau =

        2.662047097805
        2.515592373839
        1.869278025008
        0.750807688616
        0.714042268606
        0.712369842831
*/
