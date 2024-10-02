// Standard Neo-Hookean formulation  with Enzyme
#include "nh-enzyme.h"

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
  // dXdx = dX/dx = dX/dx_initial * F^{-1}
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
  // X is ref coordinate [-1,1]^3; x is physical coordinate in current configuration
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
  PullBack_symmetric(Grad_u, e_sym, E_sym);
  printf("\n\nE_sym from pull-back =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", E_sym[i]);

  double dE_sym[6];
  PullBack_symmetric(Grad_u, de_sym, dE_sym);
  printf("\n\ndE_sym from pull-back =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", dE_sym[i]);

  double S_sym_pb[6];
  PullBack_symmetric(Grad_u, tau_sym, S_sym_pb);
  printf("\n\nS_sym from pull-back =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", S_sym_pb[i]);

  printf("\n\nStrain Energy from pull-back = ");
  printf(" %.6lf", StrainEnergy_Enzyme(E_sym, lambda, mu));

  double S_sym_ad[6];
  S_sym_Enzyme(lambda, mu, E_sym, S_sym_ad);
  printf("\n\nS_sym from AD =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", S_sym_ad[i]);

  double dS_sym_ad[6];
  dS_fwd_Enzyme(lambda, mu, E_sym, dE_sym, S_sym_ad, dS_sym_ad);
  printf("\n\ndS_sym from AD =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", dS_sym_ad[i]);

  double tau_sym_pf[6];
  PushForward_symmetric(Grad_u,  S_sym_ad, tau_sym_pf);
  printf("\n\ntau_sym from push-forward =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", tau_sym_pf[i]);

  printf("\n\nStrain Energy from e = ");
  printf(" %.6lf", StrainEnergy_Enzyme(e_sym, lambda, mu));

  printf("\n\ntau from AD =");
  for (int i=0; i<6; i++) printf("\n\t%.12lf", tau_sym[i]);

  printf("\n\ndtau from AD=");
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


dtau         =

        2.662047097805
        2.515592373839
        1.869278025008
        0.750807688616
        0.714042268606
        0.712369842831
*/

/*
push-forward: a = F A F^T
pull-back: A = F^{-1} a F^{-T}
*/
