// Standard Neo-Hookean formulation in current configuration with ADOL-C

#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>
#include <adolc/adolc.h>

using namespace std;
static int m = 1, n = 6;

#include "../include/nh-current.hpp"
#include "../include/nh-current-adolc.hpp"



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

  // Compute dXdx = dX/dx = dX/dx_initial * F^{-1}
  double dXdx[3][3];
  MatMatMult(1.0, dXdx_initial, F_inv, dXdx);

  // Compute Green Euler Strain tensor (e)
  double e_sym[n];
  GreenEulerStrain(Grad_u, e_sym);

  // ------------------------------------------------------------------------
  // Automatic Differentiation
  // ------------------------------------------------------------------------
  // Initialize passive variables
  auto e_p = new double[n];
  for (int i=0; i<n; i++) e_p[i] = e_sym[i];

  // Compute gradPsi = dPsi/de
  double gradPsi_sym[6] = {0.}, gradPsi[3][3];
  ComputeGradPsi(gradPsi_sym, e_p, lambda, mu);
  SymmetricMatUnpack(gradPsi_sym, gradPsi);

  // Compute b = 2 e + I
  double b_sym[6], b[3][3];
  for (int j = 0; j < 6; j++) b_sym[j] = 2. * e_p[j] + (j < 3);
  SymmetricMatUnpack(b_sym, b);

  // Compute tau = dPsi/de * b
  double tau_sym[6], tau[3][3];
  MatMatMult(1.0, gradPsi, b, tau);
  SymmetricMatPack(tau, tau_sym);

  // Compute grad_du = ddu/dX * dX/dx
  double grad_du[3][3];
  double ddudX[3][3] = {
      {0.23617975,  0.60250516,  0.1717169},
      {0.86615524,  0.3365063,   0.17371375},
      {0.0441905,   0.16762188,  0.45047923}
  };
  MatMatMult(1.0, ddudX, dXdx, grad_du);

  // Compute de = db / 2 = (grad_du b + b (grad_du)^T) / 2
  double de_sym[6];
  GreenEulerStrain_fwd(grad_du, b, de_sym);

  // Compute the hessian of Psi (d2Psi/de2)
  double hessPsi[6][6] = {{0.}};
  ComputeHessianPsi(hessPsi, e_p, lambda, mu);

  // Compute dGradPsi = hessPsi : de
  double dGradPsi[3][3], dGradPsi_sym[6] = {0.};
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) dGradPsi_sym[i] += hessPsi[i][j] * de_sym[j];
  SymmetricMatUnpack(dGradPsi_sym, dGradPsi);

  // Compute dPsi = gradPsi : de
  double dPsi = 0.;
  for (int i=0; i<n; i++) dPsi += gradPsi_sym[i] * de_sym[i];

  // Compute dtau
  // tau = gradPsi * b => Jac_tau = dtau/de = d(gradPsi * b)/de = hessPsi * b + 2 gradPsi * I
  // dtau = Jac_tau : de = b * (hessPsi : de) + 2 (gradPsi : de) * I = b * dGradPsi + 2 dPsi * I
  double dtau[3][3], dtau_sym[6];
  MatMatMult(1.0, b, dGradPsi, dtau);
  for (int i=0; i<3; i++) dtau[i][i] += 2. * dPsi;
  SymmetricMatPack(dtau, dtau_sym);

  // ------------------------------------------------------------------------
  // Print
  // ------------------------------------------------------------------------
  cout.precision(12);
  cout.setf(ios::fixed);
  cout << "\n Strain energy = " << StrainEnergy(e_sym, lambda, mu) << endl;
  cout << "\n tau =" << endl;
  for (int i=0; i<n; i++) cout << "\t" << tau_sym[i] << endl;
  cout << endl;
  cout << "\n dtau =" << endl << endl;
  for (int i=0; i<n; i++) cout << "\t" << dtau_sym[i] << endl;
  cout << endl;

  return 0;
}

/*
 Strain energy = 0.338798323727

 tau =
        0.968543642678
        0.580221110582
        1.092965373300
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
