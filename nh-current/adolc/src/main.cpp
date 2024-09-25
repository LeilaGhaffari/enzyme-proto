// Standard Neo-Hookean formulation in current configuration with ADOL-C

#include <stdio.h>
# include <stdlib.h>
#include <math.h>
#include <string.h>
#include <adolc/adolc.h>

using namespace std;
static int m = 1, n = 6, d = 2, p = 6;

#include "../include/nh-current.hpp"
#include "../include/nh-current-adolc.hpp"



int main() {
  const double mu = 1., lambda = 1.;
  int size = binomi(p + d, d);

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
  double e_sym[n];
  GreenEulerStrain(Grad_u, e_sym);

  // ------------------------------------------------------------------------
  // Automatic Differentiation
  // ------------------------------------------------------------------------
  // Initialize passive variables
  auto e_p = new double[n];
  for (int i=0; i<n; i++) e_p[i] = e_sym[i];

  // First derivative (forward vector mode)
  auto tau_fwd = Kirchhofftau_sym_NeoHookean_AD_ADOLC(lambda, mu, e_p);

  // High order tensor
  // Tensor-AD evaluation
  auto tensor = dtau_ADOLC(e_p, lambda, mu);

  // Strain energy
  auto psi = tensor[0];

  // Populate dPsi/de
  int temp = 1;
  double dPsi_sym[6] = {0.};
  for (int i=0; i<n; i++) {
    dPsi_sym[i] = tensor[temp];
    if (i>2) dPsi_sym[i] /= 2.;
    temp += i + 2;
  }

  // Compute tau = dPsi/de * (2e + I)
  auto tau = tau_from_dPsi(dPsi_sym, e_p);

  // Populate d2Psi/de2
  temp = 1.;
  double d2Psi[6][6] = {{0.}};
  for (int i=0; i< n; i++) {
    for (int j=0; j<i+1; j++) {
      d2Psi[i][j] = tensor[temp+j+1];
      if (i != j) d2Psi[j][i] = d2Psi[i][j];
    }
    temp += i + 2;
  }
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) if (i > 2) d2Psi[i][j] /= 2.;


  // b = 2 e + I
  double b_sym[6], dPsi[3][3], b[3][3];
  for (int j = 0; j < 6; j++) b_sym[j] = 2. * e_p[j] + (j < 3);
  SymmetricMatUnpack(b_sym, b);
  SymmetricMatUnpack(dPsi_sym, dPsi);

  // Compute grad_du = ddu/dX * dX/dx
  // X is ref coordinate [-1,1]^3; x is physical coordinate in current configuration
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

  // Compute d2Psi_fwd
  double d2Psi_fwd_sym[6] = {0.}, d2Psi_fwd[3][3] = {{0.}};
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) d2Psi_fwd_sym[i] += d2Psi[i][j] * de_sym[j];
  SymmetricMatUnpack(d2Psi_fwd_sym, d2Psi_fwd);

  // Compute dtau_sym
  // tau = (dPsi / de) b => dtau = d(dPsi * b) = d2Psi * b + 2 dPsi
  double dtau[3][3] = {{0.}}, dtau_sym[6] = {0.};
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++){
      for (int k=0; k<3; k++) {
        dtau[i][j] += d2Psi_fwd[i][k] * b[k][j];
      }
      dtau[i][j] += 2. * dPsi[i][j]; // Do I need to multiply by de?
    }
  }
  SymmetricMatPack(dtau, dtau_sym);

  // ------------------------------------------------------------------------
  // Print
  // ------------------------------------------------------------------------
  cout.precision(12);
  cout.setf(ios::fixed);
  cout << "\n size = " << size << endl;
  cout << "\n Strain energy = " << psi << endl;
  cout << "\n Stress =" << endl;
  cout << "\n   Forward Vector Mode:" << endl << endl;
  for (int i=0; i<n; i++) cout << "\t" << tau_fwd[i] << endl;
  cout << endl;
  cout << "\n   Higher Order Tensor:" << endl << endl;
  for (int i=0; i<n; i++) cout << "\t" << tau[i] << endl;
  cout << endl;
  cout << "\n dtau =" << endl << endl;
  for (int i=0; i<n; i++) cout << "\t" << dtau_sym[i] << endl;
  cout << endl;

  return 0;
}

/*
 size = 28

 Strain energy = 0.338798323727

 Stress =

   Forward Vector Mode:

        0.968543642678
        0.580221110582
        1.092965373300
        0.250640282229
        0.731414543779
        0.252397811536


   Higher Order Tensor:

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
