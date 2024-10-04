#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <adolc/adolc.h>


#include "nh-initial.hpp"
#include "utils.hpp"
#include "nh-initial-adolc.hpp"

static const double mu = 1., lambda = 1.;

void f1_NeoHookeanInitial_AD_ADOLC(double dudX[3][3], double dXdx[3][3], double stored_grad_u[3][3],
                                   double stored_S_sym[6], double f1[3][3]) {
  double E_sym[6], S[3][3];

  // Compute grad_u = du/dX * dX/dx
  MatMatMult(1.0, dudX, dXdx, stored_grad_u);

  // Compute the Deformation Gradient : F = I + grad_u
  const double F[3][3] = {
      {stored_grad_u[0][0] + 1, stored_grad_u[0][1],     stored_grad_u[0][2]    },
      {stored_grad_u[1][0],     stored_grad_u[1][1] + 1, stored_grad_u[1][2]    },
      {stored_grad_u[2][0],     stored_grad_u[2][1],     stored_grad_u[2][2] + 1}
  };

  // E - Green-Lagrange strain tensor
  GreenLagrangeStrain(stored_grad_u, E_sym);

  // Compute Second Piola-Kirchhoff (S) with ADOL-C
  auto Ep = new double[6];
  for (int i=0; i<6; i++) Ep[i] = E_sym[i];
  ComputeGradPsi(stored_S_sym, Ep, lambda, mu);
  SymmetricMatUnpack(stored_S_sym, S);

  // Compute the First Piola-Kirchhoff f1 = P = F*S
  MatMatMult(1.0, F, S, f1);
}

void df1_NeoHookeanInitial_AD_ADOLC(void *ctx, double dXdx[3][3], double ddudX[3][3], double stored_grad_u[3][3],
                                    double stored_S_sym[6], double df1[3][3]) {

  double grad_du[3][3], E_sym[6], dE_sym[6], dS_sym[6], dS[3][3], S[3][3];
  HessianData *data = static_cast<HessianData *>(ctx);

  // Compute grad_du = ddu/dX * dX/dx
  MatMatMult(1.0, ddudX, dXdx, grad_du);

  // Deformation Gradient : F = I + grad_u
  const double F[3][3] = {
      {stored_grad_u[0][0] + 1, stored_grad_u[0][1],     stored_grad_u[0][2]    },
      {stored_grad_u[1][0],     stored_grad_u[1][1] + 1, stored_grad_u[1][2]    },
      {stored_grad_u[2][0],     stored_grad_u[2][1],     stored_grad_u[2][2] + 1}
  };

  // E - Green-Lagrange strain tensor
  GreenLagrangeStrain(stored_grad_u, E_sym);

  // dE - Green-Lagrange strain tensor increment
  GreenLagrangeStrain_fwd(grad_du, F, dE_sym);

  // Compute dS with ADOL-C
  auto Ep = new double[6];
  for (int i=0; i<6; i++) Ep[i] = E_sym[i];
  double hessPsi[6][6] = {{0.}};
  ComputeHessianPsi(hessPsi, Ep, lambda, mu, data);
  for (int i=0; i<6; i++) {
    dS_sym[i] = 0.;
    for (int j=0; j<6; j++) dS_sym[i] += hessPsi[i][j] * dE_sym[j];
  }
  SymmetricMatUnpack(dS_sym, dS);

  // Second Piola-Kirchhoff (S)
  SymmetricMatUnpack(stored_S_sym, S);

  // df1 = dP = dP/dF:dF = dF*S + F*dS; note dF = grad_du
  MatMatMultPlusMatMatMult(grad_du, S, F, dS, df1);

  // Print
  std::cout.precision(12);
  std::cout.setf(std::ios::fixed);

  std::cout << "\n  E:" << std::endl << std::endl;
  for (int i=0; i<6; i++) std::cout << "\t" << E_sym[i] << std::endl;
  std::cout << std::endl;

  std::cout << "\n  dE:" << std::endl << std::endl;
  for (int i=0; i<6; i++) std::cout << "\t" << dE_sym[i] << std::endl;
  std::cout << std::endl;

  std::cout << "\n Strain energy = " << StrainEnergy(E_sym, lambda, mu) << std::endl;

  std::cout << "\n  S (ADOL-C gradient):" << std::endl << std::endl;
  for (int i=0; i<6; i++) std::cout << "\t" << stored_S_sym[i] << std::endl;
  std::cout << std::endl;

  std::cout << "\n dS (ADOL-C hessian) =" << std::endl << std::endl;
  for (int i=0; i<6; i++) std::cout << "\t" << dS_sym[i] << std::endl;
  std::cout << std::endl;
}
