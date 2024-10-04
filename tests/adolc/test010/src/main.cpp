// Standard Neo-Hookean formulation in initial configuration with ADOL-C

#include "../include/qfunctions.hpp"

int main() {
  // Residual evaluation
  double dXdx[3][3] = {
      {0.555884, 0.660651, 0.233178},
      {0.871235, 0.969807, 0.904515},
      {0.287474, 0.13347, 0.843519}
  };
  double dudX[3][3] = {
      {0.510437, 0.707667, 0.782745},
      {0.921497, 0.477494, 0.655606},
      {0.425158, 0.889053, 0.166861}
  };

  double stored_grad_u[3][3], stored_S_sym[6], f1[3][3];
  f1_NeoHookeanInitial_AD_ADOLC(dXdx, dudX, stored_grad_u, stored_S_sym, f1);

  // Jacobian evaluation
  double ddudX[3][3] = {
    {0.1425560,  0.115120,  0.551640},
    {0.0591922,  0.123535,  0.166572},
    {0.1617210,  0.478828,  0.646217}
  };

  double df1[3][3];
  HessianData *data = new HessianData;
  initializeHessianData(data);
  df1_NeoHookeanInitial_AD_ADOLC(data, dXdx, ddudX, stored_grad_u, stored_S_sym, df1);

  // Free allocated memory
  freeHessianData(data);

  return 0;
}

/*
  E:

        3.165062699399
        4.595009711411
        2.045993834836
        3.272300872638
        2.625201021139
        3.716182947107


  dE:

        1.433414579285
        1.447099855442
        2.413871039923
        1.913840224073
        1.943720418915
        1.419804212546


 Strain energy = 11.332917346807

  S (ADOL-C gradient):

        4.802659104984
        5.097673027985
        9.169178355792
        -3.758231888851
        -2.040744925345
        -1.462891327345


 dS (ADOL-C hessian) =

        7.307220006812
        5.464120433286
        6.671993446558
        -2.749064108456
        -4.069974200124
        -2.593362694953
*/
