// Mocking  the function computeS() in 
// "libCEED/examples/solids/qfunctions/finite-strain-neo-hookean-initial-1.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double log1p_series_shifted(double x) {
    
    double left = sqrt(2.)/2 - 1;
    double right = sqrt(2.) - 1;
    double sum = 0;

    if (x < left) {
        sum -= log(2.) / 2;
        x = 1 + 2 * x;
    } else if (right < x) {
        sum += log(2.) / 2;
        x = (x - 1) / 2;
    }
    
    double y = x / (2. + x);
    double y2 = y*y;
    sum += y;
    y *= y2;
    sum += y / 3;
    y *= y2;
    sum += y / 5;
    y *= y2;
    sum += y / 7;

    return 2 * sum;
};

double computeDetCM1(double E2work[6]) {
  return E2work[0]*(E2work[1]*E2work[2]-E2work[3]*E2work[3]) +
         E2work[5]*(E2work[4]*E2work[3]-E2work[5]*E2work[2]) +
         E2work[4]*(E2work[5]*E2work[3]-E2work[4]*E2work[1]) +
         E2work[0] + E2work[1] + E2work[2] +
         E2work[0]*E2work[1] + E2work[0]*E2work[2] +
         E2work[1]*E2work[2] - E2work[5]*E2work[5] -
         E2work[4]*E2work[4] - E2work[3]*E2work[3];
};

int computeMatinvSym(double A[3][3], double detA, double Ainv[6]) {
    // Compute A^(-1) : A-Inverse
    double B[6] = {A[1][1]*A[2][2] - A[1][2]*A[2][1], 
                   A[0][0]*A[2][2] - A[0][2]*A[2][0], 
                   A[0][0]*A[1][1] - A[0][1]*A[1][0], 
                   A[0][2]*A[1][0] - A[0][0]*A[1][2], 
                   A[0][1]*A[1][2] - A[0][2]*A[1][1], 
                   A[0][2]*A[2][1] - A[0][1]*A[2][2]  
                  };
                      
    for ( int m = 0; m < 6; m++) Ainv[m] = B[m] / (detA);
    
    return 0;
};

int computeS(double lambda, double mu, double E2work[6], double Swork[6]) {
    
    double E2[3][3] = {{E2work[0], E2work[5], E2work[4]},
                       {E2work[5], E2work[1], E2work[3]},
                       {E2work[4], E2work[3], E2work[2]}
                      };
   
    // C : right Cauchy-Green tensor
    // C = I + 2E
    double C[3][3] = {{1 + E2[0][0], E2[0][1], E2[0][2]},
                      {E2[0][1], 1 + E2[1][1], E2[1][2]},
                      {E2[0][2], E2[1][2], 1 + E2[2][2]}
                     };
   

    // Compute C^(-1) : C-Inverse
    double Cinvwork[6];
    double detCm1 = computeDetCM1(E2work);
    computeMatinvSym(C, detCm1+1, Cinvwork);

   
    double C_inv[3][3] = {{Cinvwork[0], Cinvwork[5], Cinvwork[4]},
                          {Cinvwork[5], Cinvwork[1], Cinvwork[3]},
                          {Cinvwork[4], Cinvwork[3], Cinvwork[2]}
                         };
       

    // Compute the Second Piola-Kirchhoff (S)
    int indj[6] = {0, 1, 2, 1, 0, 0}, indk[6] = {0, 1, 2, 2, 2, 1};
    double logJ = log1p_series_shifted(detCm1) / 2.;
    for ( int m = 0; m < 6; m++) {
        Swork[m] = lambda*logJ*Cinvwork[m];
        for ( int n = 0; n < 3; n++)
            Swork[m] += mu*C_inv[indj[m]][n]*E2[n][indk[m]];
    }
    
    return 0;
};

int main() {
    
    double E = .3;
    double nu = .3;
    double TwoMu = E / (1 + nu);
    double mu = TwoMu / 2;
    double Kbulk = E / (3*(1 - 2*nu)); // Bulk Modulus
    double lambda = (3*Kbulk - TwoMu) / 3;

    double E2work[6] = {0., 0., 0., 0., 0., 0.};
    E2work[0] = 0.5895232828911128;
    E2work[1] = 0.2362491738162759;
    E2work[2] = 0.9793730522395296;
    E2work[3] = 0.2190993957421843;
    E2work[4] = 0.0126503210747925;
    E2work[5] = 0.6570956167695403;

    double Swork[6];
    computeS(lambda, mu, E2work, Swork);

    double deltaEwork[6] = {0., 0., 0., 0., 0., 0.};
    deltaEwork[0] = 0.9681576729097205;
    deltaEwork[1] = 0.7994338113484318;
    deltaEwork[2] = 0.2755183472001872;
    deltaEwork[3] = 0.6500440500146469;
    deltaEwork[4] = 0.0593948875992271;
    deltaEwork[5] = 0.6002528007029311;

    double deltaE[3][3] = {{deltaEwork[0], deltaEwork[5], deltaEwork[4]},
                           {deltaEwork[5], deltaEwork[1], deltaEwork[3]},
                           {deltaEwork[4], deltaEwork[3], deltaEwork[2]}
                          };

    double E2[3][3] = {{E2work[0], E2work[5], E2work[4]},
                       {E2work[5], E2work[1], E2work[3]},
                       {E2work[4], E2work[3], E2work[2]}
                      };   

    double C[3][3] = {{1 + E2[0][0], E2[0][1], E2[0][2]},
                              {E2[0][1], 1 + E2[1][1], E2[1][2]},
                              {E2[0][2], E2[1][2], 1 + E2[2][2]}
                             };             
    double Cinvwork[6];
    double detCm1 = computeDetCM1(E2work);
    computeMatinvSym(C, detCm1+1, Cinvwork);
    double logJ = log1p_series_shifted(detCm1) / 2.;

    double C_inv[3][3] = {{Cinvwork[0], Cinvwork[5], Cinvwork[4]},
                          {Cinvwork[5], Cinvwork[1], Cinvwork[3]},
                          {Cinvwork[4], Cinvwork[3], Cinvwork[2]}
                         };

    double Cinv_contract_E = 0;
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        Cinv_contract_E += C_inv[j][k]*deltaE[j][k];

    // -- deltaE*C_inv
    double deltaECinv[3][3];
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) {
        deltaECinv[j][k] = 0;
        for (int m = 0; m < 3; m++)
          deltaECinv[j][k] += deltaE[j][m]*C_inv[m][k];
      }
    // -- intermediate deltaS = C_inv*deltaE*C_inv
    double deltaS[3][3];
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) {
        deltaS[j][k] = 0;
        for (int m = 0; m < 3; m++)
          deltaS[j][k] += C_inv[j][m]*deltaECinv[m][k];
      }
    // -- deltaS = lambda(C_inv:deltaE)C_inv - 2(lambda*log(J)-mu)*(intermediate)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        deltaS[j][k] = lambda*Cinv_contract_E*C_inv[j][k] -
                       2.*(lambda*logJ-mu)*deltaS[j][k];


    printf("\n\nSwork =\n");
    for (int i=0; i<6; i++) printf("\t   %.6lf", Swork[i]);

    printf("\n\ndeltaS[] =\n");
    for (int i=0; i<3; i++) printf("%.6lf ", deltaS[0][i]);
    printf("\n");
    for (int i=0; i<3; i++) printf("%.6lf ", deltaS[1][i]);
    printf("\n");
    for (int i=0; i<3; i++) printf("%.6lf ", deltaS[2][i]);
    printf("\n");

    return 0;
}

/*

Output:

Swork =
	   0.098041	   0.092640	   0.104300	   0.002458	   -0.000928	   0.009383

deltaS[] =
0.169671 -0.085117 0.003474 
-0.085117 0.219554 -0.010067 
0.003474 -0.010067 0.099308 

*/
