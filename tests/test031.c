// Test for forward mode vector with for loop

#include <stdio.h>

#define VECSIZE 2

double __enzyme_fwddiff(void*, ...);
 
void fun(double x[VECSIZE], double *result) {

  double w1 = x[0];
  double w2 = x[1];
  double w3 = w1 * w2;
  double w4 = 1.0 / w1;
  double w5 = w3 + w4;
  
  *result = w5;
}

int main(int argc, char** argv) {

  for (int j = 1; j < 5; j++) {
    double x[VECSIZE]; 
    x[0] = (double)j; x[1] = x[0]*x[0]; 
    double res;
    double dres[VECSIZE];

    printf("x: %f, y: %f\n", x[0], x[1]);

    for (int i = 0; i < VECSIZE; i++) {
      double dx[VECSIZE] = {0.}; dx[i] = 1.;
      __enzyme_fwddiff(fun, x, dx, &res, &dres[i]);
      printf("dres[%d] %f\n", i, dres[i]);
    }
    printf("res %f\n", res);
    
    printf("\n");
  }

  return 0;
}

/*

x: 1.000000, y: 1.000000
dres[0] 0.000000
dres[1] 1.000000
res 2.000000

x: 2.000000, y: 4.000000
dres[0] 3.750000
dres[1] 2.000000
res 8.500000

x: 3.000000, y: 9.000000
dres[0] 8.888889
dres[1] 3.000000
res 27.333333

x: 4.000000, y: 16.000000
dres[0] 15.937500
dres[1] 4.000000
res 64.250000

*/
