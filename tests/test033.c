// Use of Enzyme's custom function

#include <stdio.h>
#include <math.h>

double __enzyme_autodiff(void*, ...);

double log1p_like_function(double x) {
  double sum = 0;
  double y   = x / (2. + x);
  double y2  = y * y;
  sum += y;
  for (int i = 0; i < 5; i++) {
    y *= y2;
    sum += y / (2 * i + 3);
  }
  return 2 * sum;
}

double test(double x) {
  return log1p_like_function(x);
}

void* __enzyme_function_like[2] = {(void*)log1p_like_function, "log1p"}; 

int main(int argc, char** argv) {

  double x = 1.45;

  double f  = log1p_like_function(x);
  double df = __enzyme_autodiff(test, x);

  printf("f = %f \ndf = %f \n", f, df);

  return 0;
}

/*
f = 0.896086 
df = 0.408163 
*/
