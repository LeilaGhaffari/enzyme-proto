#include <stdio.h>

double __enzyme_autodiff(void*, ...);

double log1p_like_function(double a) {
  return 2*a;
}

double test(double a) {
  return log1p_like_function(a);
}

void* __enzyme_function_like[2] = {(void*)log1p_like_function, "log1p"}; 

int main(int argc, char** argv) {

  double out = __enzyme_autodiff(test, 2.0); 

  printf("result = %f\n", out);

  return 0;
}

/*
Correct answer: 1/3
Wrong answer:   2
*/
