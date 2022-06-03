#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern double __enzyme_autodiff(void *, double);

double func(double x){
return x * x * x - 4 ;
}

double dfunc(double x){
// returns the derivative of func
return __enzyme_autodiff((void *) func, x);
}


// Newtons Method using enzyme Autodiff
double newton(int n_max, double target_tol, double guess){
// define variables
double tol = 1;
int n = 0;
double xs[n_max];
xs[0] = guess;

printf("|  n  |  x  |  f(x)  |\n|  %i  |  %f  |  %e  |\n", n, xs[n], func(xs[n]));
// Newton's 
while (tol > target_tol && n < n_max)
{
    xs[n + 1] = xs[n] - func(xs[n])/dfunc(xs[n]);
    tol = fabs(func(xs[n + 1]));
    n ++;
    printf("|  %i  |  %f  |  %e  |\n", n, xs[n], func(xs[n]));
}
return xs[n];
}

int main() {
double root = newton(100, 1E-10, 1.5);
}