// Replicate test017 with autodiff instead of split mode

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

void __enzyme_autodiff(void *, ...);
int enzyme_const;

void foo(const double *x, double *y, double *A) {
    double a = A[0];
    y[0] = a*a*a + a*a*x[0] + a/2.;
    y[1] = a*a*a + a*a*x[1] + a/2.;
    y[2] = a*a*a + a*a*x[2] + a/2.;
    y[3] = a*a*a + a*a*x[3] + a/2.;
    y[4] = a*a*a + a*a*x[4] + a/2.;
}

void foo_a(const double *x, double *y, double *y_,  double *a, double *a_) {  
    __enzyme_autodiff((void*)foo,  
                      enzyme_const, x,
                      y, y_,
                      a, a_
                     );
}

int main() {
    // Declarations
    const double x[] = {.5, 2.5, 5., 8., 10.};
    double y[5];
    double a[1] = {5.};
    double a_[5] = {0.};

    // AD 
    for (int i=0; i<5; i++) {
        double y_[5] = {0.}; y_[i] = 1.;
        foo_a(x, y, y_, a, &a_[i]);
    }
    
    // Print output
    printf("\na,  x,  y,  a_\n");
    for (int i=0; i<5; i++) printf("%g, %g, %g: %g\n", a[i], x[i], y[i], a_[i]);

    return 0;
}

/*
clang test018.c -Xclang -load -Xclang /home/leila/Enzyme/enzyme/build12DHB/Enzyme/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops

Output:

a  x    y    a_
5, 0.5, 140: 80.5
5, 2.5, 190: 100.5
5, 5, 252.5: 125.5
5, 8, 327.5: 155.5
5, 10, 377.5: 175.5

*/
