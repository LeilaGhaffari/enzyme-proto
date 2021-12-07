// Mock Qdot in libCEED (derivative of a vector wrt a scalar)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

int  __enzyme_augmentsize(void *, ...);
void __enzyme_augmentfwd(void *, ...);
void __enzyme_reverse(void *, ...);
int enzyme_dup, enzyme_allocated, enzyme_tape, enzyme_const, enzyme_nofree;

void square(const double *x, double *y, double *A) {
    double a = A[0];
    y[0] = a*a*a + a*a*x[0] + a/2.;
    y[1] = a*a*a + a*a*x[1] + a/2.;
    y[2] = a*a*a + a*a*x[2] + a/2.;
    y[3] = a*a*a + a*a*x[3] + a/2.;
    y[4] = a*a*a + a*a*x[4] + a/2.;
}

int getTapeSize() {
  return __enzyme_augmentsize((void *)square, enzyme_const, enzyme_dup, enzyme_dup);
}

void fwd_a(const double *x, double *y, double *a, void *tape) {  
    __enzyme_augmentfwd((void*)square, 
                        enzyme_allocated, sizeof(tape[0]), 
                        enzyme_tape, tape, 
                        enzyme_const, x,
                        y, (double *)NULL,
                        a, (double *)NULL
                       );
}

void rev_a(const double *x, double *y_, double *a_, void *tape, bool free) {
    if (free)
      __enzyme_reverse((void*)square, 
                       enzyme_allocated, sizeof(tape[0]), 
                       enzyme_tape, tape, 
                       enzyme_const, x, 
                       (double *)NULL, y_,
                       (double *)NULL, a_
                      );
    else
      __enzyme_reverse((void*)square, 
                       enzyme_allocated, sizeof(tape[0]), 
                       enzyme_tape, tape, enzyme_nofree, 
                       enzyme_const, x, 
                       (double *)NULL, y_,
                       (double *)NULL, a_
                      );                      
}

int main() {
    // Declarations
    const double x[] = {.5, 2.5, 5., 8., 10.};
    double y[5];
    void *tape;

    // Get tape size
    int tape_size = getTapeSize();

    // Allocate memory for tape
    tape = malloc(tape_size);

    // Forward mode
    double a[1] = {5.};
    fwd_a(x, y, a, tape);

    // Reverse mode
    double a_[5] = {0.};
    for (int i=0; i<5; i++) {
        double y_[5] = {0.}; y_[i] = 1.;
        rev_a(x, y_, &a_[i], tape, false);
    }
    
    // Print output
    printf("\na,  x,  y,  a_\n");
    for (int i=0; i<5; i++) printf("%g, %g, %g: %g\n", a[i], x[i], y[i], a_[i]);
    
    // Free tape
    {
        double y_[5] = {0.}; y_[0] = 1.;
        rev_a(x, y_, a_, tape, true);
        free(tape);
    }

    return 0;
}

/*
clang test017.c -Xclang -load -Xclang /home/leila/Enzyme/enzyme/build12DHB/Enzyme/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops

Output:

a  x    y    a_
5, 0.5, 140: 80.5
5, 2.5, 190: 100.5
5, 5, 252.5: 125.5
5, 8, 327.5: 155.5
5, 10, 377.5: 175.5

*/
