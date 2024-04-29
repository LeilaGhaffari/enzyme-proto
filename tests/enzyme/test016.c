// Derivative of a scalar with respect to a scalar
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

int  __enzyme_augmentsize(void *, ...);
void __enzyme_augmentfwd(void *, ...);
void __enzyme_reverse(void *, ...);
int enzyme_dup, enzyme_allocated, enzyme_tape, enzyme_const;

void square(double *x, double *y, double *a) {
    double time = a[0];
    y[0] = time*time*time + time*time*x[0] + time/2.;
}

void fwd_a(double *x, double *y, double *a, void *tape) {  
    __enzyme_augmentfwd((void*)square, 
                        enzyme_allocated, sizeof(tape[0]), 
                        enzyme_tape, tape, 
                        enzyme_const, x,
                        y, (double *)NULL,
                        a, (double *)NULL
                       );
}

void rev_a(double *a_, double *x, void *tape) {
    double y_ = 1.;
    __enzyme_reverse((void*)square, 
                     enzyme_allocated, sizeof(tape[0]), 
                     enzyme_tape, tape, 
                     enzyme_const, x, 
                     (double *)NULL, &y_,
                     (double *)NULL, a_
                    );
}

int main() {
    // Declarations
    double a[] = {1, 2, 3, 4, 5};
    double x[] = {.5, 2.5, 5., 8., 10.};
    double y[5];
    void *tape[5];

    // Get tape size
    int tape_size = __enzyme_augmentsize((void*)square, enzyme_const, enzyme_dup, enzyme_dup);

    // Allocate memory for tape
    for (int i=0; i<5; i++) tape[i] = malloc(tape_size);

    // Forward mode
    for (int i=0; i<5; i++) {
        fwd_a(&x[i], &y[i], &a[i], tape[i]);
    }

    printf("a,  x,  y,  a_\n\n");
    // Reverse mode
    for (int i=0; i<5; i++) {
        double a_ = 0.;
        rev_a(&a_, &x[i], tape[i]);
        printf("%g, %g, %g: %g\n", a[i], x[i], y[i], a_);
    }

    // Free tape
    for (int i=0; i<5; i++) free(tape[i]);

    return 0;
}

/*

a,  x,  y,  a_

1, 0.5, 2: 4.5
2, 2.5, 19: 22.5
3, 5, 73.5: 57.5
4, 8, 194: 112.5
5, 10, 377.5: 175.5

*/
