// Separate reverse and forward calls in test011

#include <stdio.h>
#include <stdlib.h>

void foo(double *u, double *x, double mu) {
    for (int i = 0; i<3; i++) u[i] = mu * (x[i] * x[i]);
}

int  __enzyme_augmentsize(void *, ...);
void __enzyme_augmentfwd(void *, ...);
void __enzyme_reverse(void *, ...);

int enzyme_dup;
int enzyme_tape;
int enzyme_const;
int enzyme_allocated;

void grad_foo_fwd(double *u, double *du, double *x, double *dx, double mu, void *data) {
    int size = __enzyme_augmentsize((void *)foo, enzyme_dup, enzyme_dup, enzyme_const);
    __enzyme_augmentfwd((void *)foo, enzyme_allocated, size, enzyme_tape, data, u, du, x, dx, mu);
}

void grad_foo_rev(double *u, double *du, double *x, double *dx, double mu, void *data) {
    int size = __enzyme_augmentsize((void *)foo, enzyme_dup, enzyme_dup, enzyme_const);
    __enzyme_reverse((void *)foo, enzyme_allocated, size, enzyme_tape, data, u, du, x, dx, mu);
}

int main() {

    double x[3] = { 0.653484990079922,
                    0.744572635699587,
                    0.23483535966444968
                  };
    double u[3]; // No need for initialization
    double dx[3][3] = { {0., 0., 0.},
                        {0., 0., 0.},
                        {0., 0., 0.}
                      };
    double mu = 2; 

    int size = __enzyme_augmentsize((void *)foo, enzyme_dup, enzyme_dup, enzyme_const);

    void *tape[3];
    for (int i=0; i<3; i++) tape[i] = malloc(size);  

    for (int i=0; i<3; i++) {
        double du[3]  = {0., 0., 0.}; du[i] = 1;
        grad_foo_fwd(u, du, x, dx[i], mu, tape[i]);
    }

    for (int i=0; i<3; i++) {
        double du[3]  = {0., 0., 0.}; du[i] = 1;
        grad_foo_rev(u, du, x, dx[i], mu, tape[i]);
    }
     
    for (int i=0; i<3; i++) free(tape[i]);

    printf("\n\ndx =\n");
    for (int i=0; i<3; i++) printf("\t%.6lf ", dx[0][i]);
    printf("\n");
    for (int i=0; i<3; i++) printf("\t%.6lf ", dx[1][i]);
    printf("\n"); 
    for (int i=0; i<3; i++) printf("\t%.6lf ", dx[2][i]);
    printf("\n\n");
    printf("\nu(x) =\n");
    for (int i=0; i<3; i++) printf("\tu(%f) = %f\n", x[i], u[i]);
    printf("\n\n");

    return 0;
}

/*
clang test012.c -S -emit-llvm -o input.ll -O2 -fno-vectorize -fno-slp-vectorize -fno-unroll-loops

opt input.ll -load=/home/linuxbrew/.linuxbrew/Cellar/enzyme/HEAD-6e45ead/lib/LLVMEnzyme-12.so -enzyme -o output.ll -S

opt output.ll -O2 -o output_opt.ll -S

clang output_opt.ll -o a_opt.exe; ./a_opt.exe

// OR

clang output.ll -o a.exe; ./a.exe 

//// A shorter path:

clang test012.c -Xclang -load -Xclang /home/linuxbrew/.linuxbrew/Cellar/enzyme/HEAD-6e45ead/lib/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops

//// Output:

dx =
	2.613940 	0.000000 	0.000000 
	0.000000 	2.978291 	0.000000 
	0.000000 	0.000000 	0.939341 


u(x) =
	u(0.653485) = 0.854085
	u(0.744573) = 1.108777
	u(0.234835) = 0.110295

*/
