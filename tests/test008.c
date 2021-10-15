// The C version of libCEED-ad/enzyme-velocity-gradient.ipynb
#include <stdio.h>

void foo(double *u, double *x, double mu) {
    for (int i = 0; i<3; i++) u[i] = mu * (x[i] * x[i]);
}

int enzyme_dup;
int enzyme_out;
int enzyme_const;

typedef void (*f_ptr)(double *, double *, double);

extern void __enzyme_autodiff(f_ptr,
    int, double *, double *,
    int, double *, double *,
    int, double);

int main() {

    double x[3] = { 0.653484990079922,
                    0.744572635699587,
                    0.23483535966444968
                  };
    double u[3] = {0., 0., 0.};
    double dx[3][3] = { {0., 0., 0.},
                        {0., 0., 0.},
                        {0., 0., 0.}
                      };
    double mu = 2;
    for (int i=0; i<3; i++) {
        double du[3]  = {0., 0., 0.}; du[i] = 1;
        __enzyme_autodiff(foo,
        enzyme_dup, u, du,
        enzyme_dup, x, dx[i],
        enzyme_const, mu);
    }

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

clang test008.c -Xclang -load -Xclang /home/linuxbrew/.linuxbrew/Cellar/enzyme/HEAD-6e45ead/lib/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops

Output:

dx =
	2.613940 	0.000000 	0.000000 
	0.000000 	2.978291 	0.000000 
	0.000000 	0.000000 	0.939341 


u(x) =
	u(0.653485) = 0.854085
	u(0.744573) = 1.108777
	u(0.234835) = 0.110295

*/