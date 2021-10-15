// The C version of libCEED-ad/enzyme-velocity-gradient-more.ipynb
#include <stdio.h>

void foo(double *u, double *x) {
    for (int i = 0; i<3; i++) u[i] = x[i] * x[i] + 2 * x[0] - x[1] + 5 * x[2];
}

int enzyme_dup;
int enzyme_out;
int enzyme_const;

typedef void (*f_ptr)(double *, double *);

extern void __enzyme_autodiff(f_ptr,
    int, double *, double *,
    int, double *, double *);

int main() {

    double x[3] = { 0.7428839934549221,
                    0.5048739222425671,
                    0.09766249467142774
                  };
    double u[3] = {0., 0., 0.};
    double dx[3][3] = { {0., 0., 0.},
                        {0., 0., 0.},
                        {0., 0., 0.}
                      };
    for (int i=0; i<3; i++) {
        double du[3]  = {0., 0., 0.}; du[i] = 1;
        __enzyme_autodiff(foo,
        enzyme_dup, u, du,
        enzyme_dup, x, dx[i]);
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

clang test009.c -Xclang -load -Xclang /home/linuxbrew/.linuxbrew/Cellar/enzyme/HEAD-6e45ead/lib/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops

Output:

dx =
	3.485768 	-1.000000 	5.000000 
	2.000000 	0.009748 	5.000000 
	2.000000 	-1.000000 	5.195325 


u(x) =
	u(0.742884) = 2.021083
	u(0.504874) = 1.724104
	u(0.097662) = 1.478745

*/