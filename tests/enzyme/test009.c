// The C version of libCEED-ad/enzyme-velocity-gradient-more.ipynb
#include <stdio.h>

void foo(double *u, double *x) {
    for (int i = 0; i<3; i++) 
      u[i] = x[i] * x[i] + 2 * x[0] - x[1] + 5 * x[2];
}

int enzyme_dup;
int enzyme_out;
int enzyme_const;

typedef void (*f_ptr)(double *, double *);

extern void __enzyme_autodiff(f_ptr,
                              int, double *, double *,
                              int, double *, double *);

void getVelocityGradient(double *x, double dx[3][3]) {
    double u[3] = {0., 0., 0.};
    for (int i=0; i<3; i++) {
        double du[3]  = {0., 0., 0.}; du[i] = 1;
        __enzyme_autodiff(foo,
                          enzyme_dup, u, du,
                          enzyme_dup, x, dx[i]);
    }
}

int main() {

    double x[3] = { 
                    0.7428839934549221,
                    0.5048739222425671,
                    0.09766249467142774
                  };

    double dx[3][3] = {{0.}};

    getVelocityGradient(x, dx);

    printf("\n\ndx =\n");
    for (int i=0; i<3; i++) printf("\t%.6lf ", dx[0][i]);
    printf("\n");
    for (int i=0; i<3; i++) printf("\t%.6lf ", dx[1][i]);
    printf("\n"); 
    for (int i=0; i<3; i++) printf("\t%.6lf ", dx[2][i]);
    printf("\n\n");

    return 0;
}

/*

Output:

dx =
	3.485768 	-1.000000 	5.000000 
	2.000000 	0.009748 	5.000000 
	2.000000 	-1.000000 	5.195325 

*/