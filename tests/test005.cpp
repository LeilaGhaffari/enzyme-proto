// Vector version of test002
#include <stdio.h>

void foo(double *x_in, double *x_out) {
    x_out[0] = x_in[0] * x_in[0] + 5 * x_in[0] + 3.;
}


int enzyme_dup;
int enzyme_out;
int enzyme_const;

typedef void (*f_ptr)(double *, double *);

extern void __enzyme_autodiff(f_ptr,
    int, double *, double *,
    int, double *, double *);

int main() {

    double x[3] = {1.4, 1.8, 2.};
    double d_x[3] = {0.0, 0.0, 0.0};
    double y[3];
    double d_y[3] = {1.0, 1.0, 1.0};

    printf("    x      f(x)      df/dx\n");
    for (int i = 0; i<3; i++) {
        __enzyme_autodiff(foo,
            enzyme_dup, &x[i], &d_x[i],
            enzyme_dup, &y[i], &d_y[i]);

        printf("%f %f %f\n", x[i], y[i], d_x[i]);
    }
}
