#include <stdio.h>

void foo(double *x_in, double *x_out) {
    x_out[0] = x_in[0] * x_in[0];
}


int enzyme_dup;
int enzyme_out;
int enzyme_const;

typedef void (*f_ptr)(double *, double *);

extern void __enzyme_autodiff(f_ptr,
    int, double *, double *,
    int, double *, double *);

int main() {

    double x = 1.4;
    double d_x = 0.0;
    double y;
    double d_y = 1.0;

    __enzyme_autodiff(foo,
        enzyme_dup, &x, &d_x,
        enzyme_dup, &y, &d_y);

    printf("%f %f\n", x, y);
    printf("%f %f\n", d_x, d_y);

}

/*

clang test2.cpp -Xclang -load -Xclang /home/linuxbrew/.linuxbrew/Cellar/enzyme/0.0.19/lib/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops

./a.out

*/