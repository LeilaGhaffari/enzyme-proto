// A split version of test002
#include <stdio.h>
#include <stdlib.h>

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

template <typename... Args>
void* __enzyme_augmentfwd(Args...);

template <typename... Args>
void __enzyme_reverse(Args...);

int main() {

    double x = 1.4;
    double d_x = 0.0;
    double y;
    double d_y = 1.0;

    __enzyme_augmentfwd((void *)foo, x, d_x, y, d_y);
    //void *tape = (void*) malloc(1 * sizeof(void *));
    void **tape = (void **) malloc(1*sizeof(void*)); 
    __enzyme_reverse((void *)foo, x, d_x, y, d_y, tape);

    printf("%f %f\n", x, y);
    printf("%f %f\n", d_x, d_y);
}

/*

clang test006.cpp -Xclang -load -Xclang /home/linuxbrew/.linuxbrew/Cellar/enzyme/0.0.19/lib/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops; 
./a.out

*/