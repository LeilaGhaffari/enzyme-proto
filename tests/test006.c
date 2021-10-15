// A simple test for split mode

#include <stdio.h>
#include <stdlib.h>

double tester(double* array) {
    return array[0] * array[0];
}

int  __enzyme_augmentsize(void*, ...);
void __enzyme_augmentfwd(void*, ...);
void __enzyme_reverse(void*, ...);

int enzyme_dup;
int enzyme_tape;
int enzyme_allocated;

void test_derivative(double* x, double* dx) {
    int size = __enzyme_augmentsize((void*)tester, enzyme_dup);
    void* data = malloc(size);
    __enzyme_augmentfwd((void*)tester, enzyme_allocated, size, enzyme_tape, data, x, dx);

    __enzyme_reverse((void*)tester, enzyme_allocated, size, enzyme_tape, data, x, dx);
    free(data);
}

int main() {

    double in[1] = {3.14};
    double grad_in[1] = {0.00};

    test_derivative(in, grad_in);

    printf("grad: %f\n", grad_in[0]);
}

/*

clang test006.c -Xclang -load -Xclang /home/linuxbrew/.linuxbrew/Cellar/enzyme/HEAD-6e45ead/lib/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops

*/
