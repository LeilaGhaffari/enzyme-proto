// Vector version of test006

#include <stdio.h>
#include <stdlib.h>

void tester(double* array) {
    for (int i = 0; i<3; i++) array[i] *= array[i];
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

    int size = 3;
    double in[3] = {3, 4, 5};
    double grad_in[3] = {1.00, 1.00, 1.00};

    test_derivative(in, grad_in);
    for (int i = 0; i<size; i++) {
        printf("%f %f\n", in[i], grad_in[i]);
    }
}

/*

clang test007.c -Xclang -load -Xclang /home/linuxbrew/.linuxbrew/Cellar/enzyme/HEAD-6e45ead/lib/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops

Output:
9.000000 6.000000
16.000000 8.000000
25.000000 10.000000

*/
