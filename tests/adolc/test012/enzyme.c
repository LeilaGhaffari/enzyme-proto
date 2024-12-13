#include <stdio.h>

extern double __enzyme_autodiff(void *, ...);
extern int enzyme_const;

double square(double *restrict x, int len) { double sum =0; for (int i=0; i<len; i++) sum += x[i] * x[i]; return sum; }

double dsquare(double *restrict x, double *restrict bx, int len) {
  return __enzyme_autodiff((void *)square, x, bx, len);
}

int main() {
    double x[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double y[5] = {0.0};
    double r = dsquare(x, y, 5);
    printf("square=%f, dsquare=%f\n", square(x, 5), r);
    printf("%f %f %f %f %f\n", y[0], y[1], y[2],y[3],y[4]);
}

/*
clang-18 -Os -fpass-plugin=/home/leila/Enzyme/enzyme/build18g/Enzyme/ClangEnzyme-18.so \
-Xclang -load -Xclang /home/leila/Enzyme/enzyme/build18g/Enzyme/ClangEnzyme-18.so -mllvm \
-enzyme-print -mllvm -enzyme-print-perf enzyme.c
*/