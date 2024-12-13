#include <iostream>
#include <adolc/adolc.h>

adouble square(const adouble *x, const int len) {
    adouble sum = 0;
    for (int i=0; i<len; i++) sum += x[i] * x[i];
    return sum;
}

void __attribute__((noinline)) dsquare(const double *x, double *xb, const int len) {
    adouble xa[3];
    adouble ya[1];
    double y[1] = {0.};

    int tag = 1;
    trace_on(tag);
    for (int i=0; i<len; i++) xa[i] <<= x[i];
    ya[0] = square(xa, len);
    ya[0] >>= y[0];
    trace_off();

    gradient(tag, len, x, xb);
}

int main() {
    const int len = 3;
    const double x[len] = {1.0, 2.0, 3.0};
    double xb[len] = {0.0};
    dsquare(x, xb, len);
    for (int i=0; i<len; i++) std::cout << xb[i] << " ";
    std::cout<<"\n";
    return 0;
}

// clang++-18 -std=c++11 -emit-llvm -S -Os -I$ADOLC_INCLUDE adolc.cpp -o adolc.ll
// clang++-18 -std=c++11 -I$ADOLC_INCLUDE -L$ADOLC_LIB -o adolc-exec adolc.cpp -ladolc