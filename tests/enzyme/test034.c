#include <stdio.h>

double __enzyme_autodiff(void*, double);

__attribute__((optnone))
void square_(const double* src, double* dest) {
    *dest = *src * *src;
}

int augment = 0;
void* augment_square_(const double* src, const double *d_src, double* dest, double* d_dest) {
    augment++;
    // intentionally incorrect for debugging
    *dest = 7.0;
    *d_dest = 11.0;
    return NULL;
}

int gradient = 0;
void gradient_square_(const double* src, double *d_src, const double* dest, const double* d_dest, void* tape) {
    gradient++;
    // intentionally incorrect for debugging
    *d_src = 13.0;
}

void* __enzyme_register_gradient_square[] = {
    (void*)square_,
    (void*)augment_square_,
    (void*)gradient_square_,
};


double square(double x) {
    double y;
    square_(&x, &y);
    return y;
}

double dsquare(double x) {
    return __enzyme_autodiff((void*)square, x);
}


int main() {
    double res = dsquare(3.0);
    printf("res=%f augment=%d gradient=%d\n", res, augment, gradient);
}
