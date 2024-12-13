#include <stdio.h>

double square(const double *restrict x, int len) {
    double sum =0;
    for (int i=0; i<len; i++)
    sum += x[i] * x[i];
    return sum; }

// Tapenade's reverse mode AD
void dsquare(const double *restrict x, double *restrict xb, int len, double squareb) {
    double sum = 0.0;
    double sumb = 0.0;
    sumb = squareb;
    for (int i = len-1; i > -1; --i)
        xb[i] = xb[i] + 2*x[i]*sumb;
}

int main() {
    double x[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double y[5] = {0.0};
    double s_ = 1.;
    dsquare(x, y, 5, s_);
    printf("square = %f\n", square(x, 5));
    printf("%f %f %f %f %f\n", y[0], y[1], y[2],y[3],y[4]);
}
