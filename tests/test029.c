// First test for forward mode vector

#include <stdio.h>

typedef struct {
  double dx, dy;
} Gradients;

extern Gradients __enzyme_fwddiff(double (*)(double, double), ...);
extern int enzyme_width;

double fun(double x1, double x2) {
  double w1 = x1;
  double w2 = x2;
  double w3 = w1 * w2;
  double w4 = 1.0 / w1;
  double w5 = w3 + w4;
  return w5;
}

Gradients dfun(double x, double y) {
  return __enzyme_fwddiff(fun, enzyme_width, 2, x, 1.0, 0.0, y, 0.0, 1.0);
}

int main() {
  double xres[] = {0.0, 3.75, 8.8888888888888893, 15.9375};
  double yres[] = {1.0, 2.0, 3.0, 4.0};

  for (int i = 1; i < 5; i++) {
    double x = (double)i;
    Gradients df = dfun(x, x * x);
    printf("dx = %f ? %f \n",   df.dx, xres[i - 1]);
    printf("dy = %f ? %f \n\n", df.dy, yres[i - 1]);
  }
}
