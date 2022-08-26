// scalar_rev.c
//  This function demonstrates the derivative of a scalar function.

#include <stdio.h>

extern double __enzyme_autodiff(void *, ...);
int enzyme_const;

void V_cylinder(double *v, double *r, double *h, const double pi) {
    *v = pi * (*r) * (*r) * (*h);
}

void dV_cylinder(double *v, double *r, double *dr, double *h, double *dh, const double pi) {
    double dv = 1;
    __enzyme_autodiff((void *) V_cylinder,
                      v, &dv,
                      r, dr,
                      h, dh,
                      enzyme_const, pi);
}

int main() {
    const double pi = 3.141593;
    double r = 3, h = 2, volume_cylinder = 0, surface_lateral = 0, surface_cross_sectional = 0;
    dV_cylinder(&volume_cylinder, &r, &surface_lateral, &h, &surface_cross_sectional, pi);
    printf("Cylinder with r=%f and h=%f \n", r, h);
    printf("----------------------------------------\n");
    printf("Volume                       = %f \n", volume_cylinder);
    printf("Lateral Surface Area         = %f \n", surface_lateral);
    printf("Cross-sectional Surface Area = %f \n", surface_cross_sectional);

    return 0;
}
