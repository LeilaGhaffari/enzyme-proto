// constant_fwd.c
// Forward mode
// This function demonstrates the derivative of a constant function.

#include <stdio.h>
extern double __enzyme_fwddiff(void *);
double Constant() { return 2.; }
double dConstant() { return __enzyme_fwddiff((void *) Constant); }
int main() {
    printf(" Constant  = %f \n dConstant = %f \n", Constant(), dConstant());
}
