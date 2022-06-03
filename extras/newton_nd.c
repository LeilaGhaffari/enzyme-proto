#include <stdio.h> 
#include <stdlib.h>
#include <math.h>

int enzyme_const;
int enzyme_dup;

// Define a function pointer?
typedef void (f_ptr)(double*, double*);


extern void __enzyme_autodiff(f_ptr,
                                int, double*, double*,
                                int, double*, double*);

void func(double *F, double *x){

    // Define equation and size
    F[0] = 3*x[0]*x[0] - x[1]*x[1];
    F[1] = 3*x[0]*x[1] - x[0]*x[0]*x[0] - 1;


}

void get_Jacobian(double* x, double J[2][2]){
 // Pre define 
 double F[2] = {0.,0.};

// Calculate Jacobian 
for (int i=0; i<2; i++){
    double dF[2] = {0., 0.}; dF[i] = 1;
    __enzyme_autodiff(func,
    enzyme_dup, F, dF,
    enzyme_dup, x, J[i]);
}

}




int main(){

    double x[2] = {2.,3.};

    double J[2][2] = { {0.,0.},
                       {0.,0.} };

    get_Jacobian(x, J);

    for(int i=0; i<2; i++) printf("\t%f\n", J[0][i]);
    for(int i=0; i<2; i++) printf("\t%f\n", J[1][i]);

return 0;
}