// Expand test013 to get a tape[][]

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

void foo(double *u, double *x, double mu) {
    for (int i = 0; i<3; i++) u[i] = mu * (x[i] * x[i]);
}

int  __enzyme_augmentsize(void *, ...);
void __enzyme_augmentfwd(void *, ...);
void __enzyme_reverse(void *, ...);

int enzyme_dup;
int enzyme_tape;
int enzyme_const;
int enzyme_allocated;
int enzyme_nofree;

void grad_foo_fwd(double *u, double *du, double *x, double *dx, double mu, void *tape) {
    __enzyme_augmentfwd((void *)foo, enzyme_allocated, sizeof(tape[0]), enzyme_tape, tape, enzyme_nofree, u, du, x, dx, enzyme_const, mu);
}

void grad_foo_rev(double *u, double *du, double *x, double *dx, double mu, void *tape, bool no_free) {
    if (no_free)
      __enzyme_reverse((void *)foo, enzyme_allocated, sizeof(tape[0]), enzyme_tape, tape, enzyme_nofree, u, du, x, dx, enzyme_const, mu);
    else 
      __enzyme_reverse((void *)foo, enzyme_allocated, sizeof(tape[0]), enzyme_tape, tape, u, du, x, dx, enzyme_const, mu);
}

int main() {

    double x[2][3] = { {0.653484990079922, 0.744572635699587, 0.23483535966444968},
                       {2.*0.653484990079922, 2.*0.744572635699587, 2.*0.23483535966444968}
                     };
    double J[3][3][2] = {{{0.}}};
    //for (int j=0; j<2; j++) for (int i=0; i<3; i++) for (int k=0; k<3; k++) J[i][k][j] = 0.;

    double mu = 2;
    void *tape[2];

    for (int j=0; j<2; j++) {

      int size = __enzyme_augmentsize((void *)foo, enzyme_dup, enzyme_dup, enzyme_const);
      tape[j] = malloc(size); 

      double u[3];
      double dx[3][3] = {{0.}};
  
      grad_foo_fwd(u, (double *)NULL, x[j], (double *)NULL, mu, tape[j]);

      bool no_free = true;
      for (int i=0; i<3; i++) {
          double du[3]  = {0.}; du[i] = 1;
          grad_foo_rev((double *)NULL, du, (double *)NULL, dx[i], mu, tape[j], no_free);
          for (int k=0; k<3; k++) J[i][k][j] = dx[i][k];
      }
       
      
      {   
          // Call without enzyme_nofree to free up the allocated memory
          no_free = false;
          double du = 1;
          double dx[3]  = {0.};
          grad_foo_rev((double *)NULL, &du, (double *)NULL, &dx[0], mu, tape[j], no_free);
      }
      free(tape[j]);


      printf("\n\n");
      for (int i=0; i<3; i++) printf("\t%.6lf ", J[0][i][j]);
      printf("\n");
      for (int i=0; i<3; i++) printf("\t%.6lf ", J[1][i][j]);
      printf("\n"); 
      for (int i=0; i<3; i++) printf("\t%.6lf ", J[2][i][j]);
      printf("\n\n");

      printf("\nu(x) =\n");
      for (int i=0; i<3; i++) printf("\tu(%f) = %f\n", x[j][i], u[i]);
      printf("\n\n");
    }
    return 0;
}

/*

clang test014.c -Xclang -load -Xclang /home/leila/Enzyme/enzyme/build12DHB/Enzyme/ClangEnzyme-12.so -O2 -fno-vectorize -fno-unroll-loops

//// Output:

        2.613940        0.000000        0.000000 
        0.000000        2.978291        0.000000 
        0.000000        0.000000        0.939341 


u(x) =
        u(0.653485) = 0.854085
        u(0.744573) = 1.108777
        u(0.234835) = 0.110295


        5.227880        0.000000        0.000000 
        0.000000        5.956581        0.000000 
        0.000000        0.000000        1.878683 


u(x) =
        u(1.306970) = 3.416341
        u(1.489145) = 4.435107
        u(0.469671) = 0.441181

*/
