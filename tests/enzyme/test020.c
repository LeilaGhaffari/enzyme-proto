// Compute div(U)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

int  __enzyme_augmentsize(void *, ...);
void __enzyme_augmentfwd(void *, ...);
void __enzyme_reverse(void *, ...);
void __enzyme_autodiff(void *, ...);
int enzyme_dup, enzyme_allocated, enzyme_tape, enzyme_const, enzyme_nofree;

int ExactSolution(double q[], double *time, double *X, double *Y, double *Z) {

  // -- Time
  double t = time[0];
  // -- Coordinates
  double x = X[0];
  double y = Y[0];
  double z = Z[0];

  // Exact solutions
  q[0] = 10.*(t*t+t+1.) + (1.*x) + (6. *y) + (11.*z);
  q[1] = 20.*(t*t+t+1.) + (2.*x) + (7. *y) + (12.*z);
  q[2] = 30.*(t*t+t+1.) + (3.*x) + (8. *y) + (13.*z);
  q[3] = 40.*(t*t+t+1.) + (4.*x) + (9. *y) + (14.*z);
  q[4] = 50.*(t*t+t+1.) + (5.*x) + (10.*y) + (15.*z);

  return 0;
}
void allocate_tape_grad_q(void *tape[3]) {
    int tape_size__grad_q_x = __enzyme_augmentsize((void *)ExactSolution, enzyme_dup, enzyme_const, enzyme_dup, enzyme_const, enzyme_const);
    int tape_size__grad_q_y = __enzyme_augmentsize((void *)ExactSolution, enzyme_dup, enzyme_const, enzyme_const, enzyme_dup, enzyme_const);
    int tape_size__grad_q_z = __enzyme_augmentsize((void *)ExactSolution, enzyme_dup, enzyme_const, enzyme_const, enzyme_const, enzyme_dup);
    tape[0] = malloc(tape_size__grad_q_x);
    tape[1] = malloc(tape_size__grad_q_y);
    tape[2] = malloc(tape_size__grad_q_z);
}

void fwd_grad_q(double *q, double *t, double *x, double *y, double *z, void *tape[3]) {
    __enzyme_augmentfwd((void*)ExactSolution, 
                        enzyme_allocated, sizeof(tape[0][0]), 
                        enzyme_tape, tape[0],
                        q, (double *)NULL, 
                        enzyme_const, t,
                        x, (double *)NULL,
                        enzyme_const, y,
                        enzyme_const, z
                       );
    __enzyme_augmentfwd((void*)ExactSolution, 
                        enzyme_allocated, sizeof(tape[1][0]), 
                        enzyme_tape, tape[1],
                        q, (double *)NULL, 
                        enzyme_const, t,
                        enzyme_const, x,
                        y, (double *)NULL,
                        enzyme_const, z
                       );
    __enzyme_augmentfwd((void*)ExactSolution, 
                        enzyme_allocated, sizeof(tape[2][0]), 
                        enzyme_tape, tape[2],
                        q, (double *)NULL, 
                        enzyme_const, t,
                        enzyme_const, x,
                        enzyme_const, y,
                        z, (double *)NULL
                       );
}

void rev_grad_q_x(double *x_, double *t, double *y, double *z, void *tape) {
    for (int i=0; i<5; i++) {
        double q_[5] = {0.}; q_[i] = 1.;
        __enzyme_reverse((void*)ExactSolution, 
                         enzyme_allocated, sizeof(tape[0]), 
                         enzyme_tape, tape, 
                         (double *)NULL, q_,
                         enzyme_const, t, 
                         (double *)NULL, &x_[i],
                         enzyme_const, y,
                         enzyme_const, z
                        );
    }
}
void rev_grad_q_y(double *y_, double *t, double *x, double *z, void *tape) {
    for (int i=0; i<5; i++) {
        double q_[5] = {0.}; q_[i] = 1.;
        __enzyme_reverse((void*)ExactSolution, 
                         enzyme_allocated, sizeof(tape[0]), 
                         enzyme_tape, tape, 
                         (double *)NULL, q_,
                         enzyme_const, t, 
                         enzyme_const, x,
                         (double *)NULL, &y_[i],
                         enzyme_const, z
                        );
    }
}
void rev_grad_q_z(double *z_, double *t, double *x, double *y, void *tape) {
    for (int i=0; i<5; i++) {
        double q_[5] = {0.}; q_[i] = 1.;
        __enzyme_reverse((void*)ExactSolution, 
                         enzyme_allocated, sizeof(tape[0]), 
                         enzyme_tape, tape, 
                         (double *)NULL, q_,
                         enzyme_const, t, 
                         enzyme_const, x,
                         enzyme_const, y,
                         (double *)NULL, &z_[i]
                        );
    }
}

void computeQdot(double *q_dot, double *t, double *x, double *y, double *z) {
    double q[5];
    for (int i=0; i<5; i++) {
        double q_[5] = {0.}; q_[i] = 1.;
        __enzyme_autodiff((void *)ExactSolution,
                          q, q_,
                          t, &q_dot[i],
                          enzyme_const, x,
                          enzyme_const, y,
                          enzyme_const, z);
    }
}

int main() {
  // Declarations
  double X[] = {.5, 2.5, 5.};
  double x[1] = {X[0]}, y[1] = {X[1]}, z[1] = {X[2]};
  double time[1] = {.2};
  // -------------------------------------------------------------------------
  // q_dot
  // -------------------------------------------------------------------------
  double q_dot[5] = {0.};
  computeQdot(q_dot, time, x, y, z);
  // Print output
  printf("\nqdot\n");
  for (int j=0; j<5; j++) printf("%f\t", q_dot[j]);
  printf("\n\n");
  // -------------------------------------------------------------------------
  // grad_q
  // -------------------------------------------------------------------------
  // Allocate memory for tape
  void *tape[3];
  allocate_tape_grad_q(tape);
  // Forward mode
  double q[5];
  fwd_grad_q(q, time, x, y, z, tape);
  // Reverse mode
  double grad_q[3][5] = {{0.}};
  rev_grad_q_x(grad_q[0], time, y, z, tape[0]);
  rev_grad_q_y(grad_q[1], time, x, z, tape[1]);
  rev_grad_q_z(grad_q[2], time, x, y, tape[2]);
  // Print output
  for (int i=0; i<3; i++) {
      printf("\nDerivative in the %d direction:\n", i);
    for (int j=0; j<5; j++) printf("%f\t", grad_q[i][j]);
    printf("\n\n");
  }
  return 0;
}

/*

Output:

qdot
14.000000	28.000000	42.000000	56.000000	70.000000	


Derivative in the 0 direction:
1.000000	2.000000	3.000000	4.000000	5.000000	


Derivative in the 1 direction:
6.000000	7.000000	8.000000	9.000000	10.000000	


Derivative in the 2 direction:
11.000000	12.000000	13.000000	14.000000	15.000000
*/
