// ------------------------------------------------------------------------------------------------------------------ //
// Example showing the differentiation a single-input single-output function
//
// Author: Roberto Agromayor
// Date: January 2019
//
// ------------------------------------------------------------------------------------------------------------------ //


// Include libraries
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <adolc/adolc.h>


// Define namespaces
using namespace std;                // cout, endl
using namespace std::chrono;        // nanoseconds, system_clock, seconds


// Quick factorial implementation
int factorial(int n) { return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n; }


// Define the function to be differentiated: f(x) = e^[(x0+x1+...+xn)/n]
adouble my_function(adouble *x) {
  adouble u = x[0], v= x[1], w = x[2];
  return 5*u*u*v*w + 3*u*v*v*w + 6*w*w*v + 1.;
  //return 4*u*u*u + 3*v*v*v + w*w*w;
  //return u*u*v*v*w*w;
};

int main() {


    // -------------------------------------------------------------------------------------------------------------- //
    // Initialize problem variables
    // -------------------------------------------------------------------------------------------------------------- //

    // Initialize passive variables
    int m = 1, n = 3;
    auto xp = new double[n];    // Independent vector
    auto yp = new double[m];    // Dependent vector

    // Set the value of the independent variables
    for (int i=0; i<n; i++) xp[i] = i + 1.;

    // Initialize active variables
    auto x = new adouble[n];
    auto y = new adouble[m];



    // -------------------------------------------------------------------------------------------------------------- //
    // Active section for automatic differentiation
    // -------------------------------------------------------------------------------------------------------------- //

    // Set the tag for the Automatic Differentiation trace
    int tag = 0;

    // Start tracing floating point operations
    trace_on(tag);  // Start of the active section

    // Assign independent variables
    for (int i = 0; i < n; ++i) x[i] <<= xp[i];

    // Evaluate the body of the differentiated code
    y[0] = my_function(x);

    // Assign dependent variables
    y[0] >>= yp[0];

    trace_off();    // End of the active section



    // -------------------------------------------------------------------------------------------------------------- //
    // Compute the higher order derivatives in forward vector mode
    // -------------------------------------------------------------------------------------------------------------- //

    // Define the number of directions (same as the number of independent variables)
    int p = n;

    // Define the order of the highest derivative
    int degree = 2;

    // Initialize the matrix of tangent directions (identity matrix)
    double ***X = myalloc3(n, p, degree);
    for (int i = 0; i < n; ++i) {
        for (int dir = 0; dir < p; ++dir) {
            if (i == dir) {
                X[i][dir][0] = 1.00;  // Last index [0], activate only the x1 direction (see eq. 8)
            }
            else {
                X[i][dir][0] = 0.00;  // Last index [0], activate only the x1 direction (see eq. 8)
            }
        }
    }

    // Declare the matrix of derivatives
    double ***Y = myalloc(m, p, degree);

    // Compute the taylor coefficients and store them in the array Y
    hov_forward(tag, m, n, degree, p, xp, X, yp, Y);

    // Scale the Taylor coefficients to compute the derivatives
    for (int deg = 0; deg < degree; ++deg) {
        for (int dir = 0; dir < n; ++dir) {
            Y[0][dir][deg] = Y[0][dir][deg]*factorial(deg+1);
        }
    }

    // Compare the AD and the analytic derivatives
    cout << "Derivative computation in forward mode" << endl;
    cout << setw(20) << "Order" << setw(20) << "Direction" << setw(20) << "AD derivative" << setw(25) << endl;
    for (int deg = 0; deg < degree; ++deg) {
            for (int dir = 0; dir < n; ++dir) {
            cout << setw(20) << deg+1 << setw(20) << dir+1 << setw(20) << Y[0][dir][deg] << setw(25) << endl;
        }
    }



    // -------------------------------------------------------------------------------------------------------------- //
    // Run a forward sweep with keep=d+1 to prepare the trace for high-order reverse differentiation
    // -------------------------------------------------------------------------------------------------------------- //

    // Define a flag to prepare for a reverse automatic differentiation
    int keep = degree+1;

    // Declare the matrix of derivatives
    double **YY = myalloc(m, degree);

    // Initialize the matrix of tangent directions (identity matrix)
    double **XX = myalloc(n, degree);
    for (int i = 0; i < n; ++i) {
        for (int dir = 0; dir < p; ++dir) {
            if (i == dir) {
                XX[i][0] = 1.00;  // Last index [0], activate only the x1 direction (see eq. 8)
            }
            else {
                XX[i][0] = 0.00;  // Last index [0], activate only the x1 direction (see eq. 8)
            }
        }
    }

    // Run the forward mode with keep=degree+1
    hos_forward(tag, m, n, degree, keep, xp, XX, yp, YY);



    // -------------------------------------------------------------------------------------------------------------- //
    // Compute the higher order derivatives in reverse vector mode
    // -------------------------------------------------------------------------------------------------------------- //

    // Define the number of weight vectors (same as the number of independent variables)
    int q = m;

    // Declare the matrix of derivatives (adjoint matrix)
    double ***Z = myalloc(q, n, degree+1);

    // Define the weight matrix
    double **U = myalloc(q, m);
    for (int i = 0; i < q; ++i) {
        for (int j = 0; j < m; ++j) {
            if (i == j) { U[i][j] = 1.00; }
            else { U[i][j] = 0.00; }

        }
    }

    // Declare the `nonzero pattern` variable. There is no myalloc() for intergers :(
    auto **nz = new short int*[q];
    for (int k = 0; k < q; ++k) {
        nz[k] = new short int[n];
    }

    // Compute the matrix of adjoints and store them in the array Z
    hov_reverse(tag, m, n, degree, q, U, Z, nz);

    // Scale the Taylor coefficients of the adjoint matrix to compute the derivatives
    // Note that Z[..., 0] stores the function values
    for (int deg = 1; deg < degree+1; ++deg) {
        for (int dir = 0; dir < n; ++dir) {
            Z[0][dir][deg] = Z[0][dir][deg]*factorial(deg);
        }
    }

    // Compare the AD and the analytic derivatives
    cout << "Derivative computation in reverse mode" << endl;
    cout << setw(20) << "Order" << setw(20) << "Direction" << setw(20) << "AD derivative" << setw(25) << endl;
    for (int deg = 1; deg < degree+1; ++deg) {
            for (int dir = 0; dir < n; ++dir) {
            cout << setw(20) << deg << setw(20) << dir+1 << setw(20) << Z[0][dir][deg] << setw(25) << endl;
        }
    }

    return 0;

}


