
#include <adolc/adolc.h>

#include <stdlib.h>
#include <iostream>
using namespace std;

static int n = 3, m = 1, d = 2, p = 3;

double func(double x[3]) {
    double u = x[0], v= x[1], w = x[2];
    return 5*u*u*v*w + 3*u*v*v*w + 6*w*w*v + 1.;
}

adouble func(adouble x[3]) {
    adouble u = x[0], v= x[1], w = x[2];
    return 5*u*u*v*w + 3*u*v*v*w + 6*w*w*v + 1.;
}

double *dfunc(double xp[3]) {
    int tag = 1;
    auto xa = new adouble[n];
    auto ya = new adouble[m];
    auto yp = new double[m];
    auto dX = new double[n];
    double **X = myalloc(n, p);
    double **Y = myalloc(m, p);

    trace_on(tag);
    for (int i=0; i<n; i++) xa[i] <<= xp[i];
    ya[0] = func(xa);
    ya[0] >>= yp[0];
    trace_off();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) if (i == j) X[i][j] = 1.;
        else  X[i][j] = 0.00;
    }

    fov_forward(tag, m, n, p, xp, X, yp, Y);

    for (int i=0; i<n; i++) dX[i] = Y[0][i];

    return dX;
}

double **dfunc_ho(double xp[3]) {
    int tag = 1;
    auto yp = new double[m];
    adouble* x = new adouble[n];
    adouble* y = new adouble[m];
    double** S = new double*[n];
    double** dtensor; // size = m x dim

    for (int i=0; i<n; i++) { // Identity matrix
        S[i] = new double[p];
        for (int j=0; j<p; j++) S[i][j] = (i == j) ? 1.0 : 0.0;
    }

    trace_on(tag);
    for (int i=0; i<n; i++) x[i] <<= xp[i];
    y[0] = func(x);
    y[0] >>= yp[0];
    trace_off();

    int dim = binomi(p + d, d);
    dtensor = myalloc2(m, dim);
    tensor_eval(tag, m, n, d, p, xp, dtensor, S);

    return dtensor;
}


int main() {
    int dim = binomi(p + d, d), temp = 1;
    double* xp = new double[n];
    for (int i=0; i<n; i++) xp[i] = i + 1.0;

    // Function evaluation
    auto f = func(xp);

    // First Derivative Vector
    auto dX = dfunc(xp);

    // Higher Derivative Tensor
    auto dtensor = dfunc_ho(xp);
    double F = dtensor[0][0];
    double dF[n] ={0.};
    double d2F[n][n] ={{0.}};
    for (int i=0; i<n; i++) {
        dF[i] = dtensor[0][temp];
        for (int j=0; j<i+1; j++) {
            d2F[i][j] = dtensor[0][temp+j+1];
            if (i != j) d2F[j][i] = d2F[i][j];
        }
        temp += i + 2;
    }

    // Print
    cout << "\nFunction Evaluation:\n\t" << f << endl << endl;
    cout << "First Order Forward Mode Derivative:" << endl;
    for (int i=0; i<n; i++) cout << "\t" << dX[i] << endl;
    cout << endl;

    cout << "Higher Derivative Tensor:\n--------------------------\n" << endl;
    cout << "F = " << F << endl;

    cout << "\ndF = " << endl;
    for (int i=0; i<n; i++) cout << dF[i] << "\t";
    cout << endl;

    cout << "\ndF = " << endl;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) cout << d2F[i][j] << "\t";
        cout << endl;
    }
    cout << endl;

    // Clean up
    myfree2(dtensor);

    return 0;
}

/****************************************************************************/

/*

*/