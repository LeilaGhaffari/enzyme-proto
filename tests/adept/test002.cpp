#include <vector>
#include <adept.h>
#include <adept/array_shortcuts.h>

adept::aVector func(const adept::aVector& x) {
    adept::aVector y(4);
    for (int i=0; i<y.size(); i++) {
        y[i] = 1. * (i+1);
        for (int j=0; j<x.size(); j++) y[i] *= 4.* x[j];
    }
    return y;
}

void func_jacobian(const adept::Vector& x_val, adept::Vector& y_val, adept::Matrix& jac) {
    adept::Stack stack;
    adept::aVector x = x_val;
    stack.new_recording();
    adept::aVector y = func(x);
    stack.independent(x);
    stack.dependent(y);
    jac = stack.jacobian(); // m×n
    y_val = value(y);
}

int main() {
    const adept::Vector& x = {1, 2, 3};
    adept::Vector y;
    adept::Matrix jac;
    func_jacobian(x, y, jac);
    for (int i=0; i<y.size(); i++) {
        std::cout << "y[" << i << "] = " << y[i] << "\t";
        for (int j=0; j<x.size(); j++)
            std::cout << "J[" << i << "]" << "[" << j << "] = " <<  jac[i][j] << "\t";
        std::cout << std::endl;
    }
    return 0;
}

/*
g++ -std=c++11 -o exec test002.cpp -ladept

y[0] = 384      J[0][0] = 384   J[0][1] = 192   J[0][2] = 128
y[1] = 768      J[1][0] = 768   J[1][1] = 384   J[1][2] = 256
y[2] = 1152     J[2][0] = 1152  J[2][1] = 576   J[2][2] = 384
y[3] = 1536     J[3][0] = 1536  J[3][1] = 768   J[3][2] = 512
*/
