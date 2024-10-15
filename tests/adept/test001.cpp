#include <iostream>
#include <adept_arrays.h>

// example from Wikipedia (https://en.wikipedia.org/wiki/Adept_(C%2B%2B_library)#Example)

int main(int argc, const char** argv) {
  adept::Stack stack;                           // Object to store differential statements
  adept::aVector x(3);                          // Independent variables: active vector with 3 elements
  x << 1.0, 2.0, 3.0;                           // Fill vector x
  stack.new_recording();                        // Clear any existing differential statements
  adept::adouble J = cbrt(sum(abs(x * x * x))); // Compute dependent variable: 3-norm in this case
  J.set_gradient(1.0);                          // Seed the dependent variable
  stack.reverse();                              // Reverse-mode differentiation
  std::cout << "dJ/dx = "
            << x.get_gradient() << "\n";        // Print the vector of partial derivatives dJ/dx

  return 0;
}

// Build with: g++ -std=c++11 -o exec test001.cpp -ladept
