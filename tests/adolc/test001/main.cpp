#include <iostream>

class DualNumber {
private:
    double value;
    double gradient;

public:
    DualNumber(double value, double gradient = 0) : value(value), gradient(gradient) {}

    // Overload addition operator
    DualNumber operator+(const DualNumber& other) const {
        // The const qualifier in the argument means that the other object is not
        //   modified  within the function, ensuring that it remains unchanged.
        // The second "const" indicates that this operator function does not modify
        //   the state of the current object (the object on the left-hand side of the + operator).
        return DualNumber(value + other.value, gradient + other.gradient);
    }

    // Overload multiplication operator
    DualNumber operator*(const DualNumber& other) const {
        // returned value = self.value * other.value
        // returned gradient = self.value * other.gradient + self.gradient * other.value
        return DualNumber(value * other.value, value * other.gradient + gradient * other.value);
    }

    // Getter for value
    double getValue() const {
        return value;
    }

    // Getter for gradient
    double getGradient() const {
        return gradient;
    }
};

// Define a function f(x, y) = x^2 * y + y^2
DualNumber f(const DualNumber& x, const DualNumber& y) {
    return x * x * y + y * y;
}

int main() {
    // Define dual numbers representing x = 3 with derivative 1, and y = 2 with derivative 0
    DualNumber x(3, 1); // x is an independent variable, dx/dx = 1
    DualNumber y(2, 0); // y is a constant and dy/dx = 0

    // Compute f(x, y) and its gradient
    DualNumber result = f(x, y);

    // Output the function value and gradient with respect to x
    std::cout << "Result: " << result.getValue() << std::endl;
    std::cout << "Gradient wrt x: " << result.getGradient() << std::endl;

    return 0;
}
