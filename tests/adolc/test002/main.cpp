#include <iostream>

class DualNumber {
private:
    double value;
    double gradient_x;
    double gradient_y;

public:
    DualNumber(double value, double gradient_x = 0, double gradient_y = 0)
        : value(value), gradient_x(gradient_x), gradient_y(gradient_y) {}

    double getValue() const { return value; }
    double getGradientX() const { return gradient_x; }
    double getGradientY() const { return gradient_y; }

    DualNumber operator+(const DualNumber& other) const {
        return DualNumber(value + other.value, gradient_x + other.gradient_x, gradient_y + other.gradient_y);
    }

    DualNumber operator*(const DualNumber& other) const {
        return DualNumber(value * other.value,
                          gradient_x * other.value + value * other.gradient_x,
                          gradient_y * other.value + value * other.gradient_y);
    }
};

DualNumber f(const DualNumber& x, const DualNumber& y) {
    DualNumber term1 = x * x; // x^2
    DualNumber term2 = y;     // y
    DualNumber term3 = y * y; // y^2
    return term1 * term2 + term3; // x^2 * y + y^2
}

int main() {
    // Define variables x and y
    DualNumber x(3, 1, 0); // x = 3, gradient with respect to x = 1
    DualNumber y(2, 0, 1); // y = 2, gradient with respect to y = 1

    // Compute the result of f(x, y)
    DualNumber result = f(x, y);

    // Output the result and gradients
    std::cout << "Result: " << result.getValue() << std::endl;
    std::cout << "Gradient wrt x: " << result.getGradientX() << std::endl;
    std::cout << "Gradient wrt y: " << result.getGradientY() << std::endl;

    return 0;
}
