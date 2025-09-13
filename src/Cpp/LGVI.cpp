#include <iostream>
#include <Eigen/Dense>

Eigen::Matrix3d S(const Eigen::Vector3d& x) {
    Eigen::Matrix3d S;
    S << 0, -x(2), x(1),
        x(2), 0, -x(0),
        -x(1), x(0), 0;
    return S;
}

int main() {
    Eigen::VectorXd v(3); // Declare a 3-element column vector of doubles
    v(0) = 1.0;
    v(1) = 2.5;
    v(2) = 3.0;
    std::cout << "Vector v:\n" << v << std::endl;
    return 0;
}