#include <iostream>
#include <Eigen/Dense>

int main(){
    Eigen::Matrix2d A;
    A(0, 0) = 0.0;
    A(0, 1) = 1.0;
    A(1, 0) = 2.0;
    A(1, 1) = 3.0;

    std::cout << A << std::endl;
    Eigen::Vector2d x;
    x << 1.0, 1.0;

    Eigen::Vector2d y = A * x;
    std::cout << y << std::endl;
}