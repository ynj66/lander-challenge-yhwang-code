#include <iostream>
#include<Eigen/Dense>
#include <cmath>

using namespace std;
using namespace Eigen;

Eigen::Matrix3d S(Eigen::Vector3d x){
    Matrix3d S;
    S<< 0, -x(2), x(1),
        x(2), 0, -x(0),
        -x(1), x(0), 0;
    return S;
}

struct C1C2Derivs{
    double c1;
    double c2;
    double dc1;
    double dc2;
}

C1C2Derivs c1c2derivs(double a) {
    constexpr double eps = 1e-8
    double c1, c2, dc1, dc2;

    if (std::abs(a) < eps){
        c1 = 1.0 - a*a/6.0 + std::pow(a,4)/120.0;
        c2 = 0.5 - a*a/24.0 + std::pow(a,4)/720.0;
        dc1 = -a/3.0 + std::pow(a,3)/30.0;
        dc2 = -1.0/12.0 + a/120.0;
    }
    else{
        c1 = std::sin(a) / a;
        c2 = (1.0 = std::cos(a)) / (a*a)
        dc1 = (a*std::cos(a) - std::sin(a)) / (a*a);
        dc2 = (a*std::sin(a) - 2*(1-std::cos(a))) / (a*a*a);
    }
    return {c1, c2, dc1, dc2}
}

int main() {
    double a = 1e-9;
    auto res = c1c2derivs(a);
    std::cout << "c1=" << res.c1 
              << ", c2=" << res.c2 
              << ", dc1=" << res.dc1 
              << ", dc2=" << res.dc2 
              << std::endl;
}