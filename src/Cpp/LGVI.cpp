#include <iostream>
#include <Eigen/Dense>

Eigen::Matrix3d S(const Eigen::Vector3d& x) {
    Eigen::Matrix3d S;
    S << 0, -x(2), x(1),
        x(2), 0, -x(0),
        -x(1), x(0), 0;
    return S;
}