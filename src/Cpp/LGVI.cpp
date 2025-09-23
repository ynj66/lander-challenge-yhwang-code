#include <iostream>
#include<Eigen/Dense>

using namespace std;
using namespace Eigen;

Eigen::Matrix3d S(Eigen::Vector3d x){
    Matrix3d S;
    S<< 0, -x(2), x(1),
        x(2), 0, -x(0),
        -x(1), x(0), 0;
    return S;
}
