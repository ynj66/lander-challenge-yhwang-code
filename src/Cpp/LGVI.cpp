#include <iostream>
#include<Eigen/Dense>

using namespace std;
using namespace Eigen;


int main(){

    Matrix <float, 2,2> A;
    A << 1,2,
         3,4;
    cout << A << endl;
    cout << "Hello WOrld!" << endl;
    return 0;
}