#include <iostream>
#include<Eigen/Dense>
#include <cmath>

using namespace std;
using namespace Eigen;

Eigen::Matrix3d S(Eigen::Vector3d x){
    Eigen::Matrix3d S;
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
};

C1C2Derivs c1c2derivs(double a) {
    constexpr double eps = 1e-8;
    double c1, c2, dc1, dc2;

    if (std::abs(a) < eps){
        c1 = 1.0 - a*a/6.0 + std::pow(a,4)/120.0;
        c2 = 0.5 - a*a/24.0 + std::pow(a,4)/720.0;
        dc1 = -a/3.0 + std::pow(a,3)/30.0;
        dc2 = -1.0/12.0 + a/120.0;
    }
    else{
        c1 = std::sin(a) / a;
        c2 = (1.0 - std::cos(a)) / (a*a);
        dc1 = (a*std::cos(a) - std::sin(a)) / (a*a);
        dc2 = (a*std::sin(a) - 2*(1-std::cos(a))) / (a*a*a);
    }
    return {c1, c2, dc1, dc2};
}

// We need to implement a newton iteration method
// The iteration method is to solve for the Lie Algebra of F_k, f_k, which is a 3d vector in R3. 
// constants that need to be set: J, Pi_k, M_k
// We also assume that the inertia tensor J has been diagonalised

Eigen::Vector3d f_newton(Eigen::Vector3d J_diag, Eigen::Vector3d Pi_k, Eigen::Vector3d M_k, double h, double tol, int maxit){
    Eigen::Vector3d b;
    b = h * Pi_k + ((h*h)/2) * M_k;

    //Initial value of f
    Eigen::Vector3d f = b.cwiseQuotient(J_diag);

    //Inertia Tensor
    Eigen::Matrix3d J = J_diag.asDiagonal();

    //Iteration
    for (int i=0; i < maxit; i++){
        double a = f.norm();
        
        auto res = c1c2derivs(a);
        double c1 = res.c1;
        double c2 = res.c2;
        double dc1 = res.dc1;
        double dc2 = res.dc2;
        Eigen::Vector3d Jf = J_diag.cwiseProduct(f);
        Eigen::Vector3d cross = f.cross(Jf);

        // G(f)
        Eigen::Vector3d G = c1 * Jf + c2 * cross - b;

        if (G.norm() < tol){
            return f;
        }

        // Jacobian DG
        Eigen::Matrix<double,1,3> f_transpose = f.transpose();
        Eigen::Matrix3d S_terms = -S(Jf) + S(f) * J;

        Eigen::Matrix3d DG = c1 * J + c2 * S_terms;

        if (a>0){
            DG += (dc1 * (Jf * f_transpose) / a) + dc2 * (cross * f_transpose) / a; 
        }

        Eigen::Vector3d delta = DG.partialPivLu().solve(G);
        Eigen::Vector3d f_new;
        f_new = f - delta;

        if (delta.norm() < 1e-12){
            return f_new;
        }

        f = f_new;
    }
    return f;
}
