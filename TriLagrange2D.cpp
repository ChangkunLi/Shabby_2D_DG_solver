#include "headers.h"

// MatrixXd C = TriLagrange2D(p);

MatrixXd TriLagrange2D(int p){
    // calculates coeffs for full-order Lagrange basis of order p
    // reference element is a unit isoceles right triangle
    ArrayXd xi = ArrayXd::LinSpaced(p+1,0,1); 
    ArrayXd eta = xi;
    int N = (p+1)*(p+2)/2; // number of basis functions
    MatrixXd A = MatrixXd::Zero(N,N); 
    int i = 0; // build A-matrix
    for(int iy=0; iy<=p; iy++){
        for(int ix=0; ix<=(p-iy); ix++){
            int k = 0;
            for(int s=0; s<=p; s++){
                for(int r=0; r<=(p-s); r++){
                    A(i,k) = pow(xi(ix), double(r)) * pow(eta(iy), double(s));
                    k += 1;
                }
            }
            i += 1;
        }
    }
    return A.inverse();
}