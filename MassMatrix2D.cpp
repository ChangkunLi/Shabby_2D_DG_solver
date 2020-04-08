#include "headers.h"

MatrixXd MassMatrix2D(const MatrixXd& C, int p, MatrixXd& Phi){
    int N = (p+1)*(p+2)/2;
    MatrixXd M_proto;
    
    // compute Phi
    int nq;
    MatrixXd xyq_2d,wq_2d;
    quad2d(2*p+1, nq, xyq_2d, wq_2d);
    MatrixXd A(nq, N);
    for(int i=0; i<nq; i++){
        int k = 0;
        for(int s=0; s<=p; s++){
            for(int r=0; r<=(p-s); r++){
                A(i,k) = pow(xyq_2d(i,0), double(r)) * pow(xyq_2d(i,1), double(s));
                k += 1;
            }
        }
    }
    Phi = A*C;

    // compute M_proto
    M_proto = Phi.transpose();
    for(int i=0; i<nq; i++){
        M_proto.col(i) *= wq_2d(i);
    }
    M_proto *= Phi;

    return M_proto;
}