#include "headers.h"

MatrixXd Inv_M(const MatrixXd& M_proto, const MatrixXd& R, const MatrixXd& Area, int p, int s){
    MatrixXd V = R;
    int Ne = Area.rows();
    int N = (p+1)*(p+2)/2;
    for(int i=0; i<Ne; i++){
        MatrixXd tmp = V.block(i*s*N, 0, s*N, 1);
        if(N>1){
            Map<MatrixXd> B(tmp.data(), s, N); 
            MatrixXd b = B.transpose();
            // MatrixXd b(N, s);
            // for(int j=0; j<N; j++){
            //     for(int k=0; k<s; k++){
            //         b(j,k) = tmp(j*s + k);
            //     }
            // }
            tmp = M_proto.colPivHouseholderQr().solve(b);
            tmp.transposeInPlace();
            Map<MatrixXd> x(tmp.data(), s*N, 1);
            V.block(i*s*N, 0, s*N, 1) = x/Area(i)/2.0; // 2*Area(i) is the determinant of Jacobian
        }
        else{
            MatrixXd b = tmp;
            b.transposeInPlace();
            tmp = M_proto.colPivHouseholderQr().solve(b);
            tmp.transposeInPlace();
            MatrixXd x = tmp;
            V.block(i*s*N, 0, s*N, 1) = x/Area(i)/2.0; // 2*Area(i) is the determinant of Jacobian
        }
    }
    return V;
}