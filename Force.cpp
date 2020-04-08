#include "headers.h"

MatrixXd Force(const MatrixXd& U, int p, int s, const Res_data_t& res_data){
    MatrixXd force = MatrixXd::Zero(1,8);
    int N = (p+1)*(p+2)/2;
    int nq_1d = res_data.wq_1d.rows();

    for(int be=0; be<res_data.B2E.rows(); be++){
        int elemL = res_data.B2E(be, 0);
        int faceL = res_data.B2E(be, 1);
        MatrixXd Phi_L = res_data.Phi_edge_CC[faceL - 1];
        MatrixXd UL = U.block((elemL-1)*s*N, 0, s*N, 1);
        Map<MatrixXd> uL(UL.data(), s, N);
        MatrixXd uqL = Phi_L*(uL.transpose());
        MatrixXd h = uqL.col(0);

        double nx = res_data.Bn(be, 0);
        double ny = res_data.Bn(be, 1);
        
        if(res_data.B2E(be,2) > 1){
            int nbuilding = res_data.B2E(be,2) - 2;
            for(int j=0; j<nq_1d; j++){
                double prod = rho*g/2.0*pow(h(j), 2.0)*res_data.wq_1d(j)*res_data.Bl(be);
                force(2*nbuilding) += prod*nx;
                force(2*nbuilding+1) += prod*ny;
            }
        }
    }
    return force;
}