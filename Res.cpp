#include "headers.h"

MatrixXd Res(const MatrixXd& U, int p, int s, const string& btype, const Res_data_t& res_data){
    MatrixXd R = MatrixXd::Zero(U.rows(), U.cols());
    int N = (p+1)*(p+2)/2;  // number of basis
    int Ne = U.rows()/N/s;  // number of elements
    int nq_2d = res_data.wq_2d.rows();
    int nq_1d = res_data.wq_1d.rows();

    // loop over elements
    for(int i=0; i<Ne; i++){
        MatrixXd Pseudo_dxPhi = res_data.Pseudo_invJ[i](0,0) * res_data.dxPhi + res_data.Pseudo_invJ[i](1,0) * res_data.dyPhi;
        Pseudo_dxPhi.transposeInPlace();
        MatrixXd Pseudo_dyPhi = res_data.Pseudo_invJ[i](0,1) * res_data.dxPhi + res_data.Pseudo_invJ[i](1,1) * res_data.dyPhi;
        Pseudo_dyPhi.transposeInPlace();
        for(int j=0; j<nq_2d; j++){
            Pseudo_dxPhi.col(j) *= res_data.wq_2d(j);
            Pseudo_dyPhi.col(j) *= res_data.wq_2d(j);
        }
        MatrixXd Ui = U.block(i*s*N, 0, s*N, 1);
        Map<MatrixXd> ui(Ui.data(), s, N);
        MatrixXd uq = res_data.Phi*(ui.transpose());

        MatrixXd Fx(nq_2d, s);
        MatrixXd Fy(nq_2d, s);

        for(int j=0; j<nq_2d; j++){
            double h, hu, hv;
            h = uq(j,0);    hu = uq(j,1);   hv = uq(j,2);
            Fx(j,0) = hu;  
            Fx(j,1) = hu*hu/h + g/2.0*h*h;
            Fx(j,2) = hu*hv/h;
            Fy(j,0) = hv;
            Fy(j,1) = hu*hv/h;
            Fy(j,2) = hv*hv/h + g/2.0*h*h;
        }

        MatrixXd R_i = Pseudo_dxPhi*Fx + Pseudo_dyPhi*Fy;
        R_i.transposeInPlace();
        Map<MatrixXd> r_i(R_i.data(), s*N, 1);
        R.block(i*s*N, 0, s*N, 1) -= r_i;
    }

    // loop over interior edges
    for(int ie=0; ie<res_data.I2E.rows(); ie++){
        int elemL = res_data.I2E(ie, 0);
        int faceL = res_data.I2E(ie, 1);
        int elemR = res_data.I2E(ie, 2);
        int faceR = res_data.I2E(ie, 3);
        MatrixXd Phi_L = res_data.Phi_edge_CC[faceL - 1];
        MatrixXd Phi_R = res_data.Phi_edge_c[faceR - 1];
        MatrixXd UL = U.block((elemL-1)*s*N, 0, s*N, 1);
        MatrixXd UR = U.block((elemR-1)*s*N, 0, s*N, 1);
        Map<MatrixXd> uL(UL.data(), s, N);
        Map<MatrixXd> uR(UR.data(), s, N);
        MatrixXd uqL = Phi_L*(uL.transpose());
        MatrixXd uqR = Phi_R*(uR.transpose());
        MatrixXd F_hat(s, nq_1d);

        double smax;
        double nx = res_data.In(ie, 0);
        double ny = res_data.In(ie, 1);
        for(int j=0; j<nq_1d; j++){
            F_hat.col(j) = roe(uqL.row(j), uqR.row(j), nx, ny, smax);
        }
        F_hat.transposeInPlace();

        // left cell residual
        Phi_L.transposeInPlace();
        for(int j=0; j<nq_1d; j++){
            Phi_L.col(j) *= res_data.wq_1d(j);
        }
        MatrixXd RL = Phi_L*F_hat*res_data.Il(ie);
        RL.transposeInPlace();
        Map<MatrixXd> rL(RL.data(), s*N, 1);
        R.block((elemL-1)*s*N, 0, s*N, 1) += rL;

        // right cell residual
        Phi_R.transposeInPlace();
        for(int j=0; j<nq_1d; j++){
            Phi_R.col(j) *= res_data.wq_1d(j);
        }
        MatrixXd RR = Phi_R*F_hat*res_data.Il(ie);
        RR.transposeInPlace();
        Map<MatrixXd> rR(RR.data(), s*N, 1);
        R.block((elemR-1)*s*N, 0, s*N, 1) -= rR;
    }

    // loop over boundary edges
    if(btype == "Full"){
        for(int be=0; be<res_data.B2E.rows(); be++){
            int elemL = res_data.B2E(be, 0);
            int faceL = res_data.B2E(be, 1);
            MatrixXd Phi_L = res_data.Phi_edge_CC[faceL - 1];
            MatrixXd UL = U.block((elemL-1)*s*N, 0, s*N, 1);
            Map<MatrixXd> uL(UL.data(), s, N);
            MatrixXd uqL = Phi_L*(uL.transpose());
            MatrixXd F_hat(s, nq_1d);

            MatrixXd uqR = MatrixXd::Ones(nq_1d, s);
            uqR.col(1) *= 0.453;
            uqR.col(2) *= 0.769;

            double smax;
            double nx = res_data.Bn(be, 0);
            double ny = res_data.Bn(be, 1);
            for(int j=0; j<nq_1d; j++){
                F_hat.col(j) = roe(uqL.row(j), uqR.row(j), nx, ny, smax);
            }
            F_hat.transposeInPlace();

            Phi_L.transposeInPlace();
            for(int j=0; j<nq_1d; j++){
                Phi_L.col(j) *= res_data.wq_1d(j);
            }
            MatrixXd RL = Phi_L*F_hat*res_data.Bl(be);
            RL.transposeInPlace();
            Map<MatrixXd> rL(RL.data(), s*N, 1);
            R.block((elemL-1)*s*N, 0, s*N, 1) += rL;
        }
    }
    else if(btype == "Wall"){
        for(int be=0; be<res_data.B2E.rows(); be++){
            int elemL = res_data.B2E(be, 0);
            int faceL = res_data.B2E(be, 1);
            MatrixXd Phi_L = res_data.Phi_edge_CC[faceL - 1];
            MatrixXd UL = U.block((elemL-1)*s*N, 0, s*N, 1);
            Map<MatrixXd> uL(UL.data(), s, N);
            MatrixXd uqL = Phi_L*(uL.transpose());
            MatrixXd F_hat(nq_1d, s);

            double nx = res_data.Bn(be, 0);
            double ny = res_data.Bn(be, 1);
            for(int j=0; j<nq_1d; j++){
                double h = uqL(j, 0);
                F_hat(j,0) = 0.0;
                F_hat(j,1) = nx/2.0*g*h*h;
                F_hat(j,2) = ny/2.0*g*h*h;
            }

            Phi_L.transposeInPlace();
            for(int j=0; j<nq_1d; j++){
                Phi_L.col(j) *= res_data.wq_1d(j);
            }
            MatrixXd RL = Phi_L*F_hat*res_data.Bl(be);
            RL.transposeInPlace();
            Map<MatrixXd> rL(RL.data(), s*N, 1);
            R.block((elemL-1)*s*N, 0, s*N, 1) += rL;
        }
    }
    else{
        cout << "Unsupported boundary type" << endl;
        return 0.0*U;
    }

    return R;
}