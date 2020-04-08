#include "headers.h"

void precompute(int p, MatrixXd& C, MatrixXd& M_proto, MatrixXd& Phi, MatrixXd& dxPhi, MatrixXd& dyPhi, vector<MatrixXd>& Phi_edge_CC, vector<MatrixXd>& Phi_edge_c, MatrixXd& wq_2d, MatrixXd& wq_1d){
    // C: lagrange basis coefficients
    // M_proto: mass matrix
    // Phi: basis on 2d quad points
    // dx(y)Phi: basis gradient on 2d quad points
    // Phi_edge_c(CC): basis on edge quad points, both clockwise(c) and counter-clockwise(CC) ordering are stored

    // precompute coefficients of lagrange basis on reference element
    C = TriLagrange2D(p);    
    // precompute mass matrix and values of basis functions on 2d quad points
    M_proto = MassMatrix2D(C, p, Phi);      
    
    // precompute basis function gradient on 2d quad points
    int nq;
    int N = (p+1)*(p+2)/2;
    MatrixXd xyq_2d;
    quad2d(2*p+1, nq, xyq_2d, wq_2d);
    MatrixXd A(nq, N);
    MatrixXd B(nq, N);
    for(int i=0; i<nq; i++){
        int k = 0;
        for(int s=0; s<=p; s++){
            for(int r=0; r<=(p-s); r++){
                A(i,k) = double(r) * pow(xyq_2d(i,0), double(r-1)) * pow(xyq_2d(i,1), double(s));  // notice that no quad point is on the axis (xi, eta != 0), so raising to exponent -1 is ok.
                B(i,k) = pow(xyq_2d(i,0), double(r)) * double(s) * pow(xyq_2d(i,1), double(s-1));
                k += 1;
            }
        }
    }
    dxPhi = A*C; dyPhi = B*C;

    // precompute basis function on edge 1d quad points
    MatrixXd xyq_1d;
    quad1d(2*p+1, nq, xyq_1d, wq_1d);
    MatrixXd A_CC(nq, N);
    MatrixXd A_c(nq, N);

    // face 1:
    if(Phi_edge_c.empty() && Phi_edge_CC.empty()){
        MatrixXd eta = xyq_1d;
        MatrixXd xi  = MatrixXd::Ones(nq, 1) - eta;
        for(int i=0; i<nq; i++){
            int k = 0;
            for(int s=0; s<=p; s++){
                for(int r=0; r<=(p-s); r++){
                    A_CC(i,k) = pow(xi(i), double(r)) * pow(eta(i), double(s));
                    A_c(nq-1-i,k) = pow(xi(i), double(r)) * pow(eta(i), double(s));
                    k += 1;
                }
            }
        }
        Phi_edge_CC.push_back(A_CC*C); Phi_edge_c.push_back(A_c*C);
    }
    else{
        cout << "Error: Phi_edge_c and Phi_edge_CC are not empty vectors!" << endl;
        return;
    }

    // face 2:
    MatrixXd eta = MatrixXd::Ones(nq, 1) - xyq_1d;
    MatrixXd xi  = 0.0*eta;
    for(int i=0; i<nq; i++){
        int k = 0;
        for(int s=0; s<=p; s++){
            for(int r=0; r<=(p-s); r++){
                A_CC(i,k) = pow(xi(i), double(r)) * pow(eta(i), double(s));
                A_c(nq-1-i,k) = pow(xi(i), double(r)) * pow(eta(i), double(s));
                k += 1;
            }
        }
    }
    Phi_edge_CC.push_back(A_CC*C); Phi_edge_c.push_back(A_c*C);

    // face 3:
    xi = xyq_1d;
    eta *= 0.0;
    for(int i=0; i<nq; i++){
        int k = 0;
        for(int s=0; s<=p; s++){
            for(int r=0; r<=(p-s); r++){
                A_CC(i,k) = pow(xi(i), double(r)) * pow(eta(i), double(s));
                A_c(nq-1-i,k) = pow(xi(i), double(r)) * pow(eta(i), double(s));
                k += 1;
            }
        }
    }
    Phi_edge_CC.push_back(A_CC*C); Phi_edge_c.push_back(A_c*C);
}