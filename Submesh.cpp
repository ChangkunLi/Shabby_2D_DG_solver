#include "headers.h"

MatrixXd Submesh(const MatrixXd& U, int p, int s, const Res_data_t& res_data, const mesh_t& mesh, int nplot){
    // nplot: number of additional plotting points on every edge

    int N = (p+1)*(p+2)/2;  // number of basis
    int Ne = U.rows()/N/s;  // number of elements
    int nq = nplot + 2;     // number of plotting points on every edge
    int Nq = (1+nq)*nq/2;   // total number of plotting points in each element
    int Nplot = Ne*(nq-1)*(nq-1);   // total number of plotting elements

    MatrixXd output(Nplot, 9);  // output matrix formatted this way: x1,y1,c1,x2,y2,c2,x3,y3,c3 (every row)

    ArrayXd xi = ArrayXd::LinSpaced(nq,0,1); 
    ArrayXd eta = xi;

    MatrixXd A(Nq, N);
    int k = 0;
    for(int j=0; j<nq; j++){    
        for(int i=0; i<(nq-j); i++){    // loop over plotting points
            int t = 0;
            for(int s=0; s<=p; s++){
                for(int r=0; r<=(p-s); r++){ // loop over basis terms
                    A(k,t) = pow(xi(i), double(r)) * pow(eta(j), double(s));
                    t += 1;
                }
            }
            k += 1;
        }
    }
    MatrixXd Phi_q = A*res_data.C;

    // submesh in every element
    for(int i=0; i<Ne; i++){ 
        MatrixXd Ui = U.block(i*s*N, 0, s*N, 1);
        Map<MatrixXd> ui(Ui.data(), s, N);
        MatrixXd uq = Phi_q*(ui.transpose()); // states at plotting nodes
        MatrixXd xyq(Nq, 2);        // x, y coordinates of plotting nodes
        int kk = 0;
        double x1,y1,x2,y2,x3,y3;
        x1 = mesh.Node(mesh.Elem(i,0)-1,0);   y1 = mesh.Node(mesh.Elem(i,0)-1,1);
        x2 = mesh.Node(mesh.Elem(i,1)-1,0);   y2 = mesh.Node(mesh.Elem(i,1)-1,1);
        x3 = mesh.Node(mesh.Elem(i,2)-1,0);   y3 = mesh.Node(mesh.Elem(i,2)-1,1);
        for(int jj=0; jj<nq; jj++){
            for(int ii=0; ii<(nq-jj); ii++){
                xyq(kk, 0) = x1*(1.0 - xi(ii) - eta(jj)) + x2*xi(ii) + x3*eta(jj);
                xyq(kk, 1) = y1*(1.0 - xi(ii) - eta(jj)) + y2*xi(ii) + y3*eta(jj);
                kk += 1;
            }
        }
        int nelem = (nq-1)*(nq-1);  // number of sub-element in each element
        MatrixXd output_i(nelem, 9);
        
        // counter-clockwise triangles
        int node = 0;
        int q1, q2, q3;
        int counter = 0;    // count the line number of output_i
        for(int j=(nq-1); j>0; j--){
            for(int m=0; m<j; m++){
                q1 = node;
                q2 = node+1;
                q3 = node+1+j;
                output_i(counter, 0) = xyq(q1,0); output_i(counter, 1) = xyq(q1,1); output_i(counter, 2) = uq(q1,0);
                output_i(counter, 3) = xyq(q2,0); output_i(counter, 4) = xyq(q2,1); output_i(counter, 5) = uq(q2,0);
                output_i(counter, 6) = xyq(q3,0); output_i(counter, 7) = xyq(q3,1); output_i(counter, 8) = uq(q3,0);
                counter += 1; // move to next plotting element
                node += 1;  // jump to next point
            }
            node += 1; // jump to next line
        }

        // clockwise triangles
        node = nq;
        for(int j=(nq-2); j>0; j--){
            for(int m=0; m<j; m++){
                q1 = node;
                q2 = node+1;
                q3 = node-j-1;
                output_i(counter, 0) = xyq(q1,0); output_i(counter, 1) = xyq(q1,1); output_i(counter, 2) = uq(q1,0);
                output_i(counter, 3) = xyq(q2,0); output_i(counter, 4) = xyq(q2,1); output_i(counter, 5) = uq(q2,0);
                output_i(counter, 6) = xyq(q3,0); output_i(counter, 7) = xyq(q3,1); output_i(counter, 8) = uq(q3,0);
                counter += 1; // move to next plotting element
                node += 1;  // jump to next point
            }
            node += 1; // jump to next line
        }

        output.block(i*nelem, 0, nelem, 9) = output_i;
    }

    return output;
}