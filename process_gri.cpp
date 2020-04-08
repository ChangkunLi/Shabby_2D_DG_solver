#include "headers.h"

void process_gri(const string& filename, MatrixXi& I2E, MatrixXi& B2E, MatrixXd& In, MatrixXd& Bn, MatrixXd& Area, MatrixXi& IE, MatrixXi& BE){
    // read mesh (.gri file)
    mesh_t mesh;
    read_gri(filename, mesh);

    edgehash(mesh.Elem, IE, BE);

    // process interior edge
    int IETot = IE.rows();
    I2E = MatrixXi::Zero(IETot, 4);
    In  = MatrixXd::Zero(IETot, 2);
    I2E.col(0) = IE.col(2);
    I2E.col(2) = IE.col(3);
    for(int i=0; i<IETot; i++){
        int n1 = IE(i,0);
        int n2 = IE(i,1);
        int elemL = I2E(i,0);
        int elemR = I2E(i,2);
        MatrixXi nL = mesh.Elem.row(elemL - 1);
        MatrixXi nR = mesh.Elem.row(elemR - 1);
        for(int iface=1; iface<=3; iface++){
            if((nL(iface-1)!=n1) && (nL(iface-1)!=n2)){
                I2E(i,1) = iface;
            }
        }
        for(int iface=1; iface<=3; iface++){
            if((nR(iface-1)!=n1) && (nR(iface-1)!=n2)){
                I2E(i,3) = iface;
            }
        }

        MatrixXd r1 = mesh.Node.row(n1-1);
        MatrixXd r2 = mesh.Node.row(n2-1);
        MatrixXd dr = r1 - r2;
        In(i,0) =  dr(1)/sqrt(dr(0)*dr(0) + dr(1)*dr(1));
        In(i,1) = -dr(0)/sqrt(dr(0)*dr(0) + dr(1)*dr(1));
    }

    // process boundary edge
    SparseMatrix<int> H(mesh.nNode,mesh.nNode);

    // Loop over elements and identify all edges
    for(int elem=1; elem<=mesh.nElem; elem++){
        MatrixXi nv = mesh.Elem.row(elem-1);
        for(int edge=1; edge<=3; edge++){
            int n1 = nv(edge % 3 );
            int n2 = nv((edge+1) % 3);
            if(H.coeffRef(n1-1,n2-1) == 0){ // edge hit for the first time
                // could be a boundary or interior; assume boundary
                H.coeffRef(n1-1,n2-1) = elem; H.coeffRef(n2-1,n1-1) = elem;
            }
            else{ // this is an interior edge, hit for the second time
                int oldelem = H.coeffRef(n1-1,n2-1);
                if(oldelem<0){
                    cout << "Mesh input error" << endl;
                    return;
                }
                H.coeffRef(n1-1,n2-1) = -1; H.coeffRef(n2-1,n1-1) = -1;
            }
        }
    }

    int nBFaceTot = mesh.B.nbface.sum();
    B2E = MatrixXi::Zero(nBFaceTot, 3);
    Bn  = MatrixXd::Zero(nBFaceTot, 2);
    int count = -1;
    for(int i=0; i<mesh.B.nbfgrp; i++){
        for(int j=0; j<mesh.B.nbface(i); j++){
            count += 1;
            int n1 = mesh.B.nodes[i](j,0);
            int n2 = mesh.B.nodes[i](j,1);
            B2E(count, 0) = H.coeff(n1-1,n2-1);
            MatrixXi nv = mesh.Elem.row(B2E(count,0)-1);
            for(int iface=1; iface<=3; iface++){
                if((nv(iface-1)!=n1) && (nv(iface-1)!=n2)){
                    B2E(count, 1) = iface;
                }
            }
            B2E(count, 2) = i+1;

            MatrixXd dr = mesh.Node.row(n1-1) - mesh.Node.row(n2-1);
            if(nv(B2E(count,1) % 3) == n1){
                Bn(count,0) = -dr(1)/sqrt(dr(0)*dr(0) + dr(1)*dr(1));
                Bn(count,1) =  dr(0)/sqrt(dr(0)*dr(0) + dr(1)*dr(1));
            }
            else{
                Bn(count,0) =  dr(1)/sqrt(dr(0)*dr(0) + dr(1)*dr(1));
                Bn(count,1) = -dr(0)/sqrt(dr(0)*dr(0) + dr(1)*dr(1));
            }
        }
    }

    Area = MatrixXd::Zero(mesh.nElem, 1);
    for(int i=0; i<mesh.nElem; i++){
        Area(i) = tri_area(mesh.Node(mesh.Elem(i,0)-1,0), mesh.Node(mesh.Elem(i,0)-1,1), mesh.Node(mesh.Elem(i,1)-1,0),
         mesh.Node(mesh.Elem(i,1)-1,1), mesh.Node(mesh.Elem(i,2)-1,0), mesh.Node(mesh.Elem(i,2)-1,1));
    }
}