#include "headers.h"

void edgehash(const MatrixXi& E2N, MatrixXi& IE, MatrixXi & BE){
    // This function identifies interior and boundary edges, and their
    // connectivities, in a triangular mesh given an element-to-node array.

    // INPUT : E2N = [nelem x 3] array mapping triangles to nodes
    // OUTPUT: IE = [niedge x 4] array giving (n1, n2, elem1, elem2)
    //               information for each interior edge
    //          BE = [nbedge x 3] array giving (n1, n2, elem)
    //               information for each boundary edge

    int nelem = E2N.rows();                         // number of elements
    int nnode = E2N.maxCoeff();     // number of nodes
    SparseMatrix<int> H(nnode,nnode);               // Create a hash list to identify edges
    IE = MatrixXi::Zero(int(nelem*3.0/2.0), 4);     // (over) allocate interior edge array
    int niedge = 0;                                 // number of interior edges (running total)

    // Loop over elements and identify all edges
    for(int elem=1; elem<=nelem; elem++){
        MatrixXi nv = E2N.row(elem-1);
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
                niedge += 1;
                IE(niedge-1,0) = n1;
                IE(niedge-1,1) = n2;
                IE(niedge-1,2) = oldelem;
                IE(niedge-1,3) = elem;
                H.coeffRef(n1-1,n2-1) = -1; H.coeffRef(n2-1,n1-1) = -1;
            }
        }
    }

    IE = IE.block(0,0,niedge,4); // clip IE

    // find boundary edges
    int nbedge = nelem*3 - niedge*2;
    BE = MatrixXi::Zero(nbedge,3);
    int ibedge = 0;             // counting coundary edge
    for (int k=0; k<H.outerSize(); ++k)
        for (SparseMatrix<int>::InnerIterator it(H,k); it; ++it)
        {
            if(it.col()>it.row() && it.value()>0){                   
                BE(ibedge,0) = it.row() + 1;   // row index
                BE(ibedge,1) = it.col() + 1;   // col index (here it is equal to k)
                BE(ibedge,2) = it.value();
                ibedge += 1; 
            }
        }
}