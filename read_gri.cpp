#include "headers.h"

void read_gri(const string& filename, mesh_t& mesh){
    ifstream file;
    file.open(filename);

    // read in nodes
    string line;
    getline(file,line);
    istringstream strm(line);
    strm >> mesh.nNode;
    int nelemtot;
    strm >> nelemtot;
    strm >> mesh.Dim;
    strm.clear();
    mesh.Node = MatrixXd::Ones(mesh.nNode,mesh.Dim);
    for(int i=0 ; i<mesh.nNode ; i++){
        getline(file,line);
        strm.str(line);
        for(int j=0 ; j<mesh.Dim ; j++){
            strm >> mesh.Node(i,j);
        }
        strm.clear();
    }

    // read boundary info
    boundary_t B;
    getline(file,line);
    strm.str(line);
    strm >> B.nbfgrp;
    strm.clear();
    B.nbface = MatrixXi::Ones(B.nbfgrp,1);
    B.nnode  = MatrixXi::Ones(B.nbfgrp,1);
    for(int ibfgrp=0 ; ibfgrp<B.nbfgrp ; ibfgrp++){
        getline(file,line);
        strm.str(line);
        strm >> B.nbface(ibfgrp);
        strm >> B.nnode(ibfgrp);
        string boundary_title;
        strm >> boundary_title;
        strm.clear();
        B.title.push_back(boundary_title);
        MatrixXi N = MatrixXi::Ones(B.nbface(ibfgrp),B.nnode(ibfgrp));
        for(int ibface=0; ibface<B.nbface(ibfgrp) ; ibface++){
            getline(file,line);
            strm.str(line);
            for(int inode=0 ; inode<B.nnode(ibfgrp) ; inode++){
                strm >> N(ibface,inode);
            }
            strm.clear();
        }
        B.nodes.push_back(N);
    }
    mesh.B = B;

    // read in elements
    getline(file,line);
    strm.str(line);
    strm >> mesh.nElem;
    strm >> mesh.QOrder;
    strm >> mesh.QBasis;
    strm.clear();
    int nnode;
    if(mesh.QBasis == "TriLagrange"){
        nnode = (mesh.QOrder + 1)*(mesh.QOrder + 2)/2;
    }
    else{
        cout << "element type not understood" << endl;
        return;
    }

    mesh.Elem = MatrixXi::Ones(mesh.nElem,nnode);
    for(int elem=0; elem<mesh.nElem; elem++){
        getline(file,line);
        strm.str(line);
        for(int inode=0; inode<nnode ; inode++){
            strm >> mesh.Elem(elem,inode);
        }
        strm.clear();
    }
    // close file
    file.close();
}