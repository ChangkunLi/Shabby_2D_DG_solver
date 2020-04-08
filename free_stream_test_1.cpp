#include "headers.h"

// Testing full state boundary condition

void free_strem_test_1(){
    string filename  = "city0.gri";
    string inputfile = "PARAM.in.full_state_test";

    mesh_t mesh;
    read_gri(filename, mesh);

    for(int p=0; p<3; p++){
        // Initialize state
        int N = (p+1)*(p+2)/2;
        MatrixXd U_initial = MatrixXd::Ones(3, N*mesh.nElem);
        U_initial.row(1) *= 0.453;
        U_initial.row(2) *= 0.769;
        Map<MatrixXd> u_initial(U_initial.data(), 3*N*mesh.nElem, 1);
        U_initial = u_initial;
        cout << "Free-stream full state boundary condition test (p=" << p << "):" << endl;
        dg(p ,mesh, filename, inputfile, U_initial, "fullstate_test_p" + to_string(p));
    }
}