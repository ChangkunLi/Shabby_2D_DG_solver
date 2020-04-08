#include "headers.h"

void dg(int p, const mesh_t& mesh, const string& filename, const string& inputfile, MatrixXd U_initial, const string& file_prefix){

    string btype, outputR, outputF, outputU;
    double CFL, tmin, tmax, dt_out, nTime;
    int nIteration;

    read_param(inputfile, btype, CFL, outputR, outputF, outputU, tmin, tmax, dt_out, nIteration, nTime);
    MatrixXi I2E, B2E, IE, BE;
    MatrixXd In, Bn, Area;
    process_gri(filename, I2E, B2E, In, Bn, Area, IE, BE);

    // precomputation for DG
    MatrixXd C, M_proto, Phi, dxPhi, dyPhi, wq_2d, wq_1d;
    vector<MatrixXd> Phi_edge_CC, Phi_edge_c;
    precompute(p, C, M_proto, Phi, dxPhi, dyPhi, Phi_edge_CC, Phi_edge_c, wq_2d, wq_1d);   
    Res_data_t res_data;
    res_data.C = C; res_data.dxPhi = dxPhi; res_data.dyPhi = dyPhi; res_data.Phi = Phi;
    res_data.Phi_edge_c = Phi_edge_c; res_data.Phi_edge_CC = Phi_edge_CC;
    res_data.wq_1d = wq_1d; res_data.wq_2d = wq_2d;

    vector<MatrixXd> Pseudo_invJ; // each element is det(J)J^(-1)
    for(int i=0; i<mesh.nElem; i++){
        MatrixXd mat(2,2);
        double x1,x2,x3,y1,y2,y3;
        x1 = mesh.Node(mesh.Elem(i,0)-1,0);   y1 = mesh.Node(mesh.Elem(i,0)-1,1);
        x2 = mesh.Node(mesh.Elem(i,1)-1,0);   y2 = mesh.Node(mesh.Elem(i,1)-1,1);
        x3 = mesh.Node(mesh.Elem(i,2)-1,0);   y3 = mesh.Node(mesh.Elem(i,2)-1,1);
        mat(0,0) = y3-y1;   mat(0,1) = x1-x3;
        mat(1,0) = y1-y2;   mat(1,1) = x2-x1;
        Pseudo_invJ.push_back(mat);
    }
    res_data.Pseudo_invJ = Pseudo_invJ;

    // precompute length of every face (interior & boundary)
    MatrixXd Il(In.rows(), 1);
    MatrixXd Bl(Bn.rows(), 1);
    for(int i=0; i<Il.rows(); i++){
        int n1 = IE(i,0);
        int n2 = IE(i,1);
        Il(i) = sqrt( (mesh.Node(n1-1, 0) - mesh.Node(n2-1, 0))*(mesh.Node(n1-1, 0) - mesh.Node(n2-1, 0))
             + (mesh.Node(n1-1, 1) - mesh.Node(n2-1, 1))*(mesh.Node(n1-1, 1) - mesh.Node(n2-1, 1)) );
    }
    for(int i=0; i<Bl.rows(); i++){
        int elem = B2E(i,0);
        int n1 = mesh.Elem(elem-1, (B2E(i,1) % 3));
        int n2 = mesh.Elem(elem-1, ((B2E(i,1)+1) % 3));
        Bl(i) = sqrt( (mesh.Node(n1-1, 0) - mesh.Node(n2-1, 0))*(mesh.Node(n1-1, 0) - mesh.Node(n2-1, 0))
             + (mesh.Node(n1-1, 1) - mesh.Node(n2-1, 1))*(mesh.Node(n1-1, 1) - mesh.Node(n2-1, 1)) );
    }    
    res_data.I2E = I2E; res_data.B2E = B2E;
    res_data.In  = In;  res_data.Bn  = Bn;
    res_data.Il  = Il;  res_data.Bl  = Bl;

    int Nt = 1;
    double dt;
    MatrixXd U = U_initial;
    MatrixXd R = U*0.0;
    ofstream out;
    double Tout = tmin;                      // time at which outputs states
    int counter = 0;                         // count the states output file

    // determine time step (dt) and number of iteration (Nt)

    dt = sqrt(Area.minCoeff()) * (CFL/double(2*p+1))/sqrt(g*U.maxCoeff());
    if(outputU == "True"){
        dt = dt_out/ceil(dt_out/dt);
        Nt = int(nTime/dt);
    }   
    else{
        Nt = nIteration;
    }       
    cout << "dt = " << dt << endl;  // output time step to terminal   

    // Time marching
    for(int n=0; n<Nt; n++){
        // Ouput (1)residual, (2)force exerted on buildings
        if(n == 0){
            if(outputR == "True"){
                R = Res(U, p, 3, btype, res_data);
                out.open(file_prefix + "_Res.dat"); // open(create) a file
                out << R.cwiseAbs().sum() << endl;
                out.close();
            }
            if(outputF == "True"){
                out.open(file_prefix + "_Force.dat"); // open(create) a file    
                out << Force(U, p, 3, res_data) << endl;   // format: F1x, F1y, F2x, F2y, F3x, F3y, F4x, F4y
                out.close();
            }
        }
        else{
            if(outputR == "True"){
                R = Res(U, p, 3, btype, res_data);
                out.open(file_prefix + "_Res.dat", ofstream::app);  // turn on append mode
                out << R.cwiseAbs().sum() << endl;
                out.close();
            }
            if(outputF == "True"){
                out.open(file_prefix + "_Force.dat", ofstream::app); // turn on append mode    
                out << Force(U, p, 3, res_data) << endl;
                out.close();
            }
        }

        // output states
        if((outputU == "True") && (fabs(Tout - dt*double(n)) < (1e-4)*dt)){
            if(fabs(Tout - tmax) > 0.1*dt_out){Tout += dt_out;}
            string outputname = file_prefix + "_State_" + to_string(counter) + ".dat";
            out.open(outputname);
            out << Submesh(U, p, 3, res_data, mesh, 2*p) << endl;
            out.close();
            counter += 1;
        }

        // RK4
        R = Res(U, p, 3, btype, res_data);
        MatrixXd f0 = -Inv_M(M_proto, R, Area, p, 3);
        R = Res(U+f0*dt/2.0, p, 3, btype, res_data);
        MatrixXd f1 = -Inv_M(M_proto, R, Area, p, 3);
        R = Res(U+f1*dt/2.0, p, 3, btype, res_data);
        MatrixXd f2 = -Inv_M(M_proto, R, Area, p, 3);
        R = Res(U+f2*dt, p, 3, btype, res_data);
        MatrixXd f3 = -Inv_M(M_proto, R, Area, p, 3);
        U += dt/6.0*(f0 + 2.0*f1 + 2.0*f2 + f3);
    }
}