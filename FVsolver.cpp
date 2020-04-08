#include "headers.h"

void FVsolver(const mesh_t& mesh, const string& filename, const string& inputfile, MatrixXd U_initial, const string& file_prefix){

    string btype, outputR, outputF, outputU;
    double CFL, tmin, tmax, dt_out, nTime;
    int nIteration;

    read_param(inputfile, btype, CFL, outputR, outputF, outputU, tmin, tmax, dt_out, nIteration, nTime);
    MatrixXi I2E, B2E, IE, BE;
    MatrixXd In, Bn, Area;
    process_gri(filename, I2E, B2E, In, Bn, Area, IE, BE);

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

    int Nt = 1;
    double dt;
    MatrixXd U = U_initial;
    MatrixXd R = U*0.0;
    MatrixXd s = MatrixXd::Zero(1,R.cols()); 
    MatrixXd F(3,1);                         // stores numerical flux
    MatrixXd Force = MatrixXd::Zero(1,8);    // stores the forces exerted on 4 buildings
    double smax;
    ofstream out;
    double Tout = tmin;                      // time at which outputs states
    int counter = 0;                         // count the states output file
    for(int n=0; n<Nt; n++){

        R *= 0.0; s *= 0.0; Force *= 0.0;

        // loop over interior faces
        for(int ie=0; ie<I2E.rows(); ie++){
            int elemL = I2E(ie,0); 
            int elemR = I2E(ie,2);
            F = roe(U.col(elemL-1), U.col(elemR-1), In(ie, 0), In(ie, 1), smax);
            R.col(elemL-1) += F*Il(ie); R.col(elemR-1) -= F*Il(ie);
            s(elemL-1) += smax*Il(ie);  s(elemR-1) += smax*Il(ie);
        }        

        // loop over boundary faces
        for(int be=0; be<B2E.rows(); be++){
            int elemL = B2E(be,0);
            if(btype == "Full"){
                F = roe(U.col(elemL-1), U_initial.col(elemL-1), Bn(be, 0), Bn(be, 1), smax);
                R.col(elemL-1) += F*Bl(be);
                s(elemL-1) += smax*Bl(be);

                // calculate forces exerted on buildings
                if((outputF == "True") && (B2E(be,2) > 1)){
                    int nbuilding = B2E(be,2) - 2;
                    double h = (U(0,elemL-1) + U_initial(0,elemL-1))/2.0;
                    Force(2*nbuilding) += Bl(be)*rho*g*h*h*Bn(be,0)/2.0;
                    Force(2*nbuilding + 1) += Bl(be)*rho*g*h*h*Bn(be,1)/2.0;
                }
            }
            else if(btype == "Wall"){
                double h = U(0, elemL-1);
                F(0) = 0.0; F(1) = g*h*h*Bn(be, 0)/2.0; F(2) = g*h*h*Bn(be, 1)/2.0;
                smax = sqrt(g*h);
                R.col(elemL-1) += F*Bl(be);
                s(elemL-1) += smax*Bl(be);

                // calculate forces exerted on buildings
                if((outputF == "True") && (B2E(be,2) > 1)){
                    int nbuilding = B2E(be,2) - 2;
                    Force(2*nbuilding) += F(1)*rho*Bl(be);
                    Force(2*nbuilding + 1) += F(2)*rho*Bl(be);
                }
            }
            else{
                cout << "Unsupported boundary type" << endl;
                return;
            }
        }

        // determine time step (dt) and number of iteration (Nt)
        if(n == 0){
            MatrixXd dt_local(R.cols(),1);
            for(int i=0; i<R.cols(); i++){
                dt_local(i) = 2.0*Area(i)*CFL/s(i);
            }
            dt = dt_local.minCoeff();
            if(outputU == "True"){
                dt = dt_out/ceil(dt_out/dt);
                Nt = int(nTime/dt);
            }   
            else{
                Nt = nIteration;
            }       
            cout << "dt = " << dt << endl;  // output time step to terminal  
        }

        // Ouput (1)residual, (2)force exerted on buildings
        if(n == 0){
            if(outputR == "True"){
                out.open(file_prefix + "_Res.dat"); // open(create) a file
                out << R.cwiseAbs().sum() << endl;
                out.close();
            }
            if(outputF == "True"){
                out.open(file_prefix + "_Force.dat"); // open(create) a file    
                out << Force << endl;   // format: F1x, F1y, F2x, F2y, F3x, F3y, F4x, F4y
                out.close();
            }
        }
        else{
            if(outputR == "True"){
                out.open(file_prefix + "_Res.dat", ofstream::app);  // turn on append mode
                out << R.cwiseAbs().sum() << endl;
                out.close();
            }
            if(outputF == "True"){
                out.open(file_prefix + "_Force.dat", ofstream::app); // turn on append mode    
                out << Force << endl;
                out.close();
            }
        }

        // output states
        if((outputU == "True") && (fabs(Tout - dt*double(n)) < (1e-4)*dt)){
            if(fabs(Tout - tmax) > 0.1*dt_out){Tout += dt_out;}
            string outputname = file_prefix + "_State_" + to_string(counter) + ".dat";
            out.open(outputname);
            out << U << endl;
            out.close();
            counter += 1;
        }

        // forward Euler time stepping, loop over elements
        for(int i=0; i<R.cols(); i++){
            U.col(i) -= dt/Area(i)*R.col(i);
        }    
    }
}