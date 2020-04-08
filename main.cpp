#include "headers.h"

int main(){
    // test implementation of roe flux
    flux_test();

    // free-stream / free-stream preservation test
    free_strem_test_1();
    free_strem_test_2();

    vector<string> Citys;
    Citys.push_back("city0.gri");
    Citys.push_back("city1.gri");
    for(int icity=0; icity<Citys.size(); icity++){
        string filename  = Citys[icity];
        string inputfile = "PARAM.in";

        mesh_t mesh;
        read_gri(filename, mesh);

        for(int p=0; p<=3; p++){
            // Initialize state
            int N = (p+1)*(p+2)/2;
            MatrixXd U_initial = MatrixXd::Zero(3, N*mesh.nElem);

            ArrayXd xi = ArrayXd::LinSpaced(p+1,0,1); 
            ArrayXd eta = xi;
            for(int i=0; i<mesh.nElem; i++){
                double x1,y1,x2,y2,x3,y3;
                x1 = mesh.Node(mesh.Elem(i,0)-1,0);   y1 = mesh.Node(mesh.Elem(i,0)-1,1);
                x2 = mesh.Node(mesh.Elem(i,1)-1,0);   y2 = mesh.Node(mesh.Elem(i,1)-1,1);
                x3 = mesh.Node(mesh.Elem(i,2)-1,0);   y3 = mesh.Node(mesh.Elem(i,2)-1,1);
                int k = 0;
                double x,y;
                for(int s=0; s<=p; s++){
                    for(int r=0; r<=(p-s); r++){
                        if(p==0){
                            x = (x1+x2+x3)/3.0;
                            y = (y1+y2+y3)/3.0;
                        }
                        else{
                            x = x1*(1.0 - xi(r) - eta(s)) + x2*xi(r) + x3*eta(s);
                            y = y1*(1.0 - xi(r) - eta(s)) + y2*xi(r) + y3*eta(s);
                        }
                        U_initial(0, i*N+k) = 1.0 + 0.3*exp(-50.0*((x-1.5)*(x-1.5) + (y-0.7)*(y-0.7)));
                        k += 1;
                    }
                }
            }
            Map<MatrixXd> u_initial(U_initial.data(), 3*N*mesh.nElem, 1);
            U_initial = u_initial;

            cout << filename << " p=" << to_string(p) << " starts" << endl;
            dg(p, mesh, filename, inputfile, U_initial, "city" + to_string(icity) + "_p" + to_string(p));
            cout << filename << " p=" << to_string(p) << " ends" << endl;
        }
    }
    return 0;
}
