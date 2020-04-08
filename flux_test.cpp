#include "headers.h"

void flux_test(){
    MatrixXd U(3,2);
    U << 0.5, 1.5,
         0.0, 0.0,
         0.0, 0.0;
    double nx,ny,smax;
    nx = 1.0;   ny = 0.0;
    MatrixXd UL = U.col(0)/g;
    MatrixXd UR = U.col(1)/g;
    MatrixXd F = roe(UL,UR,nx,ny,smax);
    cout << "Normal vector pointing towards +x direction:" << endl;
    cout << "Roe flux multiplied by g is \n" << F*g << endl;
    nx = 0.0;   ny = 1.0;
    F = roe(UL,UR,nx,ny,smax);
    cout << "Normal vector pointing towards +y direction:" << endl;
    cout << "Roe flux multiplied by g is \n" << F*g << endl;    
}