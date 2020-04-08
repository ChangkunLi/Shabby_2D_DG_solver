#include "headers.h"

MatrixXd roe(const MatrixXd& UL, const MatrixXd& UR, double nx, double ny, double& smax){
    // Left flux
    MatrixXd FL(3,1);
    double hL = UL(0);
    double uL = UL(1)/UL(0);
    double vL = UL(2)/UL(0);
    FL(0) = nx*UL(1) + ny*UL(2);
    FL(1) = (nx*UL(1)*UL(1) + ny*UL(1)*UL(2))/UL(0) + nx*g*UL(0)*UL(0)/2.0;
    FL(2) = (nx*UL(1)*UL(2) + ny*UL(2)*UL(2))/UL(0) + ny*g*UL(0)*UL(0)/2.0;

    // Right flux
    MatrixXd FR(3,1);
    double hR = UR(0);
    double uR = UR(1)/UR(0);
    double vR = UR(2)/UR(0);
    FR(0) = nx*UR(1) + ny*UR(2);
    FR(1) = (nx*UR(1)*UR(1) + ny*UR(1)*UR(2))/UR(0) + nx*g*UR(0)*UR(0)/2.0;
    FR(2) = (nx*UR(1)*UR(2) + ny*UR(2)*UR(2))/UR(0) + ny*g*UR(0)*UR(0)/2.0;

    // define Roe average state
    double h = (hR + hL)/2.0;
    double u = (sqrt(hR)*uR + sqrt(hL)*uL)/(sqrt(hR) + sqrt(hL));
    double v = (sqrt(hR)*vR + sqrt(hL)*vL)/(sqrt(hR) + sqrt(hL));

    // define normal velocity and wave speed
    double vn = nx*u + ny*v;
    double c  = sqrt(g*h);

    MatrixXd l(3,1);
    l(0) = fabs(vn + c);
    l(1) = fabs(vn - c);
    l(2) = fabs(vn);

    // entropy fix
    double epsilon = 0.01*c;
    for(int i=0 ; i<3 ; i++){
        if(l(i) < epsilon){
            l(i) = (epsilon*epsilon + l(i)*l(i))/(2.0*epsilon);
        }
    }

    // return maximum wave speed
    smax = l.maxCoeff();

    // define right eigenvectors
    MatrixXd r1(3,1), r2(3,1), r3(3,1); 
    r1(0) = 1.0; r1(1) = u + nx*c; r1(2) = v + ny*c;
    r2(0) = 1.0; r2(1) = u - nx*c; r2(2) = v - ny*c;
    r3(0) = 0.0; r3(1) = -ny*c;    r3(2) = nx*c;

    // Corfficients
    double dh  = UR(0) - UL(0);
    double dhu = UR(1) - UL(1);
    double dhv = UR(2) - UL(2);
    double C1  = ((c - vn)*dh + nx*dhu + ny*dhv)/(2.0*c);
    double C2  = ((c + vn)*dh - nx*dhu - ny*dhv)/(2.0*c);
    double C3  = ((ny*u - nx*v)*dh - ny*dhu + nx*dhv)/c;

    return (FL + FR)/2.0 - (C1*l(0)*r1 + C2*l(1)*r2 + C3*l(2)*r3)/2.0;
}