#include "headers.h"

void quad1d(int order, int& nq, MatrixXd& xyq, MatrixXd& wq){
    if(order == 1){
        nq = 1;
        MatrixXd xyq_1(nq,1);
        xyq_1 << 0.500000000000000;
        xyq = xyq_1;
        MatrixXd wq_1(nq,1);
        wq_1 << 1.000000000000000;
        wq = wq_1;
    }
    else if(order == 3){
        nq = 2;
        MatrixXd xyq_1(nq,1);
        xyq_1 << 0.211324865405187, 0.788675134594813;
        xyq = xyq_1;
        MatrixXd wq_1(nq,1);
        wq_1 << 0.500000000000000, 0.500000000000000;
        wq = wq_1;     
    }
    else if(order == 5){
        nq = 3;
        MatrixXd xyq_1(nq,1);
        xyq_1 << 0.112701665379258, 0.500000000000000, 0.887298334620742;
        xyq = xyq_1;
        MatrixXd wq_1(nq,1);
        wq_1 << 0.277777777777778, 0.444444444444444, 0.277777777777778;
        wq = wq_1;     
    }
    else if(order == 7){
        nq = 4;
        MatrixXd xyq_1(nq,1);
        xyq_1 << 0.069431844202974, 0.330009478207572, 0.669990521792428, 0.930568155797026;
        xyq = xyq_1;
        MatrixXd wq_1(nq,1);
        wq_1 << 0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727;
        wq = wq_1;     
    }
    else{
        cout << "Unsupported 1D quadrature order" << endl;
    }
}