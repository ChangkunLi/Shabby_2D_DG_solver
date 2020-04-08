#include "iostream"
#include "fstream"
#include "sstream"
#include "string"
#include "vector"
#include "cmath"
#include "Eigen/Dense"
#include "Eigen/Sparse"

using namespace Eigen;
using namespace std;

#define g 9.8
#define rho 1000.0

struct boundary_t {
    int nbfgrp;
    MatrixXi nbface;
    MatrixXi nnode;
    vector<string> title;
    vector<MatrixXi> nodes;
};

struct mesh_t {
    int nNode;
    int Dim;
    MatrixXd Node;
    boundary_t B;
    int nElem;
    string QBasis;
    int QOrder;
    MatrixXi Elem;
};

struct Res_data_t {
    MatrixXd C;
    MatrixXd Phi;
    MatrixXd dxPhi;
    MatrixXd dyPhi;
    MatrixXd wq_2d;
    MatrixXd wq_1d;
    vector<MatrixXd> Phi_edge_CC;
    vector<MatrixXd> Phi_edge_c;
    vector<MatrixXd> Pseudo_invJ;
    MatrixXi I2E;
    MatrixXi B2E;
    MatrixXd In;
    MatrixXd Bn;
    MatrixXd Il;
    MatrixXd Bl;
};

MatrixXd roe(const MatrixXd& UL, const MatrixXd& UR, double nx, double ny, double& smax);
void flux_test();
void read_gri(const string& filename, mesh_t& mesh);
void edgehash(const MatrixXi& E2N, MatrixXi& IE, MatrixXi & BE);
double tri_area(double dX0, double dY0, double dX1, double dY1, double dX2, double dY2);
void process_gri(const string& filename, MatrixXi& I2E, MatrixXi& B2E, MatrixXd& In, MatrixXd& Bn, MatrixXd& Area, MatrixXi& IE, MatrixXi& BE);
void read_param(const string& inputfile, string& btype, double& CFL,
        string& outputR, string& outputF, string& outputU, double& tmin,
        double& tmax, double& dt_out, int& nIteration, double& nTime);
void FVsolver(const mesh_t& mesh, const string& filename, const string& inputfile, MatrixXd U_initial, const string& file_prefix);     
void free_strem_test_1();   
void free_strem_test_2();
MatrixXd TriLagrange2D(int p);
MatrixXd Inv_M(const MatrixXd& M_proto, const MatrixXd& R, const MatrixXd& Area, int p, int s);
MatrixXd MassMatrix2D(const MatrixXd& C, int p, MatrixXd& Phi);
void quad2d(int order, int& nq, MatrixXd& xyq, MatrixXd& wq);
void quad1d(int order, int& nq, MatrixXd& xyq, MatrixXd& wq);
void precompute(int p, MatrixXd& C, MatrixXd& M_proto, MatrixXd& Phi, MatrixXd& dxPhi, MatrixXd& dyPhi, vector<MatrixXd>& Phi_edge_CC, vector<MatrixXd>& Phi_edge_c, MatrixXd& wq_2d, MatrixXd& wq_1d);
void dg(int p, const mesh_t& mesh, const string& filename, const string& inputfile, MatrixXd U_initial, const string& file_prefix);
MatrixXd Res(const MatrixXd& U, int p, int s, const string& btype, const Res_data_t& res_data);
MatrixXd Force(const MatrixXd& U, int p, int s, const Res_data_t& res_data);
MatrixXd Submesh(const MatrixXd& U, int p, int s, const Res_data_t& res_data, const mesh_t& mesh, int nplot);