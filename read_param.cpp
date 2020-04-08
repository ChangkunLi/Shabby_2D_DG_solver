#include "headers.h"

void read_param(const string& inputfile, string& btype, double& CFL,
        string& outputR, string& outputF, string& outputU, double& tmin,
        double& tmax, double& dt_out, int& nIteration, double& nTime){
    
    ifstream file;
    istringstream strm;
    string line;    
    file.open(inputfile);

    // # Boundary type (Full or Wall)
    getline(file, line); // take away header
    getline(file, line);
    strm.str(line);
    strm >> btype;
    strm.clear();
    getline(file, line); // remove empty line

    // # CFL
    getline(file, line); // take away header
    getline(file, line);
    strm.str(line);
    strm >> CFL;
    strm.clear();
    getline(file, line); // remove empty line    

    // # Residual Output
    getline(file, line); // take away header
    getline(file, line);
    strm.str(line);
    strm >> outputR;
    strm.clear();
    getline(file, line); // remove empty line  

     // # Residual Output
    getline(file, line); // take away header
    getline(file, line);
    strm.str(line);
    strm >> outputF;
    strm.clear();
    getline(file, line); // remove empty line      

    // # Residual Output
    getline(file, line); // take away header
    getline(file, line);
    strm.str(line);
    strm >> outputU;
    strm.clear();
    getline(file, line);
    strm.str(line);
    strm >> tmin;
    strm.clear();
    getline(file, line);
    strm.str(line);
    strm >> tmax;
    strm.clear();
    getline(file, line);
    strm.str(line);
    strm >> dt_out;
    strm.clear();        
    getline(file, line); // remove empty line   

    // # Simulation Time
    getline(file, line); // take away header
    getline(file, line);
    strm.str(line);
    strm >> nIteration;
    strm.clear();    
    getline(file, line);
    strm.str(line);
    strm >> nTime;
    strm.clear();      

    file.close();
}