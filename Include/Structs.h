#ifndef STRUCTS_H
#define STRUCTS_H

#include "Vector.h"
#include "Concentrations.h"
//#include "BTC.h"
#include "Matrix_arma.h"
#include "Matrix_arma_sp.h"
#include "Matrix.h"

class Interface;

enum class object_type {none, grid, distribution, timeseries, pathwayset, timeseriesset};

struct FunctionOutPut
{
    bool success;
    Interface *output = nullptr;
};

enum class Direction {x, y};

struct geometrical_parameters
{
    int nx, ny;
    double dx, dy;
    int nx_data, ny_data;
};

struct distributed_property
{
    CVector K;
    CVector K_gauss;
    CVector V;
    double Vtx;
    double Vbx;
    double Vfy;
    double Vby;
    double u = 0;
    double omega = 0;
    bool k_det = false;
    Concentrations C;
};

struct correl_mat_vec
{
#ifdef Use_Armadillo
    CMatrix_arma M_22;
    CVector_arma V_21;
    CVector_arma V_RHS;
#else
    CMatrix M_22;
    CVector V_21;
    CVector V_RHS;
#endif //  arma
};

struct ijval
{
    int i;
    int j;
    double val;
};



struct transportparameters
{
    CMatrix_arma_sp Kv;
    CMatrix_arma_sp KD;
    CMatrix_arma_sp Kt;
    double dt;
    double time_weight=1;
    double D=0;
    int numberofspecies=1;
    vector<double> leftboundary_C;

};
#endif // STRUCTS_H
