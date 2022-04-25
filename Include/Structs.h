#ifndef STRUCTS_H
#define STRUCTS_H

#endif // STRUCTS_H
#include "Vector.h"

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
    //Concentrations C;
};
