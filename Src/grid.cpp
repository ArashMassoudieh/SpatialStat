#include "grid.h"

Grid::Grid():Interface()
{

}

bool Grid::CreateGrid(const map<string,string> &Arguments)
{
    GeometricParameters.nx = atoi(Arguments.at("nx").c_str());
    GeometricParameters.nx = atoi(Arguments.at("ny").c_str());
    GeometricParameters.nx = atof(Arguments.at("dx").c_str());
    GeometricParameters.ny = atof(Arguments.at("dy").c_str());
    p.resize(GeometricParameters.nx);
    for (int i = 0; i < GeometricParameters.nx; i++)
    {
        p[i].resize(GeometricParameters.ny);
        for (int j = 0; j < GeometricParameters.ny; j++)
        {
            p[i][j].V = CVector(2);
            p[i][j].K = CVector(2);
            p[i][j].K_gauss = CVector(2);
        }
    }

}

