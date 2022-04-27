#include "grid.h"
#include "Utilities.h"

vector<string> list_of_commands = vector<string>({"CreateGrid"});

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
    return true;
}

vector<string> Grid::commands()
{
    return Commands();
}

vector<string> Grid::Commands()
{
    return list_of_commands;
}

bool Grid::HasCommand(const string &cmd)
{
    if (aquiutils::lookup(Commands(),cmd)!=-1)
        return true;
    else
        return false;
}
