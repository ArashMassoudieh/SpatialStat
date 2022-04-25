#ifndef GRID_H
#define GRID_H

#include "interface.h"
#include "Structs.h"
#include "string"
#include "map"

using namespace std;

class Grid : public Interface
{
public:
    Grid();
    bool CreateGrid(const map<string,string> &Arguments);
private:
    geometrical_parameters GeometricParameters;
    vector<vector<distributed_property> > p;
};

#endif // GRID_H
