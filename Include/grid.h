#ifndef GRID_H
#define GRID_H

#include "interface.h"
#include "Structs.h"
#include "string"
#include "map"
#include "Matrix.h"
#include "Matrix_arma.h"
#include "BTC.h"
#include "Distribution.h"

using namespace std;

struct correl_mat_vec
{
#ifdef  arma
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

struct field_gen_params
{
    int max_correl_n;
    double k_correlation_lenght_scale_x;
    double k_correlation_lenght_scale_y;
    int n_filled=0;
    CTimeSeries<double> inversecdf;
};

class Grid : public Interface
{
public:
    Grid();
    bool CreateGrid(const map<string,string> &Arguments);
    static bool HasCommand(const string &cmd);
    vector<string> commands();
    static vector<string> Commands();
    static vector<string> list_of_commands;
    bool Execute(const string &cmd, const map<string,string> &arguments);
    void CreateRandomKField(map<string,string> Arguments, CDistribution *dist);
private:
    geometrical_parameters GeometricParameters;
    vector<vector<distributed_property> > p;

    void AssignNewK(int i, int j,field_gen_params *FieldGeneratorParameters);
    void Clear();
    correl_mat_vec GetCorrellMatrixVec(int i, int j, field_gen_params *FieldGeneratorParameters);
    vector<ijval> GetClosestDeteminedCells(int i, int j, int n, field_gen_params *FieldGeneratorParameters);
};
vector<ijval> GetClosestCells(vector<ijval> vec, int n);
void SetProgressValue(double s);
#endif // GRID_H
