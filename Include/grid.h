#ifndef GRID_H
#define GRID_H

#include "interface.h"
#include "Structs.h"
#include "string"
#include "map"
#include "Matrix.h"
#include "Matrix_arma.h"
#include "Matrix_arma_sp.h"
#include "BTC.h"
#include "Distribution.h"
#include "vtk.h"

using namespace std;

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

struct field_gen_params
{
    int max_correl_n = 10;
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
    bool AssignKFieldToGrid(map<string,string> Arguments);
    bool RenormalizeKField(map<string,string> Arguments);
    void RenormalizeK(CDistribution *dist,int k=0);
    bool WriteKFieldToVTP(const map<string,string> &Arguments);
    bool SolveHydro(const map<string,string> &Arguments);
    CTimeSeries<double> GetKValuesToTimeSeries(int k=0);
    void RemapKFieldBasedonMarginalDistribution(CDistribution *dist,int k=0);
    double MapToMarginalDistribution(const double &u, CDistribution *dist);
    bool WriteKFieldToVTP(const string &filename="surface.vtp", const double &z_factor=0.5, bool _log = false);
    bool SolveHydro(const double &leftboundary=0, const double &rightboundary=1);
    void SetProgressValue(const double &s);

private:
    geometrical_parameters GeometricParameters;
    vector<vector<distributed_property> > p;
    double max_K();
    double min_K();
    double max_vx();
    double min_vx();
    double max_vy();
    double min_vy();
    double min_v_x = 0;
    double max_v_x=0;
    void AssignNewK(int i, int j,field_gen_params *FieldGeneratorParameters);
    void Clear();
    correl_mat_vec GetCorrellMatrixVec(int i, int j, field_gen_params *FieldGeneratorParameters);
    vector<ijval> GetClosestDeteminedCells(int i, int j, int n, field_gen_params *FieldGeneratorParameters);
    CMatrix_arma_sp CreateStiffnessMatrixHydro();
    int get_cell_no(int i, int j);
    CMatrix H;
    CMatrix vx;
    CMatrix vy;
    CVector CreateRHSHydro(const double &leftboundary, const double &rightboundary);
    CVector_arma CreateRHSHydro_ARMA(const double &leftboundary, const double &rightboundary);

};
vector<ijval> GetClosestCells(vector<ijval> vec, int n);
void SetProgressValue(double s);
#endif // GRID_H
