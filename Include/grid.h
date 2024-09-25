#pragma once
#ifndef GRID_H
#define GRID_H

#include "interface.h"
#include "string"
#include "map"
#include "Matrix.h"
#include "Matrix_arma.h"
#include "Matrix_arma_sp.h"
#include "BTC.h"
#include "vtk.h"
#include "Position.h"
#include "gsl/gsl_rng.h"

using namespace std;

struct field_gen_params
{
    int max_correl_n = 10;
    double k_correlation_lenght_scale_x;
    double k_correlation_lenght_scale_y;
    int n_filled=0;
    CTimeSeries<double> inversecdf;
};

class CDistribution;
class TimeSeriesD;
class TimeSeriesSetD;
class CPathwaySet;
class CPathway;

class Grid : public Interface
{
public:
    Grid();
    bool CreateGrid(const map<string,string> &Arguments);
    static bool HasCommand(const string &cmd);
    vector<string> commands();
    static vector<string> Commands();
    FunctionOutPut Execute(const string &cmd, const map<string,string> &arguments);
    bool AssignKFieldToGrid(const map<string,string> &Arguments);
    bool RenormalizeKField(const map<string,string> &Arguments);
    bool SolveTransport(const map<string,string> &Arguments);
    void RenormalizeK(CDistribution *dist,int k=0);
    bool WriteKFieldToVTP(const map<string,string> &Arguments);
    bool WriteConcentrationToVTP(const map<string,string> &Arguments);
    bool SolveHydro(const map<string,string> &Arguments);
    bool WriteHydroSolutionToVTP(const map<string,string> &Arguments);
    CPathwaySet CreateTrajectories(const map<string,string> &Arguments);
    TimeSeriesD GetConcentrationBTCAtX(const map<string,string> &Arguments);
    TimeSeriesD GetMarginalVelocityDistribution(const map<string,string> &Arguments);
    TimeSeriesD GetVelocityDistributionAtXSection(const map<string,string> &Arguments);
    CTimeSeries<double> GetKValuesToTimeSeries(int k=0);
    void RemapKFieldBasedonMarginalDistribution(CDistribution *dist,int k=0);
    double MapToMarginalDistribution(const double &u, CDistribution *dist);
    bool WriteKFieldToVTP(const string &filename="surface.vtp", const double &z_factor=0.5, bool _log = false);
    bool WriteToVTP(const map<string,string> &Arguments);
    bool WriteHydroSolutionToVTP(const string &filename="solution.vtp", const double &z_factor=0.5, bool _log = false);
    bool SolveHydro(const double &leftboundary=0, const double &rightboundary=1);
    transportparameters TransportParameters;
    FunctionOutPut AssignStandardNormal2ndDegree(const map<string,string> &Arguments);
    bool SetStateVariable(const string &prop, const CVector_arma &v);
    CVector_arma GetStateVariable(const string &prop);
    CDistribution GetMarginalDistribution(const string &prop, int nbins);
    CDistribution GetMarginalDistribution(const map<string,string> &Arguments);
    bool SetStateVariable(int i, int j, const string &prop, const double &val);
    double GetStateVariable(int i, int j, const string &prop);
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
    bool SolveTransport(const double &t_end, const vector<double> &decay_coeff, const vector<double> &decay_order);
    void CreateTransportKMatrix(const double &dt, const double &D, const double &weight);
    CVector_arma CreateTransportRHS(int species_counter, const double &dt, const double &weight, const double &D, const double &decay_coefficient, const double &decay_order);
    vector<CMatrix> C;
    bool WriteConcentrationToVTP(int species_counter, const string &filename, const double &z_factor, bool _log, const vector<double> &t);
    bool WriteConcentrationToVTP(int species_counter, const string &filename, const double &z_factor, bool _log, const double &t);
    TimeSeriesD GetConcentrationBTCAtX(int species_counter, const double &x, const string &filename, const string &filename_d);
    double GetConcentrationAtX(int species_counter, const double &x, int timestep);
    vector<CPosition> InitializeTrajectories(int numpoints, const double &x_0);
    TimeSeriesD GetVelocityDistributionAtXSection(const double &x, Direction dir, int nbins, const double &smoothing_factor);
    TimeSeriesD GetVelocitiesAtXSection(const double &x, Direction dir);
    CVector GetVelocity(const CPosition &pp);
    CPathwaySet BuildTrajectories(const vector<CPosition> pts, const double &dx, const double &x_end, const double &tol, const double &diffusion);
    CPathway CreateSingleTrajectoryFixDx(const CPosition &pp, const double &dx0, const double &x_end, const double &D, const double &tol);
    gsl_rng *rng_ptr;
    TimeSeriesD GetMarginalVelocityDistribution(Direction dir, int nbins, const double &smoothing_factor);
    TimeSeriesSetD GetVelocityCorrelationsBasedOnRandomSamples(int nsamples, double dx0, double x_inc, bool magnitude);
    CPosition GetARandomPoint();
    CVector GetPairVelocity(const CPosition &pp, double dx0, double x_inc, bool magnitude);
    CVector GetVelocityAt(const CPosition &pp);

};
vector<ijval> GetClosestCells(vector<ijval> vec, int n);
#endif // GRID_H
