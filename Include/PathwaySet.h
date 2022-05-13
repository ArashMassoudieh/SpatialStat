#pragma once
#include "Pathway.h"
#include <vector>
#include <string>
#include "Copula.h"
#include "interface.h"


using namespace std;
class CPathwaySet : public Interface
{
public:
    bool weighted;
    CPathwaySet();
    CPathwaySet(int n);
    CPathwaySet(unsigned int number_of_paths);
    CPathwaySet(const CPathwaySet &P);
    CPathwaySet &operator=(const CPathwaySet &P);
    ~CPathwaySet();
    bool getfromMODflowfile(const string &filename);
    bool getfromShermanfile(const string &filename, const string &reactionfilename="", int columnnumber=0);
    bool getfromShermanfile_v(const string &filename);
    vector<CPathway> paths;
    void write(string filename, int interval=1);
    void append(const CPathway& P, double weight=1);
    int max_num_points();
    int n() {return paths.size();}
    void create_ou_paths(int n, CDistribution *dist, double x_min, double x_max, double kappa, double dx, double weight=1);
    void create_copula_paths(int n, CDistribution * dist, double x_min, double x_max, double epsilon, double r, double dx, double weight=1);
    void create_copula_paths(int n, CDistribution * dist, double x_min, double x_max, double epsilon, CCopula *copula, double dx, double weight=1);
    void write_vtk(vtkSmartPointer<vtkPolyDataMapper>, string filename);
    vtkSmartPointer<vtkPolyDataMapper> pathways_vtk_pdt_vtp(double z_factor=1, double offset=0);
    CPathway snapshotattime(double t);
    CPathway snapshotatlocation(double x);
    void make_uniform_at_x(double dx);
    void make_uniform_at_t(double dt);
    CTimeSeries<double> sample_velocities();
    CPosition get_pair_v_pos(int increment, int num_seq=2);
    CTimeSeriesSet<double> get_pair_v(int increment, int n, int num_seq=2);
    CTimeSeries<double> get_BTC(double x, int n_bins, bool vel_inv_weighted = true, double smoothing_factor=0);
    CTimeSeries<double> get_BTC_points(double x, bool vel_inv_weighted = true);
    void set_progress_value(double s);
    void set_progress_value(string s);
    bool AssignVelocities();
    void show_in_window(string s);
    FunctionOutPut Execute(const string &cmd, const map<string,string> &arguments);
    static bool HasCommand(const string &cmd);
    static vector<string> Commands();

};

