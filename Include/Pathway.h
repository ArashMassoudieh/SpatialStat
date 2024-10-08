#pragma once
#include <vector>
#include "Position.h"
#include "vtk.h"
#include <string>
#include <timeseriesd.h>
#include <BTCSet.h>
#include "Copula.h"

using namespace std;

class CDistribution;

class CPathway
{
public:
	CPathway();
	CPathway(const CPathway & P);
	void create_range_x(double x_min, double x_max, double dx);
	CPathway& operator=(const CPathway & P);
	~CPathway();
	vector<CPosition> positions;
	void append(CPosition pos);
	CVector get_velocity_at_x(double x);
	CPosition get_position_at_x(double x);
	CPosition get_position_at_t(double t);
	CVector get_velocity_at_t(double t);
	CPathway make_uniform_x(double dx);
	CPathway make_uniform_t(double dt);
	int size() {return positions.size();}
	double max_x();
	double min_x();
	double max_t();
	double min_t();
	double maxx=-0.1234567;
	double minx=-0.1234567;
	double maxt = -0.1234567;
	double mint = -0.1234567;
	void write(string filename);
	void create_ou(CDistribution *dist, double x_min, double x_max, double kappa, double dx, double weight=1);
	void create_copula(CDistribution *dist, double x_min, double x_max, double epsilon, double r, double dx, double weight=1);
	void create_copula(CDistribution *dist, double x_min, double x_max, double epsilon, CCopula *copula, double dx, double _weight=1);
	vtkSmartPointer<vtkPolyData> pathway_vtk_pdt_vtp(double z_factor = 1, double offset = 0);
    CTimeSeriesSet<double> get_distribution(bool _log, int n_bins);
    TimeSeriesD get_distribution(string var, bool _log = false, int nbins = 100);
    TimeSeriesD get_distribution(string var, string weight_var, bool _log, int n_bins=100);
	vector<double> minmax(string var);
	bool uniform = false;
	double get_cross_time(double x);
    vector<double> get_cross_time_vx(double x);
    double weight = 1;
    vtkSmartPointer<vtkPolyData> TovtkPolyData(const double &z_factor, const double &offset, bool _log, bool _color);

};

