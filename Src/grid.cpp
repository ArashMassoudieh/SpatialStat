#include "grid.h"
#include "Utilities.h"
#include "NormalDist.h"


vector<string> Grid::list_of_commands = vector<string>({"CreateGrid"});

Grid::Grid():Interface()
{

}

bool Grid::CreateGrid(const map<string,string> &Arguments)
{
    GeometricParameters.nx = atoi(Arguments.at("nx").c_str());
    GeometricParameters.ny = atoi(Arguments.at("ny").c_str());
    GeometricParameters.dx = atof(Arguments.at("dx").c_str());
    GeometricParameters.dy = atof(Arguments.at("dy").c_str());
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
    //return vector<string>();
    return list_of_commands;
}

bool Grid::HasCommand(const string &cmd)
{
    if (aquiutils::lookup(Commands(),cmd)!=-1)
        return true;
    else
        return false;
}

bool Grid::Execute(const string &cmd, const map<string,string> &arguments)
{
    if (cmd=="CreateGrid")
        return CreateGrid(arguments);
    return false;
}

void Grid::CreateRandomKField(map<string,string> Arguments, CDistribution *dist)
{
    Clear();
    int n_filled = 0;
    field_gen_params Field_Generator_Parameters;
    Field_Generator_Parameters.inversecdf = dist->InverseCumulative;
    srand(time(NULL));
    while (n_filled<GeometricParameters.nx*GeometricParameters.ny)
    {
        if (n_filled%100==0)
        SetProgressValue(double(n_filled) / double(GeometricParameters.nx * GeometricParameters.nx));
        int i = unitrandom()*(GeometricParameters.nx-1) + 0.5;
        int j = unitrandom()*(GeometricParameters.ny-1) + 0.5;
        if (!p[i][j].k_det)
        {
            AssignNewK(i, j, &Field_Generator_Parameters);
            n_filled++;
        }
    }
    cout<<endl;

}

void Grid::Clear()
{
    p.resize(GeometricParameters.nx);
    for (int i = 0; i < GeometricParameters.nx; i++)
    {
        p[i].resize(GeometricParameters.ny);
        for (int j = 0; j < GeometricParameters.ny; j++)
        {
            p[i][j].k_det = false;
            p[i][j].V = CVector(2);
            p[i][j].K = CVector(2);
            p[i][j].K_gauss = CVector(2);
        }
    }
}

void Grid::AssignNewK(int i, int j, field_gen_params *FieldGeneratorParameters)
{
    CNormalDist ND;
    correl_mat_vec M = GetCorrellMatrixVec(i, j, FieldGeneratorParameters);
    double mu;
    double sigma;
    if (M.V_RHS.num == 0)
    {
        mu = 0;
        sigma = 1;
    }
    else
    {
            CMatrix_arma M_inv = inv(M.M_22);
            mu = dotproduct(M_inv*M.V_21, M.V_RHS);
            sigma = 1.0 - dotproduct(M_inv*M.V_21, M.V_21);
    }

    double K_gauss = getnormalrand(mu, sigma);
    p[i][j].k_det = true;
    p[i][j].K_gauss[0] = K_gauss;
    p[i][j].K[0] = FieldGeneratorParameters->inversecdf.interpol(K_gauss);
}

void SetProgressValue(double s)
{
#ifdef QT_version
    main_window->get_ui()->progressBar->setValue(s*100);
    QApplication::processEvents();
#endif // QT_version
    cout << "\r Progress: " << s*100 << "%                              ";
}

correl_mat_vec Grid::GetCorrellMatrixVec(int i, int j,field_gen_params *FieldGeneratorParameters)
{
    correl_mat_vec M;
    vector<ijval> ij = GetClosestCells(GetClosestDeteminedCells(i, j, min(FieldGeneratorParameters->n_filled,FieldGeneratorParameters->max_correl_n)+1,FieldGeneratorParameters), min(FieldGeneratorParameters->n_filled, FieldGeneratorParameters->max_correl_n)+1);
#ifdef Use_Armadillo
    M.M_22 = CMatrix_arma(ij.size() - 1);
    M.V_21 = CVector_arma(ij.size() - 1);
    M.V_RHS = CVector_arma(ij.size() - 1);
#else
    M.M_22 = CMatrix(ij.size() - 1);
    M.V_21 = CVector(ij.size() - 1);
    M.V_RHS = CVector(ij.size() - 1);
#endif //  arma
    for (int ii = 1; ii < int(ij.size()); ii++)
    {
        M.V_21[ii-1] = exp(-sqrt((i - ij[ii].i)*GeometricParameters.dx/ FieldGeneratorParameters->k_correlation_lenght_scale_x*(i - ij[ii].i)*GeometricParameters.dx/ FieldGeneratorParameters->k_correlation_lenght_scale_x + (j - ij[ii].j)*GeometricParameters.dy/ FieldGeneratorParameters->k_correlation_lenght_scale_y*(j - ij[ii].j)*GeometricParameters.dy/ FieldGeneratorParameters->k_correlation_lenght_scale_y));
        M.V_RHS[ii - 1] = p[ij[ii].i][ij[ii].j].K_gauss[0];
        for (int jj = 1; jj < int(ij.size()); jj++)
        {
#ifdef Use_Armadillo
        M.M_22(ii-1,jj-1) = exp(-sqrt((ij[jj].i - ij[ii].i)*GeometricParameters.dx/ FieldGeneratorParameters->k_correlation_lenght_scale_x*(ij[jj].i - ij[ii].i)*GeometricParameters.dx/ FieldGeneratorParameters->k_correlation_lenght_scale_x + (ij[jj].j - ij[ii].j)*GeometricParameters.dy/ FieldGeneratorParameters->k_correlation_lenght_scale_y*(ij[jj].j - ij[ii].j)*GeometricParameters.dy/ FieldGeneratorParameters->k_correlation_lenght_scale_y));
#else
        M.M_22[ii - 1][jj - 1] = exp(-sqrt((ij[jj].i - ij[ii].i)*GeometricParameters.dx/ FieldGeneratorParameters->k_correlation_lenght_scale_x*(ij[jj].i - ij[ii].i)*GeometricParameters.dx/ FieldGeneratorParameters->k_correlation_lenght_scale_x + (ij[jj].j - ij[ii].j)*GeometricParameters.dy/ FieldGeneratorParameters->k_correlation_lenght_scale_y*(ij[jj].j - ij[ii].j)*GeometricParameters.dy/ FieldGeneratorParameters->k_correlation_lenght_scale_y));
#endif // arma
        }
    }
    return M;
}

vector<ijval> Grid::GetClosestDeteminedCells(int i, int j, int n, field_gen_params *FieldGeneratorParameters)
{
    vector<ijval> out;
    double max_dist = 0;
    for (int k = 1; k < max(GeometricParameters.nx, GeometricParameters.ny); k++)
    {
        double min_dist = 1e12;
        int jj = j - k;
        if (jj>=0)
            for (int ii = max(i - k,0); ii <= min(i + k,GeometricParameters.nx-1); ii++)
            {
                double dist2 = ((i - ii)*GeometricParameters.dx/FieldGeneratorParameters->k_correlation_lenght_scale_x*(i - ii)*GeometricParameters.dx/FieldGeneratorParameters->k_correlation_lenght_scale_x + (j - jj)*GeometricParameters.dy/FieldGeneratorParameters->k_correlation_lenght_scale_y*(j - jj)*GeometricParameters.dy/FieldGeneratorParameters->k_correlation_lenght_scale_y);
                min_dist = min(min_dist, dist2);
                if (p[ii][jj].k_det)
                {
                    if ((dist2 < max_dist) || (int(out.size()) < n))
                    {
                        ijval pp;
                        pp.i = ii;
                        pp.j = jj;
                        pp.val = sqrt(dist2);
                        out.push_back(pp);
                        max_dist = max(dist2, max_dist);
                    }
                }

            }
        jj = j + k;
        if (jj < GeometricParameters.ny)
            for (int ii = max(i - k, 0); ii <= min(i + k, GeometricParameters.nx - 1); ii++)
            {
                double dist2 = ((i - ii)*GeometricParameters.dx/FieldGeneratorParameters->k_correlation_lenght_scale_x*(i - ii)*GeometricParameters.dx/FieldGeneratorParameters->k_correlation_lenght_scale_x + (j - jj)*GeometricParameters.dy/FieldGeneratorParameters->k_correlation_lenght_scale_y*(j - jj)*GeometricParameters.dy/FieldGeneratorParameters->k_correlation_lenght_scale_y);
                min_dist = min(min_dist, dist2);
                if (p[ii][jj].k_det)
                {
                    if ((dist2 < max_dist) || (int(out.size()) < n))
                    {
                        ijval pp;
                        pp.i = ii;
                        pp.j = jj;
                        pp.val = sqrt(dist2);
                        out.push_back(pp);
                        max_dist = max(dist2, max_dist);
                    }
                }

            }
        int ii = i - k;
        if (ii >= 0)
            for (int jj = max(j - k+1, 0); jj <= min(j + k-1, GeometricParameters.ny - 1); jj++)
            {
                double dist2 = ((i - ii)*GeometricParameters.dx/FieldGeneratorParameters->k_correlation_lenght_scale_x*(i - ii)*GeometricParameters.dx/FieldGeneratorParameters->k_correlation_lenght_scale_x + (j - jj)*GeometricParameters.dy/FieldGeneratorParameters->k_correlation_lenght_scale_y*(j - jj)*GeometricParameters.dy/FieldGeneratorParameters->k_correlation_lenght_scale_y);
                min_dist = min(min_dist, dist2);
                if (p[ii][jj].k_det)
                {
                    if ((dist2 < max_dist) || (int(out.size()) < n))
                    {
                        ijval pp;
                        pp.i = ii;
                        pp.j = jj;
                        pp.val = sqrt(dist2);
                        out.push_back(pp);
                        max_dist = max(dist2, max_dist);
                    }
                }

            }
        ii = i + k;
        if (ii < GeometricParameters.nx)
            for (int jj = max(j - k + 1, 0); jj <= min(j + k - 1, GeometricParameters.ny - 1); jj++)
            {
                double dist2 = ((i - ii)*GeometricParameters.dx/FieldGeneratorParameters->k_correlation_lenght_scale_x*(i - ii)*GeometricParameters.dx/FieldGeneratorParameters->k_correlation_lenght_scale_x + (j - jj)*GeometricParameters.dy/FieldGeneratorParameters->k_correlation_lenght_scale_y*(j - jj)*GeometricParameters.dy/FieldGeneratorParameters->k_correlation_lenght_scale_y);
                min_dist = min(min_dist, dist2);
                if (p[ii][jj].k_det)
                {
                    if ((dist2 < max_dist) || (int(out.size()) < n))
                    {
                        ijval pp;
                        pp.i = ii;
                        pp.j = jj;
                        pp.val = sqrt(dist2);
                        out.push_back(pp);
                        max_dist = max(dist2, max_dist);
                    }
                }

            }
        if ((int(out.size()) >= n) && min_dist > max_dist)
        {
            k = max(GeometricParameters.nx, GeometricParameters.ny);
            ijval pp;
            pp.i = i;
            pp.j = j;
            pp.val = 0;
            out.push_back(pp);
            return out;
        }
    }
    ijval pp;
    pp.i = i;
    pp.j = j;
    pp.val = 0;
    out.push_back(pp);
    return out;
}

vector<ijval> GetClosestCells(vector<ijval> vec, int n)
{
    vector<ijval> out;
    vector<bool> extracted(vec.size());
    for (int i = 0; i < int(vec.size()); i++) extracted[i] = false;
    int smallest_dist = -1;

    for (int i = 0; i < n; i++)
    {
        double min_dist = 1e12;
        for (int j = 0; j < int(vec.size()); j++)
            if ((vec[j].val < min_dist) && (extracted[j]==false))
            {
                smallest_dist = j;
                min_dist = vec[j].val;
            }
        out.push_back(vec[smallest_dist]);
        extracted[smallest_dist] = true;
    }

    return out;
}
