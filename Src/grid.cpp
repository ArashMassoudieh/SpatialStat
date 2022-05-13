#include "grid.h"
#include "Utilities.h"
#include "NormalDist.h"
#include "environment.h"
#include "Distribution.h"
#include "Copula.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"

Grid::Grid():Interface()
{
    rng_ptr = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rng_ptr, 0);
}

bool Grid::CreateGrid(const map<string,string> &Arguments)
{
    rng_ptr = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rng_ptr, 0);
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
    vector<string> cmds;
    for (map<string,command_parameters>::iterator i=Command::Command_Structures.begin(); i!=Command::Command_Structures.end(); i++)
    {
        if (i->second.Object==object_type::grid)
            cmds.push_back(i->first);
    }
    return cmds;
}

bool Grid::HasCommand(const string &cmd)
{
    if (aquiutils::lookup(Commands(),cmd)!=-1)
        return true;
    else
        return false;
}

FunctionOutPut Grid::Execute(const string &cmd, const map<string,string> &arguments)
{
    FunctionOutPut output;
    if (cmd=="CreateGrid")
    {   output.success = CreateGrid(arguments);
        output.output = this;
    }
    if (cmd=="AssignKField")
    {   output.success = AssignKFieldToGrid(arguments);
        output.output = nullptr;
    }
    if (cmd=="RenormalizeKField")
    {
        output.success = RenormalizeKField(arguments);
        output.output = nullptr;
    }
    if (cmd=="WriteKFieldToVTP")
    {
        output.success = WriteKFieldToVTP(arguments);
        output.output = nullptr;
    }
    if (cmd=="SolveHydro")
    {   output.success = SolveHydro(arguments);
        output.output = nullptr;
    }
    if (cmd=="WriteHydroSolutionToVTP")
    {   output.success = WriteHydroSolutionToVTP(arguments);
        output.output = nullptr;
    }
    if (cmd=="SolveTransport")
    {   output.success = SolveTransport(arguments);
        output.output = nullptr;
    }
    if (cmd=="WriteConcentrationToVTP")
    {   output.success = WriteConcentrationToVTP(arguments);
        output.output = nullptr;
    }
    if (cmd=="GetConcentrationBTCAtX")
    {   output.output = new TimeSeriesD(GetConcentrationBTCAtX(arguments));
        if (dynamic_cast<TimeSeriesD*>(output.output)->n>0)
            output.success = true;
        else
            output.success = false;
    }

    return output;
}

bool Grid::AssignKFieldToGrid(const map<string,string> &Arguments)
{
    if (!parent->Object(Arguments.at("Distribution")))
        return false;



    CDistribution* dist = dynamic_cast<CDistribution*>(parent->Object(Arguments.at("Distribution")));
    Clear();

    field_gen_params Field_Generator_Parameters;

    if (Arguments.count("Maximum_neighboring_nodes")>0)
    {
        Field_Generator_Parameters.max_correl_n = aquiutils::atof(Arguments.at("Maximum_neighboring_nodes"));
    }
    Field_Generator_Parameters.inversecdf = dist->inverse_cumulative;
    Field_Generator_Parameters.k_correlation_lenght_scale_x = aquiutils::atof(Arguments.at("correlation_length_x"));
    Field_Generator_Parameters.k_correlation_lenght_scale_y = aquiutils::atof(Arguments.at("correlation_length_y"));
    srand(time(NULL));
    while (Field_Generator_Parameters.n_filled<GeometricParameters.nx*GeometricParameters.ny)
    {
        if (Field_Generator_Parameters.n_filled%100==0)
        SetProgressValue(double(Field_Generator_Parameters.n_filled) / double(GeometricParameters.nx * GeometricParameters.nx));
        int i = unitrandom()*(GeometricParameters.nx-1) + 0.5;
        int j = unitrandom()*(GeometricParameters.ny-1) + 0.5;
        if (!p[i][j].k_det)
        {
            AssignNewK(i, j, &Field_Generator_Parameters);
            Field_Generator_Parameters.n_filled++;
        }
    }
    cout<<endl;
    return true;

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

bool Grid::RenormalizeKField(const map<string,string> &Arguments)
{
    CDistribution* dist = dynamic_cast<CDistribution*>(parent->Object(Arguments.at("Distribution")));
    RenormalizeK(dist);
    return true;
}

void Grid::RenormalizeK(CDistribution *dist,int k)
{
    CTimeSeries<double> K = GetKValuesToTimeSeries(k);
    double mu = K.mean();
    double std = K.std();
    for (int i = 0; i < GeometricParameters.nx; i++)
        for (int j = 0; j < GeometricParameters.ny; j++)
            p[i][j].K_gauss[k] = (p[i][j].K_gauss[k] - mu) / std;
    RemapKFieldBasedonMarginalDistribution(dist,k);
}

void Grid::RemapKFieldBasedonMarginalDistribution(CDistribution *dist,int k)
{
    for (int i = 0; i < GeometricParameters.nx; i++)
        for (int j = 0; j < GeometricParameters.ny; j++)
            p[i][j].K[k] = MapToMarginalDistribution(getnormalcdf(p[i][j].K_gauss[k]),dist);

}

CTimeSeries<double> Grid::GetKValuesToTimeSeries(int k)
{
    CTimeSeries<double> out;
    for (int i = 0; i < GeometricParameters.nx; i++)
        for (int j = 0; j < GeometricParameters.ny; j++)
            out.append(i + j * 10000, p[i][j].K_gauss[k]);

    return out;
}

double Grid::MapToMarginalDistribution(const double &u, CDistribution *dist)
{
    return dist->InverseCumulativeValue(u);
}

bool Grid::WriteKFieldToVTP(const string &filename, const double &z_factor, bool _log)
{
    vtkSmartPointer<vtkPoints> points_3 =
        vtkSmartPointer<vtkPoints>::New();

    double maxk = max_K();
    double xx, yy, zz;
    vtkSmartPointer<vtkFloatArray> values =
        vtkSmartPointer<vtkFloatArray>::New();

    values->SetNumberOfComponents(1);
    if (_log)
        values->SetName("Hydraulic Conductivity");
    else
        values->SetName("Log Hydraulic Conductivity");

    for (unsigned int x = 0; x < GeometricParameters.nx; x++)
    {
        for (unsigned int y = 0; y < GeometricParameters.ny; y++)
        {
            xx = x*GeometricParameters.dx;
            yy = y*GeometricParameters.dy;
            zz = p[x][y].K[0];
            if (!_log)
            {
                float t[1] = { float(zz) };
                points_3->InsertNextPoint(xx, yy, zz / maxk*z_factor);
                values->InsertNextTupleValue(t);
            }
            else
            {
                float t[1] = { float(zz)};
                points_3->InsertNextPoint(xx, yy, log(*t/maxk)*z_factor);
                values->InsertNextTupleValue(t);
            }
        }
    }

    // Add the grid points to a polydata object
    vtkSmartPointer<vtkPolyData> inputPolyData =
        vtkSmartPointer<vtkPolyData>::New();
    inputPolyData->SetPoints(points_3);

    // Triangulate the grid points
    vtkSmartPointer<vtkDelaunay2D> delaunay =
        vtkSmartPointer<vtkDelaunay2D>::New();
#if VTK_MAJOR_VERSION <= 5
    delaunay->SetInput(inputPolyData);
#else
    delaunay->SetInputData(inputPolyData);
#endif
    delaunay->Update();
    vtkPolyData* outputPolyData = delaunay->GetOutput();

    double bounds[6];
    outputPolyData->GetBounds(bounds);

    // Find min and max z
    double minz = bounds[4];
    double maxz = bounds[5];

    std::cout << "minz: " << minz << std::endl;
    std::cout << "maxz: " << maxz << std::endl;

    // Create the color map
    vtkSmartPointer<vtkLookupTable> colorLookupTable =
        vtkSmartPointer<vtkLookupTable>::New();
    colorLookupTable->SetTableRange(minz, maxz);
    colorLookupTable->Build();

    // Generate the colors for each point based on the color map
    vtkSmartPointer<vtkUnsignedCharArray> colors_2 =
        vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors_2->SetNumberOfComponents(3);
    colors_2->SetName("Colors");

//	std::cout << "There are " << outputPolyData->GetNumberOfPoints()
//		<< " points." << std::endl;

    for (int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
    {
        double p[3];
        outputPolyData->GetPoint(i, p);

        double dcolor[3];
        colorLookupTable->GetColor(p[2], dcolor);
        //std::cout << "dcolor: "
        //	<< dcolor[0] << " "
        //	<< dcolor[1] << " "
        //	<< dcolor[2] << std::endl;
        unsigned char color[3];
        for (unsigned int j = 0; j < 3; j++)
        {
            color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
        }
        //std::cout << "color: "
        //	<< (int)color[0] << " "
        //	<< (int)color[1] << " "
        //	<< (int)color[2] << std::endl;

        colors_2->InsertNextTupleValue(color);
    }

    outputPolyData->GetPointData()->SetScalars(values);


    //Append the two meshes
    vtkSmartPointer<vtkAppendPolyData> appendFilter =
        vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
    appendFilter->AddInputConnection(input1->GetProducerPort());
    appendFilter->AddInputConnection(input2->GetProducerPort());
#else
    //appendFilter->AddInputData(polydata);
    //appendFilter->AddInputData(polydata_1);
    appendFilter->AddInputData(outputPolyData);
#endif
    appendFilter->Update();


    // Visualization
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
    mapper->SetInputConnection(polydata->GetProducerPort());
#else
    mapper->SetInputConnection(appendFilter->GetOutputPort());
    //mapper->SetInputData(polydata_1);
#endif

    vtkSmartPointer<vtkActor> actor =
        vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(5);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(mapper->GetInput());
    // This is set so we can see the data in a text editor.
    writer->SetDataModeToAscii();
    writer->Write();
    return true;
}

double Grid::max_K()
{
    double max_k = -1e23;
    for (int i = 0; i < GeometricParameters.nx; i++)
        for (int j = 0; j < GeometricParameters.ny; j++)
            max_k = max(max_k, p[i][j].K[0]);
    return max_k;

}
double Grid::min_K()
{
    double min_k = 1e23;
    for (int i = 0; i < GeometricParameters.nx; i++)
        for (int j = 0; j < GeometricParameters.ny; j++)
            min_k = min(min_k, p[i][j].K[0]);
    return min_k;
}

double Grid::max_vx()
{
    double max_k = -1e23;
    for (int i = 0; i < GeometricParameters.nx; i++)
        for (int j = 0; j < GeometricParameters.ny; j++)
            max_k = max(max_k, p[i][j].V[0]);
    return max_k;
}

double Grid::min_vx()
{
    double min_k = 1e23;
    for (int i = 0; i < GeometricParameters.nx; i++)
        for (int j = 0; j < GeometricParameters.ny; j++)
            min_k = min(min_k, p[i][j].V[0]);
    return min_k;
}
double Grid::max_vy()
{
    double max_k = -1e23;
    for (int i = 0; i < GeometricParameters.nx; i++)
        for (int j = 0; j < GeometricParameters.ny; j++)
            max_k = max(max_k, p[i][j].V[1]);
    return max_k;

}
double Grid::min_vy()
{
    double min_k = 1e23;
    for (int i = 0; i < GeometricParameters.nx; i++)
        for (int j = 0; j < GeometricParameters.ny; j++)
            min_k = min(min_k, p[i][j].V[1]);
    return min_k;
}

bool Grid::WriteKFieldToVTP(const map<string,string> &Arguments)
{
    if (Arguments.count("filename")==0)
        return false;
    double z_factor=1;
    bool log_scale = false;
    if (Arguments.count("z_factor")>0)
        z_factor = aquiutils::atof(Arguments.at("z_factor"));
    if (Arguments.count("log_scale")>0)
        log_scale = aquiutils::atoi(Arguments.at("log_scale"));
    return WriteKFieldToVTP(Arguments.at("filename"),z_factor,log_scale);
}

bool Grid::SolveHydro(const map<string,string> &Arguments)
{
    double l_boundary = 1;
    double r_boundary = 0;
    if (Arguments.count("l_boundary")>0)
        l_boundary = aquiutils::atof(Arguments.at("l_boundary"));
    if (Arguments.count("r_boundary")>0)
        r_boundary = aquiutils::atof(Arguments.at("r_boundary"));

    return SolveHydro(l_boundary,r_boundary);
}

bool Grid::WriteHydroSolutionToVTP(const map<string,string> &Arguments)
{
    if (Arguments.count("filename")==0)
        return false;
    double z_factor=1;
    bool log_scale = false;
    if (Arguments.count("z_factor")>0)
        z_factor = aquiutils::atof(Arguments.at("z_factor"));
    if (Arguments.count("log_scale")>0)
        log_scale = aquiutils::atoi(Arguments.at("log_scale"));
    return WriteHydroSolutionToVTP(Arguments.at("filename"),z_factor,log_scale);
}

bool Grid::SolveHydro(const double &leftboundary, const double &rightboundary)
{
    SetProgressValue(0);
    CMatrix_arma_sp K = CreateStiffnessMatrixHydro();
    CVector_arma V = CreateRHSHydro_ARMA(leftboundary,rightboundary);
    CVector_arma S = solve_ar(K, V);

    H = CMatrix(GeometricParameters.nx+1,GeometricParameters.ny+1);
    vx = CMatrix(GeometricParameters.nx, GeometricParameters.ny-1);
    vy = CMatrix(GeometricParameters.nx - 1, GeometricParameters.ny);
    for (int i = 0; i < GeometricParameters.nx+1; i++)
        for (int j = 0; j < GeometricParameters.ny+1; j++)
            H[i][j] = S[get_cell_no(i, j)];

    for (int i = 0; i < GeometricParameters.nx; i++)
    for (int j = 0; j < GeometricParameters.ny; j++)
    {
            double Kx1 = 0.5*(p[i][max(j - 1,0)].K[0] + p[i][j].K[0]);
            double Kx2 = 0.5*(p[i][j].K[0] + p[i][min(j+1,GeometricParameters.ny-1)].K[0]);
            double Ky1 = 0.5*(p[max(i - 1, 0)][j].K[0] + p[i][j].K[0]);
            double Ky2 = 0.5*(p[i][j].K[0] + p[min(i+1,GeometricParameters.nx-1)][j].K[0]);
            p[i][j].V[0] = -(0.5*Kx1 * (H[i + 1][j] - H[i][j]) / GeometricParameters.dx + 0.5*Kx2 * (H[i + 1][j+1] - H[i][j+1]) / GeometricParameters.dx);
            p[i][j].V[1] = -(0.5*Ky1 * (H[i][j+1] - H[i][j]) / GeometricParameters.dy + 0.5*Ky2 * (H[i + 1][j + 1] - H[i+1][j]) / GeometricParameters.dy);
            p[i][j].Vbx = -Kx1*(H[i+1][j] - H[i][j]) / GeometricParameters.dx;
            p[i][j].Vtx = -Kx2*(H[i+1][j+1] - H[i][j+1]) / GeometricParameters.dx;
            p[i][j].Vby = -Ky1*(H[i][j+1] - H[i][j]) / GeometricParameters.dy;
            p[i][j].Vtx = -Ky2*(H[i+1][j+1] - H[i+1][j]) / GeometricParameters.dy;
    }

    for (int i=0; i<GeometricParameters.nx; i++)
        for (int j = 0; j < GeometricParameters.ny - 1; j++)
        {
            double Kx = 0.5*(p[i][j].K[0] + p[i][j + 1].K[0]);
            vx[i][j] = -Kx * (H[i + 1][j + 1] - H[i][j + 1]) / GeometricParameters.dx;
        }

    SetProgressValue(0.5);
    for (int i = 0; i<GeometricParameters.nx-1; i++)
        for (int j = 0; j < GeometricParameters.ny; j++)
        {
            double Ky = 0.5*(p[i][j].K[0] + p[i+1][j].K[0]);
            vy[i][j] = -Ky * (H[i + 1][j+1] - H[i+1][j]) / GeometricParameters.dy;
        }

    SetProgressValue(1);

    max_v_x = max_vx();
    min_v_x = min_vx();
    return true;
    cout<<endl;
}



CMatrix_arma_sp Grid::CreateStiffnessMatrixHydro()
{
    string averaging = "arithmetic";
//	qDebug() << "Creating stiffness matrix" << endl;
        CMatrix_arma_sp K((GeometricParameters.nx+1)*(GeometricParameters.ny+1), (GeometricParameters.nx + 1)*(GeometricParameters.ny + 1));
//	qDebug() << "Stiffness matrix created" << endl;
    for (int i = 1; i < GeometricParameters.nx; i++)
    {
        for (int j = 1; j < GeometricParameters.ny; j++)
        {
            K.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) = aquiutils::avg(p[i-1][j-1].K[0], p[i-1][j].K[0], averaging) / (GeometricParameters.dx*GeometricParameters.dx);
            K.matr(get_cell_no(i, j), get_cell_no(i, j)) = -(aquiutils::avg(p[i-1][j-1].K[0], p[i-1][j].K[0], averaging) + aquiutils::avg(p[i][j-1].K[0], p[i][j].K[0], averaging)) / (GeometricParameters.dx*GeometricParameters.dx);
            K.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) = aquiutils::avg(p[i][j-1].K[0], p[i][j].K[0], averaging) / (GeometricParameters.dx*GeometricParameters.dx);

            K.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) = aquiutils::avg(p[i-1][j-1].K[0], p[i][j-1].K[0], averaging) / (GeometricParameters.dy*GeometricParameters.dy);
            K.matr(get_cell_no(i, j), get_cell_no(i, j)) += -(aquiutils::avg(p[i-1][j-1].K[0], p[i][j-1].K[0], averaging) + aquiutils::avg(p[i-1][j].K[0], p[i][j].K[0], averaging)) / (GeometricParameters.dy*GeometricParameters.dy);
            K.matr(get_cell_no(i, j), get_cell_no(i, j + 1)) = aquiutils::avg(p[i-1][j].K[0], p[i][j].K[0], averaging) / (GeometricParameters.dy*GeometricParameters.dy);
        }
        // top boundary
        int j = GeometricParameters.ny;
        K.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) = 1;
        K.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;

        // bottom boundary
        j = 0;
        K.matr(get_cell_no(i, j), get_cell_no(i, j + 1)) = 1;
        K.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;
    }

    //left boundary
    int i = 0;
    for (int j = 0; j < GeometricParameters.ny + 1; j++)
    {
        K.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
        K.matr(get_cell_no(i, j), get_cell_no(i+1, j)) = 1;

    }

    //right boundary
    i = GeometricParameters.nx;
    for (int j = 0; j < GeometricParameters.ny + 1; j++)
    {
        K.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
        K.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) = 1;
        }

    return K;

}

int Grid::get_cell_no(int i, int j)
{
    return i*(GeometricParameters.ny+1) + j;
}

CVector_arma Grid::CreateRHSHydro_ARMA(const double &leftboundary, const double &rightboundary)
{
    CVector_arma V((GeometricParameters.nx + 1)*(GeometricParameters.ny + 1));

    //left boundary
    int i = 0;
    for (int j = 0; j < GeometricParameters.ny + 1; j++)
        V[get_cell_no(i, j)] = 2*leftboundary;


    //right boundary
    i = GeometricParameters.nx;
    for (int j = 0; j < GeometricParameters.ny + 1; j++)
        V[get_cell_no(i, j)] = 2*rightboundary;


    return V;
}

CVector Grid::CreateRHSHydro(const double &leftboundary, const double &rightboundary)
{
    CVector V((GeometricParameters.nx + 1)*(GeometricParameters.ny + 1));

    //left boundary
    int i = 0;
    for (int j = 0; j < GeometricParameters.ny + 1; j++)
        V[get_cell_no(i, j)] = 2 * leftboundary;


    //right boundary
    i = GeometricParameters.nx;
    for (int j = 0; j < GeometricParameters.ny + 1; j++)
        V[get_cell_no(i, j)] = 2 * rightboundary;


    return V;

}

bool Grid::WriteHydroSolutionToVTP(const string &filename, const double &z_factor, bool _log)
{
    vtkSmartPointer<vtkPoints> points_3 =
        vtkSmartPointer<vtkPoints>::New();

    double maxk = max_K();
    double xx, yy, zz;
    vtkSmartPointer<vtkFloatArray> values =
        vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkFloatArray> _H =
        vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkFloatArray> v =
        vtkSmartPointer<vtkFloatArray>::New();

    values->SetNumberOfComponents(1);
    _H->SetNumberOfComponents(1);
    v->SetNumberOfComponents(2);

    if (_log)
        values->SetName("Hydraulic Conductivity");
    else
        values->SetName("Log Hydraulic Conductivity");

    _H->SetName("Hydraulic Head");
    v->SetName("Velocity");

    for (unsigned int x = 0; x < GeometricParameters.nx; x++)
    {
        for (unsigned int y = 0; y < GeometricParameters.ny; y++)
        {
            xx = x*GeometricParameters.dx;
            yy = y*GeometricParameters.dy;
            zz = p[x][y].K[0];
            if (!_log)
            {
                float t[1] = { float(zz) };
                float vv[2] = { float(p[x][y].V[0]), float(p[x][y].V[1]) };
                float HH[1] = { float(H[x][y]) };
                points_3->InsertNextPoint(xx, yy, zz / maxk*z_factor);
                values->InsertNextTupleValue(t);
                v->InsertNextTupleValue(vv);
                _H->InsertNextTupleValue(HH);

            }
            else
            {
                float t[1] = { float(zz) };
                float vv[2] = { float(p[x][y].V[0]), float(p[x][y].V[1]) };
                float HH[1] = { float(H[x][y]) };
                points_3->InsertNextPoint(xx, yy, log(*t / maxk)*z_factor);
                values->InsertNextTupleValue(t);
                v->InsertNextTupleValue(vv);
                _H->InsertNextTupleValue(HH);

            }
        }
    }

    // Add the grid points to a polydata object
    vtkSmartPointer<vtkPolyData> inputPolyData =
        vtkSmartPointer<vtkPolyData>::New();
    inputPolyData->SetPoints(points_3);

    // Triangulate the grid points
    vtkSmartPointer<vtkDelaunay2D> delaunay =
        vtkSmartPointer<vtkDelaunay2D>::New();
#if VTK_MAJOR_VERSION <= 5
    delaunay->SetInput(inputPolyData);
#else
    delaunay->SetInputData(inputPolyData);
#endif
    delaunay->Update();
    vtkPolyData* outputPolyData = delaunay->GetOutput();

    double bounds[6];
    outputPolyData->GetBounds(bounds);

    // Find min and max z
    double minz = bounds[4];
    double maxz = bounds[5];

    //std::cout << "minz: " << minz << std::endl;
    //std::cout << "maxz: " << maxz << std::endl;

    // Create the color map
    vtkSmartPointer<vtkLookupTable> colorLookupTable =
        vtkSmartPointer<vtkLookupTable>::New();
    colorLookupTable->SetTableRange(minz, maxz);
    colorLookupTable->Build();

    // Generate the colors for each point based on the color map
    vtkSmartPointer<vtkUnsignedCharArray> colors_2 =
        vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors_2->SetNumberOfComponents(3);
    colors_2->SetName("Colors");

//	std::cout << "There are " << outputPolyData->GetNumberOfPoints()
//		<< " points." << std::endl;

    for (int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
    {
        double p[3];
        outputPolyData->GetPoint(i, p);

        double dcolor[3];
        colorLookupTable->GetColor(p[2], dcolor);
        //std::cout << "dcolor: "
        //	<< dcolor[0] << " "
        //	<< dcolor[1] << " "
        //	<< dcolor[2] << std::endl;
        unsigned char color[3];
        for (unsigned int j = 0; j < 3; j++)
        {
            color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
        }
//		std::cout << "color: "
//			<< (int)color[0] << " "
//			<< (int)color[1] << " "
//			<< (int)color[2] << std::endl;

        colors_2->InsertNextTupleValue(color);
    }

    outputPolyData->GetPointData()->SetScalars(values);
    outputPolyData->GetPointData()->AddArray(_H);
    outputPolyData->GetPointData()->AddArray(v);

    //Append the two meshes
    vtkSmartPointer<vtkAppendPolyData> appendFilter =
        vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
    appendFilter->AddInputConnection(input1->GetProducerPort());
    appendFilter->AddInputConnection(input2->GetProducerPort());
#else
    //appendFilter->AddInputData(polydata);
    //appendFilter->AddInputData(polydata_1);
    appendFilter->AddInputData(outputPolyData);
#endif
    appendFilter->Update();


    // Visualization
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
    mapper->SetInputConnection(polydata->GetProducerPort());
#else
    mapper->SetInputConnection(appendFilter->GetOutputPort());
    //mapper->SetInputData(polydata_1);
#endif

    vtkSmartPointer<vtkActor> actor =
        vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(5);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(mapper->GetInput());
    // This is set so we can see the data in a text editor.
    writer->SetDataModeToAscii();
    writer->Write();
    return true;

}

bool Grid::SolveTransport(const double &t_end, const vector<double> &decay_coeff, const vector<double> &decay_order)
{
    CreateTransportKMatrix(TransportParameters.dt, TransportParameters.D, TransportParameters.time_weight);
    C.resize(TransportParameters.numberofspecies);
    for (int species_counter=0; species_counter<TransportParameters.numberofspecies; species_counter++)
    {
        C[species_counter] = CMatrix(GeometricParameters.nx + 1, GeometricParameters.ny + 1);// = leftboundary_C;
    }
    for (int i = 0; i < GeometricParameters.nx; i++)
        for (int j = 0; j < GeometricParameters.ny; j++)
            p[i][j].C.resize(int(t_end/TransportParameters.dt), TransportParameters.numberofspecies);
    CMatrix_arma_sp K = TransportParameters.KD + TransportParameters.Kt + TransportParameters.Kv;
    //Kt.writetofile(pathout + "Kt_matrix.txt");
    //KD.writetofile(pathout + "KD_matrix.txt");
    //Kv.writetofile(pathout + "Kv_matrix.txt");
    //K.writetofile(pathout + "transport_matrix.txt");
    SetProgressValue(0);
        int counter=0;
        for (double t = 0; t < t_end; t += TransportParameters.dt)
        {
            for (int species_counter=0; species_counter<TransportParameters.numberofspecies; species_counter++)
            {   CVector_arma RHS = CreateTransportRHS(species_counter, TransportParameters.dt, TransportParameters.time_weight, TransportParameters.D, decay_coeff[species_counter], decay_order[species_counter]);
                CVector_arma S = solve_ar(K, RHS);

                for (int i=0; i<GeometricParameters.nx+1; i++)
                    for (int j=0; j<GeometricParameters.ny+1; j++)
                        C[species_counter][i][j] = S[get_cell_no(i, j)];
            }
            for (int i = 0; i < GeometricParameters.nx; i++)
                for (int j = 0; j < GeometricParameters.ny; j++)
                {
                    vector<double> cc;
                    for (int species_counter = 0; species_counter<TransportParameters.numberofspecies; species_counter++)
                        p[i][j].C.setvalue(counter,species_counter,C[species_counter][i+1][j]);
                }
            counter++;
            SetProgressValue(t / t_end);

        }
        return true;
}

void Grid::CreateTransportKMatrix(const double &dt, const double &D, const double &weight)
{
    string averaging = "arithmetic";
    TransportParameters.Kv = CMatrix_arma_sp((GeometricParameters.nx + 1)*(GeometricParameters.ny + 1), (GeometricParameters.nx + 1)*(GeometricParameters.ny + 1));
    TransportParameters.KD = CMatrix_arma_sp((GeometricParameters.nx + 1)*(GeometricParameters.ny + 1), (GeometricParameters.nx + 1)*(GeometricParameters.ny + 1));
    TransportParameters.Kt = CMatrix_arma_sp((GeometricParameters.nx + 1)*(GeometricParameters.ny + 1), (GeometricParameters.nx + 1)*(GeometricParameters.ny + 1));

    for (int i = 1; i < GeometricParameters.nx; i++)
    {
        for (int j = 1; j < GeometricParameters.ny; j++)
        {
            TransportParameters.Kv.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) += -weight*aquiutils::Pos(vx[i-1][j-1]) / (GeometricParameters.dx);
            TransportParameters.Kv.matr(get_cell_no(i, j), get_cell_no(i, j)) += weight*(aquiutils::Neg(vx[i-1][j-1])/GeometricParameters.dx + aquiutils::Pos(vx[i][j-1])/GeometricParameters.dx);
            TransportParameters.Kv.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) += -weight*aquiutils::Neg(vx[i][j - 1]) / (GeometricParameters.dx);

            TransportParameters.Kv.matr(get_cell_no(i, j), get_cell_no(i, j-1)) += -weight*aquiutils::Pos(vy[i - 1][j - 1]) / GeometricParameters.dy;
            TransportParameters.Kv.matr(get_cell_no(i, j), get_cell_no(i, j)) += weight*(aquiutils::Neg(vy[i - 1][j - 1]) / GeometricParameters.dy + aquiutils::Pos(vy[i-1][j]) / GeometricParameters.dy);
            TransportParameters.Kv.matr(get_cell_no(i, j), get_cell_no(i, j+1)) += -weight*aquiutils::Neg(vy[i-1][j]) / GeometricParameters.dy;

            TransportParameters.KD.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) += -weight*D / (GeometricParameters.dx*GeometricParameters.dx);
            TransportParameters.KD.matr(get_cell_no(i, j), get_cell_no(i, j)) += 2*weight*D / (GeometricParameters.dx*GeometricParameters.dx);
            TransportParameters.KD.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) += -weight*D / (GeometricParameters.dx*GeometricParameters.dx);

            TransportParameters.KD.matr(get_cell_no(i, j), get_cell_no(i, j-1)) += -weight*D / (GeometricParameters.dy*GeometricParameters.dy);
            TransportParameters.KD.matr(get_cell_no(i, j), get_cell_no(i, j)) += 2*weight*D / (GeometricParameters.dy*GeometricParameters.dy);
            TransportParameters.KD.matr(get_cell_no(i, j), get_cell_no(i, j+1)) += -weight*D / (GeometricParameters.dy*GeometricParameters.dy);

            TransportParameters.Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) += 1.0 / dt;

        }
        // top boundary
        int j = GeometricParameters.ny;
        TransportParameters.Kt.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) = 1;
        TransportParameters.Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;

        // bottom boundary
        j = 0;
        TransportParameters.Kt.matr(get_cell_no(i, j), get_cell_no(i, j + 1)) = 1;
        TransportParameters.Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;
    }

    //left boundary
    int i = 0;
    for (int j = 0; j < GeometricParameters.ny + 1; j++)
    {
        TransportParameters.Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
        TransportParameters.Kt.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) = 1;
    }

    //right boundary
    i = GeometricParameters.nx;
    for (int j = 0; j < GeometricParameters.ny + 1; j++)
    {
        TransportParameters.Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
        TransportParameters.Kt.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) = -1;
    }


}

CVector_arma Grid::CreateTransportRHS(int species_counter, const double &dt, const double &weight, const double &D, const double &decay_coefficient, const double &decay_order)
{
    CVector_arma RHS((GeometricParameters.nx + 1)*(GeometricParameters.ny + 1));
    for (int i = 1; i < GeometricParameters.nx; i++)
    {
        for (int j = 1; j < GeometricParameters.ny; j++)
        {
            double rhs = 0;
            rhs += (1-weight)*aquiutils::Pos(vx[i - 1][j - 1]) / (GeometricParameters.dx)*C[species_counter][i-1][j];
            rhs += -(1-weight)*(aquiutils::Neg(vx[i - 1][j - 1]) / GeometricParameters.dx + aquiutils::Pos(vx[i][j-1]) / GeometricParameters.dx)*C[species_counter][i][j];
            rhs += (1-weight)*aquiutils::Neg(vx[i][j - 1]) / (GeometricParameters.dx)*C[species_counter][i+1][j];

            rhs += (1-weight)*aquiutils::Pos(vy[i - 1][j - 1]) / GeometricParameters.dy*C[species_counter][i][j-1];
            rhs += -(1 - weight)*(aquiutils::Neg(vy[i - 1][j - 1]) / GeometricParameters.dy + aquiutils::Pos(vy[i - 1][j]) / GeometricParameters.dy)*C[species_counter][i][j];
            rhs += (1-weight)*aquiutils::Neg(vy[i - 1][j]) / GeometricParameters.dy*C[species_counter][i][j+1];

            rhs += (1-weight)*D / (GeometricParameters.dx*GeometricParameters.dx)*C[species_counter][i-1][j];
            rhs += -2 * (1-weight)*D / (GeometricParameters.dx*GeometricParameters.dx)*C[species_counter][i][j];
            rhs += (1-weight)*D / (GeometricParameters.dx*GeometricParameters.dx)*C[species_counter][i + 1][j];

            rhs += (1-weight)*D / (GeometricParameters.dy*GeometricParameters.dy)*C[species_counter][i][j-1];
            rhs += -2 * (1-weight)*D / (GeometricParameters.dy*GeometricParameters.dy)*C[species_counter][i][j];
            rhs += (1-weight)*D / (GeometricParameters.dy*GeometricParameters.dy)*C[species_counter][i][j+1];

            if (TransportParameters.numberofspecies==1)
                rhs += 1.0 / dt*C[species_counter][i][j] - decay_coefficient*pow(C[species_counter][i][j], decay_order);
            else if (TransportParameters.numberofspecies==2)
                rhs += 1.0 / dt*C[species_counter][i][j] - decay_coefficient*C[0][i][j]*C[1][i][j];
            else if (TransportParameters.numberofspecies==3)
            {
                if (species_counter<2)
                {
                    rhs += 1.0 / dt*C[species_counter][i][j] - decay_coefficient*C[0][i][j]*C[1][i][j];
                }
                else
                    rhs += 1.0 / dt*C[species_counter][i][j] + decay_coefficient*C[0][i][j]*C[1][i][j];
            }
            RHS[get_cell_no(i, j)] = rhs;
        }
        // top boundary
        int j = GeometricParameters.ny;
        RHS[get_cell_no(i, j)] = 0;


        // bottom boundary
        j = 0;
        RHS[get_cell_no(i, j)] = 0;

    }

    //left boundary
    int i = 0;
    for (int j = 0; j < GeometricParameters.ny + 1; j++)
        RHS[get_cell_no(i, j)] = 2*TransportParameters.leftboundary_C[species_counter];


    //right boundary
    i = GeometricParameters.nx;
    for (int j = 0; j < GeometricParameters.ny + 1; j++)
        RHS[get_cell_no(i, j)] = 0;

    return RHS;

}

bool Grid::SolveTransport(const map<string,string> &Arguments)
{
    TransportParameters.numberofspecies=1;
    if (Arguments.count("nspecies")>0)
        TransportParameters.numberofspecies = aquiutils::atoi(Arguments.at("nspecies"));
    vector<double> decay_coeff(TransportParameters.numberofspecies);
    vector<double> decay_order(TransportParameters.numberofspecies);

    if (Arguments.count("time_weight")>0)
        TransportParameters.time_weight = aquiutils::atof(Arguments.at("time_weight"));
    if (Arguments.count("l_boundary")>0)
        TransportParameters.leftboundary_C = aquiutils::ATOF(aquiutils::split(Arguments.at("l_boundary"),';'));
    if (Arguments.count("diffusion")>0)
        TransportParameters.D = aquiutils::atof(Arguments.at("diffusion"));
    if (Arguments.count("dt")>0)
        TransportParameters.dt = aquiutils::atof(Arguments.at("dt"));
    if (Arguments.count("decay_coeff")>0)
        decay_coeff = aquiutils::ATOF(aquiutils::split(Arguments.at("decay_coeff"),';'));
    if (Arguments.count("decay_order")>0)
        decay_order = aquiutils::ATOF(aquiutils::split(Arguments.at("decay_order"),';'));
    double t_end = 1;
    if (Arguments.count("t_end")>0)
        t_end = aquiutils::atof(Arguments.at("t_end"));

    return SolveTransport(t_end,decay_coeff,decay_order);
}

bool Grid::WriteConcentrationToVTP(int species_counter, const string &filename, const double &z_factor, bool _log, const vector<double> &t)
{
    bool result=true;
    for (int i = 0; i < t.size(); i++)
    {
        string filename_1 = aquiutils::split(filename, '.')[0] + "_" + aquiutils::numbertostring(i) +"." + aquiutils::split(filename, '.')[1];
        result &= WriteConcentrationToVTP(species_counter, filename_1, z_factor, _log, t[i]);
        SetProgressValue(double(i)/double(t.size()));
    }
    return result;
}

bool Grid::WriteConcentrationToVTP(int species_counter, const string &filename, const double &z_factor, bool _log, const double &t)
{
    vtkSmartPointer<vtkPoints> points_3 =
        vtkSmartPointer<vtkPoints>::New();


    double xx, yy, zz;
    vtkSmartPointer<vtkFloatArray> values =
        vtkSmartPointer<vtkFloatArray>::New();


    values->SetNumberOfComponents(1);

    if (_log)
        values->SetName("Log Concentration");
    else
        values->SetName("Concentration");

    for (unsigned int x = 0; x < GeometricParameters.nx; x++)
    {
        for (unsigned int y = 0; y < GeometricParameters.ny; y++)
        {
            xx = x*GeometricParameters.dx;
            yy = y*GeometricParameters.dy;
            zz = p[x][y].C[int(t/TransportParameters.dt)][species_counter];
            if (!_log)
            {
                float CC[1] = { float(zz) };
                points_3->InsertNextPoint(xx, yy, zz * z_factor);
                values->InsertNextTupleValue(CC);
            }
            else
            {
                float CC[1] = { float(zz) };
                points_3->InsertNextPoint(xx, yy, log(max(zz,1e-25))*z_factor);
                values->InsertNextTupleValue(CC);
            }
        }
    }

    // Add the grid points to a polydata object
    vtkSmartPointer<vtkPolyData> inputPolyData =
        vtkSmartPointer<vtkPolyData>::New();
    inputPolyData->SetPoints(points_3);

    // Triangulate the grid points
    vtkSmartPointer<vtkDelaunay2D> delaunay =
        vtkSmartPointer<vtkDelaunay2D>::New();
#if VTK_MAJOR_VERSION <= 5
    delaunay->SetInput(inputPolyData);
#else
    delaunay->SetInputData(inputPolyData);
#endif
    delaunay->Update();
    vtkPolyData* outputPolyData = delaunay->GetOutput();

    double bounds[6];
    outputPolyData->GetBounds(bounds);

    // Find min and max z
    double minz = bounds[4];
    double maxz = bounds[5];

    //std::cout << "minz: " << minz << std::endl;
    //std::cout << "maxz: " << maxz << std::endl;

    // Create the color map
    vtkSmartPointer<vtkLookupTable> colorLookupTable =
        vtkSmartPointer<vtkLookupTable>::New();
    colorLookupTable->SetTableRange(minz, maxz);
    colorLookupTable->Build();

    // Generate the colors for each point based on the color map
    vtkSmartPointer<vtkUnsignedCharArray> colors_2 =
        vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors_2->SetNumberOfComponents(3);
    colors_2->SetName("Colors");

//	std::cout << "There are " << outputPolyData->GetNumberOfPoints()
//		<< " points." << std::endl;

    for (int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
    {
        double p[3];
        outputPolyData->GetPoint(i, p);

        double dcolor[3];
        colorLookupTable->GetColor(p[2], dcolor);
        //std::cout << "dcolor: "
        //	<< dcolor[0] << " "
        //	<< dcolor[1] << " "
        //	<< dcolor[2] << std::endl;
        unsigned char color[3];
        for (unsigned int j = 0; j < 3; j++)
        {
            color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
        }
        //std::cout << "color: "
        //	<< (int)color[0] << " "
        //	<< (int)color[1] << " "
        //	<< (int)color[2] << std::endl;

        colors_2->InsertNextTupleValue(color);
    }

    outputPolyData->GetPointData()->SetScalars(values);

    //Append the two meshes
    vtkSmartPointer<vtkAppendPolyData> appendFilter =
        vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
    appendFilter->AddInputConnection(input1->GetProducerPort());
    appendFilter->AddInputConnection(input2->GetProducerPort());
#else
    //appendFilter->AddInputData(polydata);
    //appendFilter->AddInputData(polydata_1);
    appendFilter->AddInputData(outputPolyData);
#endif
    appendFilter->Update();


    // Visualization
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
    mapper->SetInputConnection(polydata->GetProducerPort());
#else
    mapper->SetInputConnection(appendFilter->GetOutputPort());
    //mapper->SetInputData(polydata_1);
#endif

    vtkSmartPointer<vtkActor> actor =
        vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(5);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(mapper->GetInput());
    // This is set so we can see the data in a text editor.
    writer->SetDataModeToAscii();
    writer->Write();
    return true;
}

bool Grid::WriteConcentrationToVTP(const map<string,string> &Arguments)
{
    if (Arguments.count("filename") == 0) return false;
    int species_counter;
    if (TransportParameters.numberofspecies==1)
        species_counter = 0;
    else if (Arguments.count("species")==0)
        species_counter = 0;
    else
        species_counter = aquiutils::atoi(Arguments.at("species"));
    double concentration_interval = -1;

    if (Arguments.count("interval") > 0)
        concentration_interval = aquiutils::atof(Arguments.at("interval"));
    else
        concentration_interval = 1;
    vector<double> intervals;
    for (int i = 0; i < p[0][0].C.size(); i+=concentration_interval)
    {
            intervals.push_back(i*TransportParameters.dt);
    }

    double z_factor=1;
    if (Arguments.count("z_factor") > 0)
        z_factor = aquiutils::atof(Arguments.at("z_factor"));

    bool _log = false;
    if (Arguments.count("log") > 0)
        _log = aquiutils::atof(Arguments.at("log"));

    return WriteConcentrationToVTP(species_counter,Arguments.at("filename"), z_factor, _log, intervals);

}

TimeSeriesD Grid::GetConcentrationBTCAtX(int species_counter, const double &x, const string &filename, const string &filename_d)
{
    TimeSeriesD output;
    for (int tt=0; tt<p[0][0].C.size(); tt++)
    {
        output.append(tt*TransportParameters.dt,GetConcentrationAtX(species_counter, x,tt));
    }
    if (filename!="")
        output.writefile(filename);
    if (filename_d!="")
        output.derivative().writefile(filename_d);
    return output;
}

double Grid::GetConcentrationAtX(int species_counter, const double &x, int timestep)
{
    int i=x/GeometricParameters.dx;
    double output = 0;
    for (int j=0; j<GeometricParameters.ny; j++)
        output += p[i][j].C[timestep][species_counter]/GeometricParameters.ny;

    return output;
}

TimeSeriesD Grid::GetConcentrationBTCAtX(const map<string,string> &Arguments)
{
    int species_id=0;
    if (Arguments.count("species")>0)
        species_id = aquiutils::atoi(Arguments.at("species"));
    if (Arguments.count("x")==0)
        return TimeSeriesD();
    string filename, filename_d;
    if (Arguments.count("filename")>0)
        filename = Arguments.at("filename");

    if (Arguments.count("filename_d")>0)
        filename_d = Arguments.at("filename_d");

    return GetConcentrationBTCAtX(species_id, aquiutils::atof(Arguments.at("x")),filename,filename_d);
}

CPathwaySet Grid::CreateTrajectories(const map<string,string> &Arguments)
{
    int n=1;
    double x_0=GeometricParameters.dx/2;
    if (Arguments.count("n")!=0)
        n = aquiutils::atoi(Arguments.at("n"));
    if (Arguments.count("x_0")!=0)
        x_0 = aquiutils::atof(Arguments.at("x_0"));

    vector<CPosition> pts = InitializeTrajectories(n,x_0);
}

vector<CPosition> Grid::InitializeTrajectories(int numpoints, const double &x_0)
{
    vector<CPosition> pts;
    pts.clear();
    int burnout = 0;
    TimeSeriesD boundary_v_dist = GetVelocityDistributionAtXSection(x_0,0);
    double v_max = boundary_v_dist.maxC();

    double y_0 = unitrandom()*GeometricParameters.dy*(GeometricParameters.ny-1);
    CPosition pt_0; pt_0.x = x_0; pt_0.y = y_0;
    double v_x = GetVelocity(pt_0)[0];
    pts.push_back(pt_0);
    for (int i = 1; i < numpoints; i++)
    {
        bool accepted = false;
        while (!accepted)
        {
            y_0 = unitrandom()*GeometricParameters.dy*(GeometricParameters.ny-1);
            pt_0.x = x_0; pt_0.y = y_0;
            v_x = GetVelocity(pt_0)[0];
            double u = unitrandom();
            if (u < (v_x / v_max/5)) accepted = true;
        }
        pts.push_back(pt_0);

        SetProgressValue(double(i) / double(numpoints));
    }

    return pts;

}

TimeSeriesD Grid::GetVelocityDistributionAtXSection(const double &x, int direction)
{
    TimeSeriesD out;
    int i = int(x / GeometricParameters.dx);
        for (int j = 0; j < GeometricParameters.ny; j++)
            out.append(i + j * 10000, p[i][j].V[direction]);
    return out;

}

CVector Grid::GetVelocity(const CPosition &pp)
{
    int i_floar_x = int(pp.x / GeometricParameters.dx);
    int j_floar_x = int(pp.y / GeometricParameters.dy+0.5);
    int i_floar_y = int(pp.x / GeometricParameters.dx+0.5);
    int j_floar_y = int(pp.y / GeometricParameters.dy);
    if (pp.x<=0 || pp.x>=(GeometricParameters.nx-1)*GeometricParameters.dx || pp.y<=0 || pp.y>=GeometricParameters.dy*(GeometricParameters.ny-1))
    {
        return CVector();
    }

    double vx1 = vx[i_floar_x][max(j_floar_x-1,0)] + 1.0/GeometricParameters.dx*(pp.x - GeometricParameters.dx*i_floar_x)*(vx[min(i_floar_x+1,GeometricParameters.nx-1)][max(j_floar_x-1,0)]-vx[i_floar_x][max(j_floar_x-1,0)]);
    double vx2 = vx[i_floar_x][min(j_floar_x,GeometricParameters.ny-2)] + 1.0/GeometricParameters.dx*(pp.x - GeometricParameters.dx*i_floar_x)*(vx[min(i_floar_x+1,GeometricParameters.nx-1)][min(j_floar_x,GeometricParameters.ny-2)]-vx[i_floar_x][min(j_floar_x,GeometricParameters.ny-2)]);
    double vx_interp = vx1 + (vx2-vx1)/GeometricParameters.dy*(pp.y-GeometricParameters.dy*(double(j_floar_x)-0.5));

    double vy1 = vy[max(i_floar_y-1,0)][j_floar_y] + 1.0/GeometricParameters.dy*(pp.y - GeometricParameters.dy*j_floar_y)*(vy[max(i_floar_y-1,0)][min(j_floar_y+1,GeometricParameters.nx-1)]-vy[max(i_floar_y-1,0)][j_floar_y]);
    double vy2 = vy[min(i_floar_y,GeometricParameters.nx-2)][j_floar_y] + 1.0/GeometricParameters.dy*(pp.y - GeometricParameters.dy*j_floar_y)*(vy[min(i_floar_y,GeometricParameters.nx-2)][min(j_floar_y+1,GeometricParameters.ny-1)]-vy[min(i_floar_y,GeometricParameters.nx-2)][j_floar_y]);
    double vy_interp = vy1 + (vy2-vy1)/GeometricParameters.dx*(pp.x-GeometricParameters.dx*(double(i_floar_y)-0.5));

    CVector V(2);
    V[0] = vx_interp;
    V[1] = vy_interp;
    return V;
}

CPathwaySet Grid::BuildTrajectories(const vector<CPosition> pts, const double &dx, const double &x_end, const double &tol, const double &diffusion)
{
    CPathwaySet X(int(pts.size()));
    unsigned int counter = 0;
    #pragma omp parallel for
    for (int i = 0; i < int(pts.size()); i++)
    {
        CPathway X1 = CreateSingleTrajectoryFixDx(pts[i], dx, x_end, diffusion, tol);
        {   X.weighted = false;
            X.paths[i] = X1;
        }
        #pragma omp critical
        {
            counter++;
            SetProgressValue(double(counter) / double(pts.size()));
        }
    }

    return X;
}

CPathway Grid::CreateSingleTrajectoryFixDx(const CPosition &pp, const double &dx0, const double &x_end, const double &D, const double &tol)
{
    double x0 = pp.x;
    double backward = 1;
    if (x_end<x0)
        backward = -1;

    CPosition pt = pp;
    CPathway Trajectory;
    Trajectory.weight = 1;
    double t = 0;

    bool ex = false;
    int counter = 0;
    while ( backward*x0<=backward*pt.x && backward*pt.x<=backward*x_end && ex==false)
    {
        counter ++;
        CVector V = GetVelocity(pt);
        if (V.getsize() == 0)
            {
                return Trajectory;
            }
        if (V[0]==0)
        {
            cout<< "Vx = 0!" <<endl;
        }

        double dx = fabs(dx0/(sqrt(pow(V[0],2)+pow(V[1],2)))*V[0]);

        bool changed_sign = true;
        while (changed_sign)
        {   if (V.num == 2)
            {
                if (V[0]==0)
                {
                    cout<<"Vx = 0!"<<endl;
                }

                double dt0 = fabs(dx/V[0]);
                if (dt0<0)
                {
                    cout<<"negative dt0!"<<endl;
                }

                CPosition p_new;
                p_new.x = pt.x + backward*dt0*V[0];//+diffusion[0];
                p_new.y = pt.y + backward*dt0*V[1];//+diffusion[1];
                CVector V_new = GetVelocity(p_new);
                if (V_new.getsize() == 0)
                {
                    return Trajectory;
                }
                if (V_new[0]*V[0]<0)
                {
                    cout<<"V changed sign!"<<endl;
                    dx/=2;
                }
                else
                {   changed_sign = false;
                    double err = norm(V-V_new);
                    double dt1;
                    while (err>=tol)
                    {

                        if (V[0]==V_new[0])
                            dt1 = fabs(dx/V[0]);
                        else
                            dt1 = fabs(dx*(log(fabs(V_new[0]))-log(fabs(V[0])))/(fabs(V_new[0])-fabs(V[0])));

                        p_new.y = pt.y + 0.5*(V[1]+V_new[1])*dt1*backward;// + diffusion[1]*sqrt(dt1/dt0);
                        p_new.x = pt.x + dx*backward; // + diffusion[0]*sqrt(dt1/dt0);
                        CVector V_star = GetVelocity(p_new);
                        if (V_star.getsize()!=2)
                            return Trajectory;
                        err = norm(V_star - V_new);
                        V_new = V_star;
                    }
                    p_new.v = V_new;
                    p_new.t = t+dt1;
                    CVector diffusion(2);

                    diffusion[0] = gsl_ran_gaussian(rng_ptr,sqrt(D*dt1));
                    diffusion[1] = gsl_ran_gaussian(rng_ptr,sqrt(D*dt1));
                    p_new.x = p_new.x + diffusion[0];
                    p_new.y = p_new.y + diffusion[1];

                    if (isnan(p_new.x))
                    {
                        cout<< "nan!" << endl;
                    }
                    Trajectory.append(p_new);
                    pt = p_new;

                    t += dt1;
                }
            }
            else ex = true;
        }
    }

    return Trajectory;
}




