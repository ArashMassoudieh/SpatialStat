#include "grid.h"
#include "Utilities.h"
#include "NormalDist.h"
#include "environment.h"
#include "Distribution.h"


vector<string> Grid::list_of_commands = vector<string>({"CreateGrid","AssignKField","WriteKFieldToVTP","RenormalizeKField","SolveHydro","WriteHydroSolutionToVTP"});

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
    if (cmd=="AssignKField")
        return AssignKFieldToGrid(arguments);
    if (cmd=="RenormalizeKField")
        return RenormalizeKField(arguments);
    if (cmd=="WriteKFieldToVTP")
        return WriteKFieldToVTP(arguments);
    if (cmd=="SolveHydro")
        return SolveHydro(arguments);
    if (cmd=="WriteHydroSolutionToVTP")
        return WriteHydroSolutionToVTP(arguments);
    return false;
}

bool Grid::AssignKFieldToGrid(map<string,string> Arguments)
{
    if (!parent->Object(Arguments["Distribution"]))
        return false;



    CDistribution* dist = dynamic_cast<CDistribution*>(parent->Object(Arguments["Distribution"]));
    Clear();

    field_gen_params Field_Generator_Parameters;

    if (Arguments.count("Maximum_neighboring_nodes")>0)
    {
        Field_Generator_Parameters.max_correl_n = aquiutils::atof(Arguments.at("Maximum_neighboring_nodes"));
    }
    Field_Generator_Parameters.inversecdf = dist->inverse_cumulative;
    Field_Generator_Parameters.k_correlation_lenght_scale_x = aquiutils::atof(Arguments["correlation_length_x"]);
    Field_Generator_Parameters.k_correlation_lenght_scale_y = aquiutils::atof(Arguments["correlation_length_y"]);
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

bool Grid::RenormalizeKField(map<string,string> Arguments)
{
    CDistribution* dist = dynamic_cast<CDistribution*>(parent->Object(Arguments["Distribution"]));
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

void Grid::SetProgressValue(const double &s)
{
#ifdef QT_version
    main_window->get_ui()->progressBar->setValue(s*100);
    QApplication::processEvents();
#endif // QT_version
    cout << "\r Progress: " << s*100 << "%                                     ";
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


