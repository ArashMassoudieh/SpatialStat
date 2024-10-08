#include "PathwaySet.h"
#include "Utilities.h"
#include "grid.h"
#include "command.h"


CPathwaySet::CPathwaySet()
{
}

CPathwaySet::CPathwaySet(unsigned int n)
{
    paths.resize(n);
}


CPathwaySet::CPathwaySet(const CPathwaySet & P)
{
	paths = P.paths;
    weighted = P.weighted;

}

CPathwaySet::CPathwaySet(int number_of_paths)
{
    paths.resize(number_of_paths);
}

CPathwaySet &CPathwaySet::operator=(const CPathwaySet & P)
{
	paths = P.paths;
	return *this;
}


CPathwaySet::~CPathwaySet()
{
}

void CPathwaySet::write(string filename, int interval)
{
	ofstream file;
	file.open(filename.c_str());
	for (int j = 0; j < paths.size(); j+=interval)
	{
        file << "t_" << aquiutils::numbertostring(j) << ",x_" << aquiutils::numbertostring(j) << ",y_" << aquiutils::numbertostring(j) << ",u_" << aquiutils::numbertostring(j) << ",v_" << aquiutils::numbertostring(j) << ",u_" << aquiutils::numbertostring(j) << ",z_" << aquiutils::numbertostring(j) << ",";
	}
	file << endl;
	for (int i = 0; i < max_num_points(); i++)
	{
		for (int j = 0; j < paths.size(); j+=interval)
		{
			if (i < paths[j].positions.size())
				file << paths[j].positions[i].t << "," << paths[j].positions[i].x << "," << paths[j].positions[i].y << "," << paths[j].positions[i].v[0] << "," << paths[j].positions[i].v[1] << "," << paths[j].positions[i].u << "," << paths[j].positions[i].z << ",";
			else
				file << "," << "," << "," << "," << "," << "," << ",";
		}
		file << endl;
	}

file.close();
}

void CPathwaySet::append(const CPathway & P, double weight)
{
   paths.push_back(P);
}

int CPathwaySet::max_num_points()
{
	int max_np = 0;
	for (int i = 0; i < paths.size(); i++)
		max_np = max(max_np, int(paths[i].positions.size()));

	return max_np;
}

void CPathwaySet::create_ou_paths(int n, CDistribution * dist, double x_min, double x_max, double kappa, double dx, double weight)
{
	for (int i = 0; i < n; i++)
	{
            CPathway P;
            P.create_ou(dist, x_min, x_max, kappa, dx);
            P.weight = weight;
            append(P);
	}
}

void CPathwaySet::create_copula_paths(int n, CDistribution * dist, double x_min, double x_max, double epsilon, double r, double dx, double weight)
{
	for (int i = 0; i < n; i++)
	{
            CPathway P;
            P.create_copula(dist, x_min, x_max, epsilon, r, dx);
            P.weight = weight;
            append(P);
	}
}

void CPathwaySet::create_copula_paths(int n, CDistribution * dist, double x_min, double x_max, double epsilon, CCopula *copula, double dx, double weight)
{
	for (int i = 0; i < n; i++)
	{
            CPathway P;
            P.create_copula(dist, x_min, x_max, epsilon, copula, dx);
            P.weight = weight;
            append(P);
	}
}


bool CPathwaySet::write_vtk(vtkSmartPointer<vtkPolyDataMapper> mapper, string filename)
{
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(mapper->GetInput());
	// This is set so we can see the data in a text editor.
	writer->SetDataModeToAscii();
	writer->Write();
    return true;

}

vtkSmartPointer<vtkPolyDataMapper> CPathwaySet::pathways_vtk_pdt_vtp(double z_factor, double offset)
{
	vector<vtkSmartPointer<vtkPolyData>> outarray;
	for (int i = 0; i < paths.size(); i++)
		outarray.push_back(paths[i].pathway_vtk_pdt_vtp(z_factor, offset));

	vtkSmartPointer<vtkAppendPolyData> appendFilter =
		vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	appendFilter->AddInputConnection(input1->GetProducerPort());
	appendFilter->AddInputConnection(input2->GetProducerPort());
#else
	for (int i = 0; i < outarray.size(); i++)
		appendFilter->AddInputData(outarray[i]);
#endif
	appendFilter->Update();

	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputConnection(appendFilter->GetOutputPort());
#endif

	return mapper;
}

CPathway CPathwaySet::snapshotattime(double t)
{
    CPathway Ptwy;
    for (int i = 0; i < paths.size(); i++)
	Ptwy.append(paths[i].get_position_at_t(t));

    return Ptwy;
}

CPathway CPathwaySet::snapshotatlocation(double x)
{
	CPathway Ptwy;
	for (int i = 0; i < paths.size(); i++)
            Ptwy.append(paths[i].get_position_at_x(x));

	return Ptwy;
}

bool CPathwaySet::make_uniform_at_x(double dx)
{
	for (int i = 0; i < paths.size(); i++)
		paths[i] = paths[i].make_uniform_x(dx);
    return true;
}

CTimeSeries<double> CPathwaySet::sample_velocities()
{
	cout<<"Getting velocities from the trajectories"<<endl;
    CTimeSeries<double> out;
	for (int i = 0; i < paths.size(); i++)
    {
       for (int j = 0; j < paths[i].size(); j++)
            out.append(paths[i].positions[j].v[0]);
    set_progress_value(i);
    }
    return out;
}

void CPathwaySet::make_uniform_at_t(double dt)
{
	for (int i = 0; i < paths.size(); i++)
		paths[i] = paths[i].make_uniform_t(dt);

}



CPosition CPathwaySet::get_pair_v_pos(int increment, int num_seq)
{
    CPosition p(num_seq);
    int i = int(unitrandom()*paths.size());
    int j = int((paths[i].positions.size()- (num_seq-1)*increment)*unitrandom());
    for (int ii=0 ; ii<num_seq; ii++)
       p.v[ii] = paths[i].positions[j+ii*increment].v[0];

    p.weight = paths[i].weight;
    return p;
}

CTimeSeriesSet<double> CPathwaySet::get_pair_v(int increment, int n, int num_seq)
{
    CTimeSeriesSet<double> out(num_seq);
	for (int i = 0; i < n; i++)
	{
		CPosition p = get_pair_v_pos(increment, num_seq);
		if (weighted)
            out.append(double(i), p.v.vec); //used to have weight
        else
            out.append(double(i), p.v.vec);
		cout << "\r" << float(i)/float(n)*100 << "%" << std::flush;
	}
	return out;

}

CTimeSeries<double> CPathwaySet::get_BTC(double x, int n_bins, bool velweight, double smoothing_factor)
{
    CTimeSeries<double> BTC;

    for (int i = 0; i < paths.size(); i++)
        BTC.append(i, paths[i].get_cross_time(x));

    return BTC.distribution(n_bins,(BTC.maxC()-BTC.minC())*smoothing_factor, 0);

}

CTimeSeries<double> CPathwaySet::get_BTC_points(double x, bool vel_inv_weighted)
{
    CTimeSeries<double> BTC;

    for (int i = 0; i < paths.size(); i++)
        BTC.append(i, paths[i].get_cross_time(x));

    return BTC;

}

bool CPathwaySet::AssignVelocities()
{
    show_in_window("Calculating velocities for trajectories");
    for (int i=0; i<paths.size(); i++)
    {
        set_progress_value((i*100)/paths.size());
        paths[i].positions[0].v[0] = (paths[i].positions[1].x-paths[i].positions[0].x)/(paths[i].positions[1].t-paths[i].positions[0].t);
        paths[i].positions[0].v[1] = (paths[i].positions[1].y-paths[i].positions[0].y)/(paths[i].positions[1].t-paths[i].positions[0].t);
        paths[i].positions[0].v[2] = (paths[i].positions[1].z-paths[i].positions[0].z)/(paths[i].positions[1].t-paths[i].positions[0].t);

        for (int j=1; j<paths[i].positions.size()-1; j++)
        {
            paths[i].positions[j].v[0] = (paths[i].positions[j+1].x-paths[i].positions[j-1].x)/(paths[i].positions[j+1].t-paths[i].positions[j-1].t);
            paths[i].positions[j].v[1] = (paths[i].positions[j+1].y-paths[i].positions[j-1].y)/(paths[i].positions[j+1].t-paths[i].positions[j-1].t);
            paths[i].positions[j].v[2] = (paths[i].positions[j+1].z-paths[i].positions[j-1].z)/(paths[i].positions[j+1].t-paths[i].positions[j-1].t);

            if (paths[i].positions[j].v[0]==0)
            {
                show_in_window("Velocity zero for path: " + aquiutils::numbertostring(i) + "@ position: " + aquiutils::numbertostring(j));
                paths[i].positions[j].v[0] = 1e-6;
            }
        }

        paths[i].positions[paths[i].positions.size()-1].v[0] = (paths[i].positions[paths[i].positions.size()-1].x-paths[i].positions[paths[i].positions.size()-2].x)/(paths[i].positions[paths[i].positions.size()-1].t-paths[i].positions[paths[i].positions.size()-2].t);
        paths[i].positions[paths[i].positions.size()-1].v[1] = (paths[i].positions[paths[i].positions.size()-1].y-paths[i].positions[paths[i].positions.size()-2].y)/(paths[i].positions[paths[i].positions.size()-1].t-paths[i].positions[paths[i].positions.size()-2].t);
        paths[i].positions[paths[i].positions.size()-1].v[2] = (paths[i].positions[paths[i].positions.size()-1].z-paths[i].positions[paths[i].positions.size()-2].z)/(paths[i].positions[paths[i].positions.size()-1].t-paths[i].positions[paths[i].positions.size()-2].t);
        set_progress_value(double(i)/double(paths.size()) );
    }
    show_in_window("Calculating velocities for trajectories, done!");
    return true;
}

void CPathwaySet::show_in_window(string s)
{
    #ifdef QT_version
    qDebug()<<QString::fromStdString(s);
    main_window->get_ui()->ShowOutput->append(QString::fromStdString(s));
    QApplication::processEvents();
    #else
    cout<<s<<endl;
    #endif // Qt_version
}

bool CPathwaySet::getfromMODflowfile(const string &filename)
{
    ifstream file;
    file.open (filename, std::fstream::in);
    if (!file.good()) return false;
    int rownum = 0;
    double age=0;
    while (!file.eof())
    {
        vector<double> s1 = aquiutils::ATOF(aquiutils::getline(file,' '));
        if (n() == 0)
            if (s1[0]!=0)
            {   paths.resize(s1[0]);
                cout<<"Number of paths = " << s1[0]<< endl;

            }
        if (s1.size()==2)
        {
            int numberofpoints = s1[0];

            age = s1[1];

            for (int i=0; i<numberofpoints; i++)
            {
                //cout <<i<< endl;
                vector<double> s2 = aquiutils::ATOF(aquiutils::getline(file,' '));

                if (s2.size()>5);
                {
                    //cout <<i<<","<<s2[0]<<","<<s2[1]<<","<<s2[2]<<","<<s2[3]<<","<<s2[4]<<","<<s2[5]<< endl;
                    CPosition P;
                    P.x = s2[1];
                    P.y = s2[2];
                    P.z = s2[3];
                    P.v = CVector(3);
                    P.t = s2[4];

                    P.weight = s2[5];
                    paths[(int)s2[0]-1].append(P);
                    //cout <<"Done!"<<endl;
                }
            }
            //cout << "Number of points = " << numberofpoints << "," << endl;
            set_progress_value(aquiutils::numbertostring(age)+", number of points:" + aquiutils::numbertostring(numberofpoints));
        }
        else
        {
            cout<<"!"<<s1.size();
        }

        rownum++;

    }
    AssignVelocities();
    file.close();
    //cout<<"Reading Trajectories Done!"<<endl;
    return true;
}


bool CPathwaySet::getfromShermanfile(const string &filename, const string &reactionfile, int columnnumber)
{
    ifstream file, filereaction;
    file.open (filename, std::fstream::in);
    if (!file.good())
    {
        show_in_window("The program was not able to open " + filename);
        return false;
    }
    if (reactionfile!="")
    {
        filereaction.open(reactionfile, std::fstream::in);
        if (!filereaction.good())
        {
            show_in_window("The program was not able to open " + reactionfile);
            return false;
        }
    }

    int rownum = 0;
    double age=0;
    int maxparticlecount=0;
    while (!file.eof())
    {
        vector<double> s1 = aquiutils::ATOF(aquiutils::getline(file,' '));
        if (s1.size()>4)
            maxparticlecount = max(maxparticlecount,int(s1[0]));

    }
    show_in_window("maximum number of trajectories:" + aquiutils::numbertostring(maxparticlecount));
    paths.resize(maxparticlecount);
    file.close();
    file.open (filename, std::fstream::in);

    while (!file.eof())
    {
        vector<double> s1 = aquiutils::ATOF(aquiutils::getline(file,' '));
        vector<int> s_react;
        if (reactionfile!="")
            s_react = aquiutils::ATOI(aquiutils::getline(filereaction, ' '));
        if (s1.size()>5)
        {
            set_progress_value(s1[1]);
            //cout <<s1[0]<<","<<s1[1]<<","<<s1[2]<<","<<s1[3]<<","<<s1[4]<< endl;
            CPosition P;
            P.x = s1[2];
            P.y = s1[3];
            P.z = s1[4];
            P.v = CVector(3);
            P.t = s1[1];
            if (reactionfile!="")
            {
                P.reacted = s_react[columnnumber];
            }
            if (int(s1[0])>0 && int(s1[0])<=maxparticlecount)
                paths[(int)s1[0]-1].append(P);
            //cout <<"Done!"<<endl;
        }

    }
    show_in_window("Reading trajectories done!");
    AssignVelocities();
    file.close();
    //cout<<"Reading Trajectories Done!"<<endl;
    return true;
}

bool CPathwaySet::getfromShermanfile_v(const string &filename)
{
    ifstream file;
    file.open (filename, std::fstream::in);
    if (!file.good()) return false;
    int rownum = 0;
    double age=0;
    int maxparticlecount=0;
    while (!file.eof())
    {
        vector<double> s1 = aquiutils::ATOF(aquiutils::getline(file,' '));
        if (s1.size()>4)
            maxparticlecount = max(maxparticlecount,int(s1[0]));

    }
    show_in_window("maximum number of trajectories:" + aquiutils::numbertostring(maxparticlecount));
    paths.resize(maxparticlecount);
    file.close();
    file.open (filename, std::fstream::in);

    while (!file.eof())
    {
        vector<double> s1 = aquiutils::ATOF(aquiutils::getline(file,' '));
        if (s1.size()>5)
        {
            set_progress_value(s1[1]);
            //cout <<s1[0]<<","<<s1[1]<<","<<s1[2]<<","<<s1[3]<<","<<s1[4]<< endl;
            CPosition P;
            P.x = s1[2];
            P.y = s1[3];
            P.z = s1[4];
            P.v = CVector(3);
            P.v[0] = s1[6];
            P.v[1] = s1[7];
            P.v[2] = s1[8];
            P.t = s1[1];

            if (int(s1[0])>0 && int(s1[0])<=maxparticlecount)
                paths[(int)s1[0]-1].append(P);
            //cout <<"Done!"<<endl;
        }

    }
    show_in_window("Reading trajectories done!");
    //AssignVelocities();
    file.close();
    //cout<<"Reading Trajectories Done!"<<endl;
    return true;
}


void CPathwaySet::set_progress_value(double s)
{
#ifdef QT_version
	main_window->get_ui()->progressBar->setValue(s*100);
	QApplication::processEvents();
#endif // QT_version
    cout << "\r Progress: " << s << "                         ";
}

void CPathwaySet::set_progress_value(string s)
{
#ifdef QT_version
	main_window->get_ui()->progressBar->setValue(s*100);
	QApplication::processEvents();
#endif // QT_version
    cout << "\r Progress: " << s << "                         ";
}

FunctionOutPut CPathwaySet::Execute(const string &cmd, const map<string,string> &arguments)
{
    FunctionOutPut output;
    if (cmd=="Uniformize")
    {   output.success = Uniformize(arguments);
        output.output = nullptr;
    }
    if (cmd=="WritePathwayToVTP")
    {
        output.success = WriteToVTP(arguments);
        output.output = nullptr;
    }
    return output;
}

bool CPathwaySet::HasCommand(const string &cmd)
{
    if (aquiutils::lookup(Commands(),cmd)!=-1)
        return true;
    else
        return false;
}

vector<string> CPathwaySet::Commands()
{
    vector<string> cmds;
    for (map<string,command_parameters>::iterator i=Command::Command_Structures.begin(); i!=Command::Command_Structures.end(); i++)
    {
        if (i->second.Object==object_type::pathwayset)
            cmds.push_back(i->first);
    }
    return cmds;
}

bool CPathwaySet::Uniformize(const map<string,string> &Arguments)
{
    if (Arguments.count("dx")==0)
        return false;

    double dx = aquiutils::atof(Arguments.at("dx"));

    return make_uniform_at_x(dx);
}

bool CPathwaySet::WriteToVTP(const map<string,string> &Arguments)
{

    string filename;
    if (Arguments.count("filename")==0)
        return false;
    else
        filename = Arguments.at("filename");

    double z_factor=0;
    double offset=0;
    bool _log=false;
    bool _color=false;
    int interval=1;
    if (Arguments.count("z_factor")!=0)
        z_factor = aquiutils::atof(Arguments.at("z_factor"));

    if (Arguments.count("offset")!=0)
        offset = aquiutils::atof(Arguments.at("offset"));

    if (Arguments.count("log")!=0)
        _log = aquiutils::atoi(Arguments.at("log"));

    if (Arguments.count("color")!=0)
        _color = aquiutils::atoi(Arguments.at("color"));

    if (Arguments.count("interval")!=0)
        interval = aquiutils::atoi(Arguments.at("interval"));

    WriteToVTP(filename,z_factor,offset,_log,_color, interval);

    return true;
}

void CPathwaySet::WriteToVTP(const string &filename, const double &z_factor, const double &offset, bool _log, bool _color, int interval)
{
    vector<vtkSmartPointer<vtkPolyData>> outputmappers;
        double max_v_x = 1;
        set_progress_value(0);
        for (int i = 0; i < paths.size(); i+=interval)
        {	outputmappers.push_back(paths[i].TovtkPolyData(z_factor, offset, _log, _color));
            set_progress_value((double)i/(double)n());
        }

        vtkSmartPointer<vtkAppendPolyData> appendFilter =
            vtkSmartPointer<vtkAppendPolyData>::New();
    #if VTK_MAJOR_VERSION <= 5
        appendFilter->AddInputConnection(input1->GetProducerPort());
        appendFilter->AddInputConnection(input2->GetProducerPort());
    #else
        for (int i=0; i<outputmappers.size(); i++)
            appendFilter->AddInputData(outputmappers[i]);
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

        #ifdef QT_version
        main_window->get_ui()->ShowOutput->append("Writing vtp file... ");
        #endif // QT_version
        vtkSmartPointer<vtkXMLPolyDataWriter> writer =
            vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName(filename.c_str());
        writer->SetInputData(mapper->GetInput());
        // This is set so we can see the data in a text editor.
        writer->SetDataModeToAscii();
        writer->Write();

}



