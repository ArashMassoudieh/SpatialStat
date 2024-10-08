#include "timeseriessetd.h"
#include "timeseriesd.h"
#include "command.h"
#include "Distribution.h"


TimeSeriesSetD::TimeSeriesSetD():CTimeSeriesSet<double>()
{

}

TimeSeriesSetD::TimeSeriesSetD(int n):CTimeSeriesSet<double>(n)
{

}

TimeSeriesSetD::TimeSeriesSetD(const CTimeSeriesSet<double> &timeseries):CTimeSeriesSet<double>(timeseries)
{

}

FunctionOutPut TimeSeriesSetD::Execute(const string &cmd, const map<string,string> &arguments)
{
    FunctionOutPut output;
    if (cmd=="WriteTimeSeriesToFile")
    {   output.success = WriteToFile(arguments);

    }
    return output;
}

bool TimeSeriesSetD::WriteToFile(const map<string,string> &arguments)
{
    string filename;
    if (arguments.count("filename")==0)
        return false;
    else
        filename = arguments.at("filename");
    int interval = 1;
    if (arguments.count("interval")!=0)
        interval = aquiutils::atoi(arguments.at("interval"));

    return writetofile(filename);

}

bool TimeSeriesSetD::HasCommand(const string &cmd)
{
    if (aquiutils::lookup(Commands(),cmd)!=-1)
        return true;
    else
        return false;
}

vector<string> TimeSeriesSetD::Commands()
{
    vector<string> cmds;
    for (map<string,command_parameters>::iterator i=Command::Command_Structures.begin(); i!=Command::Command_Structures.end(); i++)
    {
        if (i->second.Object==object_type::timeseries)
            cmds.push_back(i->first);
    }
    return cmds;
}

