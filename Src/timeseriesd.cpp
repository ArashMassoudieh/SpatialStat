#include "timeseriesd.h"
#include "command.h"


TimeSeriesD::TimeSeriesD():CTimeSeries<double>()
{

}

TimeSeriesD::TimeSeriesD(int n):CTimeSeries<double>(n)
{

}

TimeSeriesD::TimeSeriesD(const CTimeSeries<double> &timeseries):CTimeSeries<double>(timeseries)
{

}

FunctionOutPut TimeSeriesD::Execute(const string &cmd, const map<string,string> &arguments)
{
    FunctionOutPut output;
    if (cmd=="WriteTimeSeriesToFile")
    {   output.success = WriteToFile(arguments);

    }
    return output;
}

bool TimeSeriesD::WriteToFile(const map<string,string> &arguments)
{
    string filename;
    if (arguments.count("filename")==0)
        return false;
    else
        filename = arguments.at("filename");
    int interval = 1;
    if (arguments.count("interval")!=0)
        interval = aquiutils::atoi(arguments.at("interval"));

    return writefile(filename);

}

bool TimeSeriesD::HasCommand(const string &cmd)
{
    if (aquiutils::lookup(Commands(),cmd)!=-1)
        return true;
    else
        return false;
}

vector<string> TimeSeriesD::Commands()
{
    vector<string> cmds;
    for (map<string,command_parameters>::iterator i=Command::Command_Structures.begin(); i!=Command::Command_Structures.end(); i++)
    {
        if (i->second.Object==object_type::timeseries)
            cmds.push_back(i->first);
    }
    return cmds;
}

