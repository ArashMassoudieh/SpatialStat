#include "timeseriesd.h"
#include "command.h"


TimeSeriesD::TimeSeriesD():CTimeSeries<double>()
{

}

TimeSeriesD::TimeSeriesD(int n):CTimeSeries<double>(n)
{

}

FunctionOutPut TimeSeriesD::Execute(const string &cmd, const map<string,string> &arguments)
{
    FunctionOutPut output;
    return output;
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
