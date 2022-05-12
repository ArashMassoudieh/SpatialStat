#include "timeseriesd.h"

vector<string> TimeSeriesD::list_of_commands = vector<string>({"CreateDistribution","WriteToFile", "SetInverseCumulative", "WriteInverseCumulativeToFile"});

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
    //return vector<string>();
    return list_of_commands;
}
