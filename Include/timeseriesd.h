#ifndef TIMESERIESD_H
#define TIMESERIESD_H

#include "BTC.h"
#include "interface.h"

class TimeSeriesD : public CTimeSeries<double>, public Interface
{
public:
    TimeSeriesD();
    TimeSeriesD(int n);
    TimeSeriesD(const CTimeSeries<double> &timeseries);
    FunctionOutPut Execute(const string &cmd, const map<string,string> &arguments);
    static bool HasCommand(const string &cmd);
    vector<string> commands();
    static vector<string> Commands();
    CDistribution GetDistribution(int nbins);
private:
    bool WriteToFile(const map<string,string> &arguments);
};

#endif // TIMESERIESD_H
