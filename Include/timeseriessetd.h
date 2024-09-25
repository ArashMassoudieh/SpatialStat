#ifndef TIMESERIESSETD_H
#define TIMESERIESSETD_H

#include "BTCSet.h"
#include "interface.h"

class CDistribution;

class TimeSeriesSetD : public CTimeSeriesSet<double>, public Interface
{
public:
    TimeSeriesSetD();
    TimeSeriesSetD(int n);
    TimeSeriesSetD(const CTimeSeriesSet<double> &timeseries);
    FunctionOutPut Execute(const string &cmd, const map<string,string> &arguments);
    static bool HasCommand(const string &cmd);
    static vector<string> Commands();
private:
    bool WriteToFile(const map<string,string> &arguments);
};

#endif // TIMESERIESSETD_H
