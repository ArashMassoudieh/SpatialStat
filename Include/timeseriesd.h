#ifndef TIMESERIESD_H
#define TIMESERIESD_H

#include "BTC.h"
#include "interface.h"

class TimeSeriesD : public CTimeSeries<double>, public Interface
{
public:
    TimeSeriesD();
    TimeSeriesD(int n);
    FunctionOutPut Execute(const string &cmd, const map<string,string> &arguments);
    static bool HasCommand(const string &cmd);
    static vector<string> Commands();
};

#endif // TIMESERIESD_H
