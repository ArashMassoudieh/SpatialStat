#ifndef TIMESERIESD_H
#define TIMESERIESD_H

#include "BTC.h"
#include "interface.h"

class TimeSeriesD : public CTimeSeries<double>, public Interface
{
public:
    TimeSeriesD();
    FunctionOutPut Execute(const string &cmd, const map<string,string> &arguments);
};

#endif // TIMESERIESD_H
