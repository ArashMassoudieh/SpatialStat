#include "MapAsTimeSeriesSet.h"
#include "Utilities.h""

MapAsTimeSeriesSet::MapAsTimeSeriesSet():CTimeSeriesSet<double>()
{
    //ctor
}

MapAsTimeSeriesSet::~MapAsTimeSeriesSet()
{
    //dtor
}

MapAsTimeSeriesSet::MapAsTimeSeriesSet(const MapAsTimeSeriesSet& other):CTimeSeriesSet<double>(other)
{
    x = other.x;
}

MapAsTimeSeriesSet& MapAsTimeSeriesSet::operator=(const MapAsTimeSeriesSet& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    CTimeSeriesSet<double>::operator=(rhs);
    x = rhs.x;
    //assignment operator

    return *this;
}

void MapAsTimeSeriesSet::append(CTimeSeries<double>& _BTC, double xx)
{
    CTimeSeriesSet<double>::append(_BTC, aquiutils::numbertostring(xx));
    x.push_back(xx);
}

bool MapAsTimeSeriesSet::writetofile(const string &filename)
{
    CTimeSeriesSet<double>::writetofile(filename,1);
    return true;
}
