#ifndef INTERFACE_H
#define INTERFACE_H

#include <string>
#include <vector>

enum class object_type {PathwaySet, D2Grid, Distribution};

using namespace std;

class Interface
{
public:
    Interface();
    object_type ObjectType;
    virtual vector<string> commands();
    bool HasCommand(const string &cmd);

};

#endif // INTERFACE_H
