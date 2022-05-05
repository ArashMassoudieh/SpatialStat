#ifndef INTERFACE_H
#define INTERFACE_H

#include <string>
#include <vector>
#include <map>

enum class object_type {PathwaySet, D2Grid, Distribution};

class Environment;

using namespace std;

class Interface
{
public:
    Environment *parent;
    Interface(Environment* _parent = nullptr);
    object_type ObjectType;
    virtual vector<string> commands();
    bool HasCommand(const string &cmd);
    virtual bool Execute(const string &cmd, const map<string,string> &arguments)=0;

};

#endif // INTERFACE_H
