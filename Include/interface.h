#pragma once
#ifndef INTERFACE_H
#define INTERFACE_H

#include <string>
#include <vector>
#include <map>
#include "Structs.h"


class Environment;
class QTableWidget;


using namespace std;

class Interface
{
public:
    Environment *parent;
    Interface(Environment* _parent = nullptr);
    object_type ObjectType;
    virtual vector<string> commands();
    bool HasCommand(const string &cmd);
    virtual FunctionOutPut Execute(const string &cmd, const map<string,string> &arguments)=0;
    QTableWidget *outputwindow();
    void SetProgressValue(const double &x);


};

#endif // INTERFACE_H
