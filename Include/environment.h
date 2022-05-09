#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <vector>
#include <map>
#include <interface.h>
#include "command.h"
#include <QTableWidget>

using namespace std;

class Environment
{
public:
    QTableWidget *outputwindow = nullptr;
    Environment(QTableWidget *outputwindow = nullptr);
    map<string, Interface*> Objects;
    Interface* Object(const string &objectname)
    {
        if (Objects.count(objectname)>0)
            return Objects.at(objectname);
        else
            return nullptr;

    }
    bool Execute(const Command &cmd);

};

#endif // ENVIRONMENT_H
