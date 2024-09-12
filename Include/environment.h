#pragma once
#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <vector>
#include <map>
#include "command.h"
#include <QTableWidget>
#include <QProgressBar>
using namespace std;

class Interface;

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
    QProgressBar *current_progress_bar;

};

#endif // ENVIRONMENT_H
