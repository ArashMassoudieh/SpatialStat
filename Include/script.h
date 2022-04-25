#ifndef SCRIPT_H
#define SCRIPT_H
#include "command.h"
#include "vector"


class Script
{
public:
    Script();
    Script(const string &filename);
    bool GetFromFile(const string &filename);
    std::vector<Command> Commands;
};

#endif // SCRIPT_H
